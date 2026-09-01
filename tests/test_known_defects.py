"""
Characterization tests for known defects in bubble detection and deduplication.

Each test below asserts the behaviour PICOTA *should* have. Those that are
currently wrong are marked `xfail(strict=True)`, so the moment a fix lands the
test XPASSes, strict mode turns that into a failure, and whoever made the fix is
forced to remove the marker. This file is therefore both the defect record and
the regression net for roadmap phases 1-4 (see docs/ROADMAP.md).

Defect IDs match the table in docs/ROADMAP.md.

Run:
    python3 -m pytest tests/test_known_defects.py -v
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "picota"))
sys.path.insert(0, os.path.dirname(__file__))

from src.cycle_finderv2 import (  # noqa: E402
    Cycle,
    GraphWork,
    cycle_info_optimized,
    cycle_match_based_on_contig_id,
    reverse_complement,
)
from src.cycle_kmer_hash import filter_cycles_with_kmer, get_kmer_hashes  # noqa: E402
from synthetic_gfa import (  # noqa: E402
    make_dna,
    mutate,
    shared_repeat_n_cargos,
    shared_repeat_two_cargos,
)

K_MER_SIM = 80
THRESHOLD_SIM = 80


def run_cycle_pipeline(gfa_path, k_mer_sim=K_MER_SIM, threshold_sim=THRESHOLD_SIM,
                       path_limit=25):
    """Run detection + both deduplication stages, mirroring cycle_analysis()."""
    graph_work = GraphWork()
    graph_work.find_all_path = True
    graph_work.path_limit = path_limit

    node_dict, edge_dict = graph_work.parse_gfa(gfa_path)
    graph_work.dfs_iterative(graph_work.generate_genome_graph(node_dict, edge_dict))

    node_lengths = {k: len(v["Sequence"]) for k, v in node_dict.items()}
    unique_paths = []
    for path in graph_work.allPaths:
        if cycle_match_based_on_contig_id(path, node_lengths, unique_paths):
            unique_paths.append(path)

    cycles = []
    for path in unique_paths:
        cycle = cycle_info_optimized(path, node_dict, edge_dict, cycles)
        if cycle in (None, "Pass"):
            continue
        cycle.name = "tmp"
        cycles.append(cycle)

    return filter_cycles_with_kmer(cycles, k_mer_sim, threshold_sim, "Cycle")


# ─── D1: node depth is discarded by parse_gfa ────────────────────────────────

def test_d1_parse_gfa_retains_node_depth(tmp_path):
    """
    Coverage is the only signal that says how many genomic copies an IS node
    represents, and assemblers hand it to us for free (SPAdes encodes it in the
    node name, MEGAHIT/gfatools as KC:i:, the GFA spec as dp:f:). parse_gfa
    used to drop every tag; resolved in roadmap phase 1.
    """
    gfa = shared_repeat_two_cargos(str(tmp_path / "shared.gfa"),
                                   repeat_depth=60.0, cargo_depth=30.0)
    node_dict, _ = GraphWork().parse_gfa(gfa)
    assert "Depth" in node_dict["repeat+"]
    assert node_dict["repeat+"]["Depth"] == pytest.approx(60.0, rel=0.05)
    assert node_dict["cargoA+"]["Depth"] == pytest.approx(30.0, rel=0.05)


def test_d1_depth_tags_are_present_in_fixture(tmp_path):
    """Guard: the fixture really does emit depth tags, so D1 is about the parser."""
    gfa = shared_repeat_two_cargos(str(tmp_path / "shared.gfa"))
    seg_lines = [ln for ln in open(gfa) if ln.startswith("S")]
    assert all("dp:f:" in ln and "KC:i:" in ln for ln in seg_lines)


# ─── D2: path dedup deletes bubbles that share a repeat node ─────────────────

def recovered_cargo_count(gfa_path, node_dict, n_cargos):
    """How many of the ground-truth cargos appear in the reported cycles."""
    cycles = run_cycle_pipeline(gfa_path)
    found = set()
    for cycle in cycles:
        both_strands = cycle.sequence + "|" + reverse_complement(cycle.sequence)
        for i in range(n_cargos):
            if node_dict[f"cargo{i}+"]["Sequence"] in both_strands:
                found.add(i)
    return len(found), len(cycles)


@pytest.mark.parametrize("n_cargos", [3, 4, 5])
@pytest.mark.xfail(strict=True, reason="D2: output saturates at 2 regardless of ground truth")
def test_d2_every_cargo_sharing_one_repeat_is_reported(tmp_path, n_cargos):
    """
    N distinct composite transposons sharing one IS. Ground truth is N.

    Measured today: 2/2, 2/3, 2/4, 2/5. Once the first plus-strand bubble is
    accepted, every remaining plus-strand bubble exceeds the 70% shared-node
    threshold and is dropped; exactly one minus-strand bubble survives because
    the similarity check compares orientation-specific node names while the
    duplicate check strips orientation. The result is capped at two candidates
    no matter how many real CTs are present.
    """
    gfa = shared_repeat_n_cargos(str(tmp_path / "n.gfa"), n_cargos)
    node_dict, _ = GraphWork().parse_gfa(gfa)
    recovered, _ = recovered_cargo_count(gfa, node_dict, n_cargos)
    assert recovered == n_cargos


def test_d2_two_cargos_are_recovered_by_strand_accident(tmp_path):
    """
    The n=2 case does succeed -- but only because the second CT survives as its
    minus-strand traversal. Pinned so a future fix cannot silently regress it.
    """
    gfa = shared_repeat_n_cargos(str(tmp_path / "n2.gfa"), 2)
    node_dict, _ = GraphWork().parse_gfa(gfa)
    recovered, reported = recovered_cargo_count(gfa, node_dict, 2)
    assert recovered == 2
    assert reported == 2


def test_d2_survives_when_repeat_is_short_relative_to_cargo(tmp_path):
    """The same topology is handled correctly when the shared repeat is smaller."""
    gfa = shared_repeat_two_cargos(str(tmp_path / "shared.gfa"),
                                   repeat_len=820, cargo_a_len=1500, cargo_b_len=1700)
    assert len(run_cycle_pipeline(gfa)) == 2


@pytest.mark.parametrize("repeat_len,cargo_len,expected_kept", [
    (2500, 1400, True),   # 64.1% shared -> below threshold
    (2500, 1000, False),  # 71.4% shared -> dropped
    (2500, 900, False),   # 73.5% shared -> dropped
])
def test_d2_threshold_boundary_is_length_driven(repeat_len, cargo_len, expected_kept):
    """
    Documents that whether a real CT survives depends only on the repeat/cargo
    length ratio -- not on any biological property of the candidate.
    """
    node_lengths = {"repeat+": repeat_len, "cargoA+": cargo_len, "cargoB+": cargo_len}
    accepted = [["repeat+", "cargoA+"]]
    kept = cycle_match_based_on_contig_id(["repeat+", "cargoB+"], node_lengths, accepted)
    assert kept is expected_kept


def test_d2_similarity_and_duplicate_checks_disagree_on_node_identity():
    """
    Root cause of D2: cycle_match_based_on_contig_id() strips orientation for the
    exact-duplicate test but keeps it for the similarity test, so the very same
    bubble is a duplicate in one comparison and unrelated in the other.
    """
    node_lengths = {"repeat+": 2500, "repeat-": 2500, "cargoB+": 1000, "cargoB-": 1000,
                    "cargoA+": 900, "cargoA-": 900}
    accepted = [["repeat+", "cargoA+"]]
    assert cycle_match_based_on_contig_id(["repeat+", "cargoB+"], node_lengths, accepted) is False
    # Same biological bubble, opposite strand -> shares no orientation-tagged node
    assert cycle_match_based_on_contig_id(["repeat-", "cargoB-"], node_lengths, accepted) is True


# ─── D3: k-mer dedup is set-based, so copy number is invisible ───────────────

@pytest.mark.xfail(strict=True, reason="D3: get_kmer_hashes returns a set, not a multiset")
def test_d3_full_composite_transposon_not_merged_into_its_half():
    """
    IS-cargo and IS-cargo-IS have almost identical k-mer *sets*. The complete
    composite transposon is the second one, and it is the one discarded --
    PICOTA keeps the partial structure and drops the real CT.
    """
    is_element = make_dna(1500, 10)
    cargo = make_dna(900, 11)
    partial = is_element + cargo
    complete = is_element + cargo + is_element

    kept = filter_cycles_with_kmer(
        [Cycle("a", partial, len(partial), 2, []),
         Cycle("b", complete, len(complete), 3, [])],
        K_MER_SIM, THRESHOLD_SIM, "Cycle")

    assert len(kept) == 2
    assert any(c.length == len(complete) for c in kept)


# ─── D4: the similarity denominator is asymmetric, so order decides ──────────

@pytest.mark.xfail(strict=True, reason="D4: common/len(new) makes containment order-dependent")
def test_d4_result_is_independent_of_input_order():
    """
    common_kmers / len(new_cycle_kmers) is a containment measure being used to
    make an identity decision: a small cycle contained in an already-accepted
    larger one always scores ~100% and is dropped, while the reverse order keeps
    both. Which composite transposon survives is decided by DFS traversal order.
    """
    shared = make_dna(1500, 10) + make_dna(900, 11)
    small = Cycle("small", shared, len(shared), 2, [])
    large_seq = shared + make_dna(2000, 12)
    large = Cycle("large", large_seq, len(large_seq), 2, [])

    def dedup(cycles):
        fresh = [Cycle(c.name, c.sequence, c.length, 2, []) for c in cycles]
        return len(filter_cycles_with_kmer(fresh, K_MER_SIM, THRESHOLD_SIM, "Cycle"))

    assert dedup([large, small]) == dedup([small, large])


# ─── D5: k-mer identity collapses into a cliff around 0.5% divergence ────────

@pytest.mark.parametrize("divergence,shared_fraction_k80", [
    (0.0001, 100.0),
    (0.001, 86.0),
    (0.005, 64.5),
    (0.01, 33.6),
    (0.02, 16.2),
])
def test_d5_kmer_identity_is_a_cliff_not_a_gradient(divergence, shared_fraction_k80):
    """
    Measured behaviour of k=80 against sequence divergence. Between 0.1% and 1%
    divergence -- exactly the range separating IS copies within one genome --
    shared k-mers fall from 86% to 34%, crossing the 80% dedup threshold with no
    usable middle ground.
    """
    seq = make_dna(4000, 20)
    mutated = mutate(seq, divergence, 7)
    kmers_a = get_kmer_hashes(seq, K_MER_SIM)
    kmers_b = get_kmer_hashes(mutated, K_MER_SIM)
    observed = len(kmers_a & kmers_b) / len(kmers_b) * 100
    assert observed == pytest.approx(shared_fraction_k80, abs=2.0)


def test_d5_shorter_k_widens_the_usable_range():
    """k=31 keeps two 0.5%-divergent copies above the threshold; k=80 does not."""
    seq = make_dna(4000, 20)
    mutated = mutate(seq, 0.005, 7)

    def shared(k):
        a, b = get_kmer_hashes(seq, k), get_kmer_hashes(mutated, k)
        return len(a & b) / len(b) * 100

    assert shared(80) < THRESHOLD_SIM
    assert shared(31) > THRESHOLD_SIM
