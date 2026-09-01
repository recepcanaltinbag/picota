"""
End-to-end test: implant composite transposons, sequence, assemble, detect.

Every other test in this suite checks a component. This one checks the claim the
tool actually makes -- put N composite transposons into a genome and PICOTA
finds them -- by running the real pipeline: simulate_ct_genome.py -> wgsim ->
SPAdes -> cycle_analysis -> BLAST against ground truth. A small genome keeps it
to a few seconds.

Skipped automatically when wgsim, spades.py or blastn are not on PATH.

Run:
    python3 -m pytest tests/test_e2e_simulated_ct.py -v
"""

import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "picota"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from score_picota_benchmark import read_ground_truth, run_blast, score  # noqa: E402
from src.cycle_finderv2 import cycle_analysis  # noqa: E402

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
PICOTA = os.path.join(os.path.dirname(__file__), "..", "picota")
IS_FASTA = os.path.join(PICOTA, "DBs", "ISes", "IS.fna")
CARGO_FASTA = os.path.join(PICOTA, "DBs", "Antibiotics",
                           "nucleotide_fasta_protein_homolog_model.fasta")

# ART is the publication-grade simulator (empirical per-cycle quality profiles);
# wgsim is accepted here only so the test can still run where ART is absent.
READ_SIMULATORS = ["art_illumina", "wgsim"]
REQUIRED_TOOLS = ["spades.py", "blastn", "makeblastdb"]


def available_simulator():
    for tool in READ_SIMULATORS:
        if shutil.which(tool):
            return tool
    return None

pytestmark = [
    pytest.mark.skipif(any(shutil.which(tool) is None for tool in REQUIRED_TOOLS),
                       reason=f"needs {', '.join(REQUIRED_TOOLS)} on PATH"),
    pytest.mark.skipif(available_simulator() is None,
                       reason=f"needs one of {', '.join(READ_SIMULATORS)} on PATH"),
    pytest.mark.skipif(not os.path.exists(IS_FASTA) or not os.path.exists(CARGO_FASTA),
                       reason="needs the bundled IS and CARD databases"),
]

# Deliberately small: one CT whose graph cycle lands under the default
# min_size_of_cycle and one comfortably above it.
N_CTS = 2
BACKBONE_LENGTH = 150000
COVERAGE = 30


@pytest.fixture(scope="module")
def simulated_assembly(tmp_path_factory):
    """Simulate, sequence and assemble once; every test reads the result."""
    work = tmp_path_factory.mktemp("e2e_ct")

    subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "simulate_ct_genome.py"),
         "--out-dir", str(work), "--backbone-length", str(BACKBONE_LENGTH),
         "--n-cts", str(N_CTS), "--shared-is", str(N_CTS),
         "--is-copies-outside", "3", "--cargo-genes", "1",
         "--spacing", "8000", "--seed", "5",
         "--is-fasta", IS_FASTA, "--cargo-fasta", CARGO_FASTA],
        check=True, capture_output=True)

    genome = work / "genome.fasta"
    length = sum(len(line.strip()) for line in open(genome)
                 if not line.startswith(">"))
    pairs = COVERAGE * length // (2 * 150)

    simulator = available_simulator()
    if simulator == "art_illumina":
        subprocess.run(
            ["art_illumina", "-ss", "HSXt", "-i", str(genome), "-p", "-l", "150",
             "-f", str(COVERAGE), "-m", "350", "-s", "50", "-rs", "1",
             "-o", str(work / "art_"), "-na"],
            check=True, capture_output=True)
        os.replace(str(work / "art_1.fq"), str(work / "r1.fq"))
        os.replace(str(work / "art_2.fq"), str(work / "r2.fq"))
    else:
        subprocess.run(
            ["wgsim", "-N", str(pairs), "-1", "150", "-2", "150", "-e", "0.001",
             "-r", "0", "-R", "0", "-X", "0", "-S", "1", str(genome),
             str(work / "r1.fq"), str(work / "r2.fq")],
            check=True, capture_output=True)

    subprocess.run(
        ["spades.py", "-1", str(work / "r1.fq"), "-2", str(work / "r2.fq"),
         "-o", str(work / "sp"), "-k", "77", "--only-assembler",
         "-t", "4", "-m", "8"],
        check=True, capture_output=True)

    gfa = work / "sp" / "assembly_graph_with_scaffolds.gfa"
    assert gfa.exists(), "SPAdes produced no assembly graph"
    return work, str(gfa)


def detect_and_score(work, gfa, dedup_mode, min_size):
    out = str(work / f"cycles_{dedup_mode}_{min_size}.fasta")
    cycle_analysis(gfa, out, True, 25, min_size, 40000, "Cycle", 1, 25, 80, 80,
                   dedup_mode=dedup_mode)
    ground_truth = read_ground_truth(str(work / "ground_truth.tsv"))
    cycle_ids = [line[1:].strip() for line in open(out) if line.startswith(">")]
    rows = run_blast(out, str(work / "ground_truth_cts.fasta"),
                     "blastn", "makeblastdb")
    return score(ground_truth, rows, cycle_ids, 95.0, 0.95)


class TestImplantedCTsAreFound:
    @pytest.mark.parametrize("dedup_mode", ["legacy", "strict"])
    def test_every_implanted_ct_is_recovered(self, simulated_assembly, dedup_mode):
        """The claim the tool makes, run end to end."""
        work, gfa = simulated_assembly
        result = detect_and_score(work, gfa, dedup_mode, min_size=1000)
        assert result["ct_recall"] == (N_CTS, N_CTS)

    @pytest.mark.parametrize("dedup_mode", ["legacy", "strict"])
    def test_cts_sharing_an_is_are_reported_separately(self, simulated_assembly,
                                                       dedup_mode):
        """
        Both implanted CTs are built on the same IS element, which the assembler
        collapses into one node. They must still come back as two candidates --
        this is the case the pipeline exists to handle.
        """
        work, gfa = simulated_assembly
        result = detect_and_score(work, gfa, dedup_mode, min_size=1000)
        assert result["copy_distinctness"] == (N_CTS, N_CTS)

    def test_recovered_cycles_are_shorter_than_the_ct(self, simulated_assembly):
        """
        Expected, not a defect: both flanking copies of the IS collapse into one
        graph node, so the cycle is IS + cargo where the CT is IS + cargo + IS.
        Pinned because it looks like a truncation bug to anyone reading output
        for the first time.
        """
        work, gfa = simulated_assembly
        result = detect_and_score(work, gfa, "strict", min_size=1000)
        ground_truth = read_ground_truth(str(work / "ground_truth.tsv"))
        for ct_id, (cycle, _) in result["matches"].items():
            assert cycle, f"{ct_id} not recovered"
            reported_length = int(cycle.split("-len")[1].split("-")[0])
            assert reported_length < int(ground_truth[ct_id]["CT_Length"])


class TestMinCycleSizeExcludesCompactCTs:
    """
    A real limitation of the shipped defaults, found by this test.

    config.yaml sets min_size_of_cycle: 2000. The graph cycle for a composite
    transposon is IS + cargo, so a compact CT -- IS26 at 820 bp plus a single
    resistance gene -- produces a cycle well under 2 kb and is dropped before
    anything scores it. Those are among the most common composite transposons
    in clinical isolates, so the default cannot be treated as neutral.
    """

    def test_compact_ct_is_missed_at_the_default_threshold(self, simulated_assembly):
        work, gfa = simulated_assembly
        result = detect_and_score(work, gfa, "strict", min_size=2000)
        recovered, total = result["ct_recall"]
        assert recovered < total, (
            "the compact CT is now found at the default threshold -- if the "
            "default changed, update this test and docs/ROADMAP.md")

    def test_lowering_the_threshold_recovers_it(self, simulated_assembly):
        work, gfa = simulated_assembly
        assert detect_and_score(work, gfa, "strict", min_size=1000)["ct_recall"] \
            == (N_CTS, N_CTS)

    def test_the_missed_cycle_is_between_the_two_thresholds(self, simulated_assembly):
        """Pin why it is missed, so a future failure is diagnosable."""
        work, gfa = simulated_assembly
        result = detect_and_score(work, gfa, "strict", min_size=1000)
        lengths = [int(cycle.split("-len")[1].split("-")[0])
                   for cycle, _ in result["matches"].values() if cycle]
        assert any(1000 <= length < 2000 for length in lengths)
