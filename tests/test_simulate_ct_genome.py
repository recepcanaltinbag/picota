"""
Tests for scripts/simulate_ct_genome.py (roadmap phase 0.5).

The benchmark's entire value rests on the ground truth being exactly right, so
the coordinate and copy-count invariants are tested hardest: if
genome[CT_Start-1:CT_End] is not the CT sequence, every recall number computed
downstream is wrong in a way nothing else would catch.
"""

import json
import os
import random
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from simulate_ct_genome import (  # noqa: E402
    GROUND_TRUTH_COLUMNS,
    build_composite_transposon,
    build_parser,
    is_clean_dna,
    mutate,
    parse_card_name,
    parse_is_name,
    random_backbone,
    read_fasta,
    reverse_complement,
    simulate,
    write_outputs,
)


@pytest.fixture
def is_fasta(tmp_path):
    """Three ISfinder-shaped records inside the usable length band."""
    path = tmp_path / "is.fna"
    rng = random.Random(7)
    lines = []
    for name, family, length in [("ISfoo1", "IS6", 900),
                                 ("ISbar2", "IS3", 1200),
                                 ("ISbaz3", "IS110", 1500)]:
        lines.append(f">{name}_{family}_grp")
        lines.append("".join(rng.choice("ACGT") for _ in range(length)))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


@pytest.fixture
def cargo_fasta(tmp_path):
    """CARD-shaped records."""
    path = tmp_path / "card.fna"
    rng = random.Random(8)
    lines = []
    for i, gene in enumerate(["blaTEM-1", "aadA1", "tetA", "sul1", "dfrA7", "ermB"]):
        lines.append(f">gb|ACC{i}|+|0-800|ARO:300000{i}|{gene} [Escherichia coli]")
        lines.append("".join(rng.choice("ACGT") for _ in range(800)))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def make_args(is_fasta, cargo_fasta, tmp_path, **overrides):
    argv = ["--out-dir", str(tmp_path / "out"),
            "--is-fasta", is_fasta, "--cargo-fasta", cargo_fasta,
            "--backbone-length", "60000", "--n-cts", "3", "--shared-is", "2",
            "--is-copies-outside", "2", "--is-divergence", "0.5",
            "--cargo-genes", "1", "--spacing", "2000", "--seed", "1"]
    for key, value in overrides.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    return build_parser().parse_args(argv)


class TestSequenceHelpers:
    def test_reverse_complement(self):
        assert reverse_complement("ACGTN") == "NACGT"

    def test_reverse_complement_is_an_involution(self):
        seq = "ACGTTGCAAC"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_is_clean_dna_rejects_ambiguity(self):
        assert is_clean_dna("ACGT") is True
        assert is_clean_dna("ACGTN") is False
        assert is_clean_dna("") is False

    def test_mutate_at_zero_is_identity(self):
        seq = "ACGT" * 100
        assert mutate(seq, 0, random.Random(1)) == seq

    def test_mutate_changes_roughly_the_requested_fraction(self):
        seq = "".join(random.Random(2).choice("ACGT") for _ in range(10000))
        mutated = mutate(seq, 1.0, random.Random(3))
        changed = sum(1 for a, b in zip(seq, mutated) if a != b)
        assert changed == pytest.approx(100, rel=0.4)

    def test_mutate_never_substitutes_a_base_for_itself(self):
        seq = "A" * 1000
        assert "A" not in mutate(seq, 100.0, random.Random(4))

    def test_random_backbone_respects_gc_content(self):
        seq = random_backbone(20000, 0.65, random.Random(5))
        gc = sum(1 for base in seq if base in "GC") / len(seq)
        assert gc == pytest.approx(0.65, abs=0.02)


class TestHeaderParsing:
    def test_isfinder_header(self):
        assert parse_is_name("IS26_IS6_IS6") == ("IS26", "IS6")

    def test_isfinder_header_with_variant(self):
        assert parse_is_name("ISBj2_B_IS5_IS5") == ("ISBj2", "IS5")

    def test_unparseable_isfinder_header(self):
        assert parse_is_name("weird") == ("weird", "unknown")

    def test_card_header(self):
        header = "gb|HQ845196.1|+|0-861|ARO:3001109|SHV-52 [Klebsiella pneumoniae]"
        assert parse_card_name(header) == "SHV-52"

    def test_unparseable_card_header(self):
        assert parse_card_name("plain_name rest") == "plain_name"


class TestReadFasta:
    def test_reads_records(self, is_fasta):
        records = list(read_fasta(is_fasta))
        assert len(records) == 3
        assert records[0][0].startswith("ISfoo1")

    def test_length_filters_applied(self, is_fasta):
        assert len(list(read_fasta(is_fasta, min_length=1000))) == 2
        assert len(list(read_fasta(is_fasta, min_length=0, max_length=1000))) == 1


class TestBuildCompositeTransposon:
    def test_structure_is_is_cargo_is(self):
        rng = random.Random(9)
        is_seq = "".join(rng.choice("ACGT") for _ in range(500))
        cargo = "".join(rng.choice("ACGT") for _ in range(300))
        built, cargo_out = build_composite_transposon(is_seq, [cargo], 0, "direct", rng)
        assert built.startswith(is_seq)
        assert built.endswith(is_seq)
        assert cargo in built
        assert cargo_out == cargo

    def test_inverted_orientation_flips_the_right_copy(self):
        rng = random.Random(10)
        is_seq = "".join(rng.choice("ACGT") for _ in range(500))
        cargo = "".join(rng.choice("ACGT") for _ in range(300))
        built, _ = build_composite_transposon(is_seq, [cargo], 0, "inverted", rng)
        assert built.startswith(is_seq)
        assert built.endswith(reverse_complement(is_seq))

    def test_divergence_makes_the_two_copies_differ(self):
        rng = random.Random(11)
        is_seq = "".join(rng.choice("ACGT") for _ in range(1000))
        cargo = "".join(rng.choice("ACGT") for _ in range(300))
        built, _ = build_composite_transposon(is_seq, [cargo], 2.0, "direct", rng)
        assert not built.startswith(is_seq)

    def test_multiple_cargo_genes_are_all_present(self):
        rng = random.Random(12)
        is_seq = "".join(rng.choice("ACGT") for _ in range(400))
        cargo = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(3)]
        built, _ = build_composite_transposon(is_seq, cargo, 0, "direct", rng)
        assert all(gene in built for gene in cargo)


class TestGroundTruth:
    def test_coordinates_locate_the_ct_in_the_genome(self, is_fasta, cargo_fasta, tmp_path):
        """
        The invariant the whole benchmark rests on. If this drifts, every recall
        number computed downstream is quietly wrong.
        """
        genome, records = simulate(make_args(is_fasta, cargo_fasta, tmp_path))
        assert records
        for record in records:
            excerpt = genome[record["CT_Start"] - 1:record["CT_End"]]
            assert excerpt == record["Sequence"]
            assert len(excerpt) == record["CT_Length"]

    def test_requested_number_of_cts_is_produced(self, is_fasta, cargo_fasta, tmp_path):
        _, records = simulate(make_args(is_fasta, cargo_fasta, tmp_path, n_cts=3))
        assert len(records) == 3

    def test_shared_is_really_is_shared(self, is_fasta, cargo_fasta, tmp_path):
        _, records = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=2))
        assert records[0]["IS_Name"] == records[1]["IS_Name"]
        assert records[2]["IS_Name"] != records[0]["IS_Name"]

    def test_cts_sharing_an_is_carry_different_cargo(self, is_fasta, cargo_fasta, tmp_path):
        """Otherwise they would be genuine duplicates and prove nothing."""
        _, records = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=2))
        assert records[0]["Cargo_Genes"] != records[1]["Cargo_Genes"]

    def test_is_copy_count_includes_free_standing_copies(self, is_fasta, cargo_fasta, tmp_path):
        """Two flanking copies per CT that uses it, plus the free-standing ones."""
        _, records = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=2, is_copies_outside=4))
        assert records[0]["IS_Genome_Copies"] == 2 * 2 + 4
        assert records[2]["IS_Genome_Copies"] == 2

    def test_ct_regions_do_not_overlap(self, is_fasta, cargo_fasta, tmp_path):
        _, records = simulate(make_args(is_fasta, cargo_fasta, tmp_path))
        spans = sorted((r["CT_Start"], r["CT_End"]) for r in records)
        for (_, end), (next_start, _) in zip(spans, spans[1:]):
            assert next_start > end

    def test_genome_reaches_the_requested_length(self, is_fasta, cargo_fasta, tmp_path):
        genome, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                       backbone_length=60000))
        assert len(genome) >= 60000

    def test_seed_makes_the_run_reproducible(self, is_fasta, cargo_fasta, tmp_path):
        first = simulate(make_args(is_fasta, cargo_fasta, tmp_path, seed=42))
        second = simulate(make_args(is_fasta, cargo_fasta, tmp_path, seed=42))
        assert first[0] == second[0]

    def test_different_seeds_give_different_genomes(self, is_fasta, cargo_fasta, tmp_path):
        first = simulate(make_args(is_fasta, cargo_fasta, tmp_path, seed=1))
        second = simulate(make_args(is_fasta, cargo_fasta, tmp_path, seed=2))
        assert first[0] != second[0]

    def test_shared_is_cannot_exceed_n_cts(self, is_fasta, cargo_fasta, tmp_path):
        with pytest.raises(SystemExit):
            simulate(make_args(is_fasta, cargo_fasta, tmp_path, n_cts=2, shared_is=3))

    def test_no_shared_is_gives_every_ct_its_own(self, is_fasta, cargo_fasta, tmp_path):
        _, records = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=0))
        assert len({r["IS_Name"] for r in records}) == 3


class TestWriteOutputs:
    def _run(self, is_fasta, cargo_fasta, tmp_path):
        args = make_args(is_fasta, cargo_fasta, tmp_path)
        genome, records = simulate(args)
        return genome, records, write_outputs(args.out_dir, genome, records, vars(args))

    def test_all_four_files_written(self, is_fasta, cargo_fasta, tmp_path):
        _, _, paths = self._run(is_fasta, cargo_fasta, tmp_path)
        assert all(os.path.exists(path) for path in paths)

    def test_genome_fasta_round_trips(self, is_fasta, cargo_fasta, tmp_path):
        genome, _, (genome_path, _, _, _) = self._run(is_fasta, cargo_fasta, tmp_path)
        written = "".join(line.strip() for line in open(genome_path)
                          if not line.startswith(">"))
        assert written == genome

    def test_tsv_has_the_documented_columns(self, is_fasta, cargo_fasta, tmp_path):
        _, records, (_, tsv_path, _, _) = self._run(is_fasta, cargo_fasta, tmp_path)
        rows = [line.rstrip("\n").split("\t") for line in open(tsv_path)]
        assert rows[0] == GROUND_TRUTH_COLUMNS
        assert len(rows) == len(records) + 1

    def test_sequence_column_is_not_leaked_into_the_tsv(self, is_fasta, cargo_fasta, tmp_path):
        """The sequences belong in the FASTA; a whole CT per TSV cell is unusable."""
        _, _, (_, tsv_path, _, _) = self._run(is_fasta, cargo_fasta, tmp_path)
        assert "Sequence" not in open(tsv_path).readline()

    def test_ct_fasta_matches_the_ground_truth_records(self, is_fasta, cargo_fasta, tmp_path):
        _, records, (_, _, cts_path, _) = self._run(is_fasta, cargo_fasta, tmp_path)
        written = {}
        name = None
        for line in open(cts_path):
            if line.startswith(">"):
                name = line[1:].strip()
                written[name] = ""
            else:
                written[name] += line.strip()
        assert len(written) == len(records)
        for record in records:
            key = next(k for k in written if k.startswith(record["CT_ID"]))
            assert written[key] == record["Sequence"]

    def test_json_records_parameters_and_length(self, is_fasta, cargo_fasta, tmp_path):
        genome, records, (_, _, _, json_path) = self._run(is_fasta, cargo_fasta, tmp_path)
        payload = json.load(open(json_path))
        assert payload["genome_length"] == len(genome)
        assert len(payload["composite_transposons"]) == len(records)
        assert payload["parameters"]["seed"] == 1
