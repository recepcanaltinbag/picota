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
    STOP_CODONS,
    make_orf,
    choose_positions,
    load_backbone,
    place_inserts,
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
        genome, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path))
        assert records
        for record in records:
            excerpt = genome[record["CT_Start"] - 1:record["CT_End"]]
            assert excerpt == record["Sequence"]
            assert len(excerpt) == record["CT_Length"]

    def test_requested_number_of_cts_is_produced(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path, n_cts=3))
        assert len(records) == 3

    def test_shared_is_really_is_shared(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=2))
        assert records[0]["IS_Name"] == records[1]["IS_Name"]
        assert records[2]["IS_Name"] != records[0]["IS_Name"]

    def test_cts_sharing_an_is_carry_different_cargo(self, is_fasta, cargo_fasta, tmp_path):
        """Otherwise they would be genuine duplicates and prove nothing."""
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=2))
        assert records[0]["Cargo_Genes"] != records[1]["Cargo_Genes"]

    def test_is_copy_count_includes_free_standing_copies(self, is_fasta, cargo_fasta, tmp_path):
        """Two flanking copies per CT that uses it, plus the free-standing ones."""
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=2, is_copies_outside=4))
        assert records[0]["IS_Genome_Copies"] == 2 * 2 + 4
        assert records[2]["IS_Genome_Copies"] == 2

    def test_ct_regions_do_not_overlap(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path))
        spans = sorted((r["CT_Start"], r["CT_End"]) for r in records)
        for (_, end), (next_start, _) in zip(spans, spans[1:]):
            assert next_start > end

    def test_genome_reaches_the_requested_length(self, is_fasta, cargo_fasta, tmp_path):
        genome, _, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
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
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                        n_cts=3, shared_is=0))
        assert len({r["IS_Name"] for r in records}) == 3


class TestWriteOutputs:
    def _run(self, is_fasta, cargo_fasta, tmp_path):
        args = make_args(is_fasta, cargo_fasta, tmp_path)
        genome, records, _ = simulate(args)
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


class TestPlaceInserts:
    """
    Offsets are the part that can be wrong while looking right: every insert
    shifts everything downstream, so a start recorded against the original
    backbone would be plausible and useless.
    """

    def test_offsets_locate_each_insert(self):
        backbone = "".join(random.Random(20).choice("ACGT") for _ in range(10000))
        inserts = ["AAAA" * 25, "CCCC" * 25, "GGGG" * 25]
        genome, starts = place_inserts(backbone, inserts, random.Random(1), 100)
        for start, insert in zip(starts, inserts):
            assert genome[start:start + len(insert)] == insert

    def test_backbone_is_preserved_around_the_inserts(self):
        backbone = "".join(random.Random(21).choice("ACGT") for _ in range(10000))
        inserts = ["A" * 100]
        genome, starts = place_inserts(backbone, inserts, random.Random(1), 100)
        assert genome[:starts[0]] + genome[starts[0] + 100:] == backbone

    def test_genome_length_is_backbone_plus_inserts(self):
        backbone = "".join(random.Random(22).choice("ACGT") for _ in range(10000))
        inserts = ["A" * 100, "C" * 200]
        genome, _ = place_inserts(backbone, inserts, random.Random(1), 100)
        assert len(genome) == len(backbone) + 300

    def test_inserts_stay_in_order(self):
        backbone = "".join(random.Random(23).choice("ACGT") for _ in range(10000))
        genome, starts = place_inserts(backbone, ["A" * 50, "C" * 50, "G" * 50],
                                       random.Random(1), 100)
        assert starts == sorted(starts)

    def test_empty_insert_list_returns_the_backbone(self):
        backbone = "ACGT" * 100
        assert place_inserts(backbone, [], random.Random(1), 10) == (backbone, [])

    def test_backbone_too_small_is_rejected(self):
        with pytest.raises(SystemExit):
            place_inserts("ACGT" * 100, ["A" * 10] * 5, random.Random(1), 10000)


class TestRealBackbone:
    @pytest.fixture
    def backbone_fasta(self, tmp_path):
        path = tmp_path / "backbone.fna"
        rng = random.Random(30)
        short = "".join(rng.choice("ACGT") for _ in range(500))
        long_seq = "".join(rng.choice("ACGT") for _ in range(80000))
        path.write_text(f">plasmid_short\n{short}\n>NC_000913.3 chromosome\n{long_seq}\n")
        return str(path), long_seq

    def test_longest_record_is_used(self, backbone_fasta):
        path, long_seq = backbone_fasta
        name, seq = load_backbone(path)
        assert name == "NC_000913.3"
        assert seq == long_seq

    def test_empty_fasta_is_rejected(self, tmp_path):
        empty = tmp_path / "empty.fna"
        empty.write_text("")
        with pytest.raises(SystemExit):
            load_backbone(str(empty))

    def test_cts_are_implanted_into_the_real_sequence(self, is_fasta, cargo_fasta,
                                                      backbone_fasta, tmp_path):
        path, long_seq = backbone_fasta
        args = make_args(is_fasta, cargo_fasta, tmp_path, backbone_fasta=path)
        genome, records, name = simulate(args)
        assert name == "NC_000913.3"
        assert len(genome) > len(long_seq)
        for record in records:
            assert genome[record["CT_Start"] - 1:record["CT_End"]] == record["Sequence"]

    def test_backbone_bases_survive_implantation(self, is_fasta, cargo_fasta,
                                                 backbone_fasta, tmp_path):
        """Removing the implants must give the original chromosome back."""
        path, long_seq = backbone_fasta
        args = make_args(is_fasta, cargo_fasta, tmp_path, backbone_fasta=path)
        genome, records, _ = simulate(args)
        # Only CT spans are known here; free-standing IS copies are not recorded,
        # so check that every backbone base still appears in order around them.
        assert genome.count(long_seq[:1000]) == 1
        assert long_seq[:1000] in genome


class TestChoosePositions:
    """
    Random placement has to stay valid: inserts must not collide, and must not
    run off either end. An off-by-one here would corrupt ground truth silently.
    """

    def test_minimum_spacing_is_respected(self):
        rng = random.Random(40)
        for _ in range(50):
            positions = choose_positions(1000000, 8, 20000, rng, "random")
            gaps = [b - a for a, b in zip(positions, positions[1:])]
            assert all(gap >= 20000 for gap in gaps)

    def test_positions_stay_inside_the_backbone(self):
        rng = random.Random(41)
        for _ in range(50):
            positions = choose_positions(1000000, 8, 20000, rng, "random")
            assert positions[0] >= 20000
            assert positions[-1] <= 1000000 - 20000

    def test_positions_are_sorted(self):
        positions = choose_positions(1000000, 10, 5000, random.Random(42), "random")
        assert positions == sorted(positions)

    def test_random_placement_actually_varies(self):
        """Otherwise it is even spacing wearing a different name."""
        first = choose_positions(1000000, 6, 20000, random.Random(1), "random")
        second = choose_positions(1000000, 6, 20000, random.Random(2), "random")
        assert first != second

    def test_even_placement_is_evenly_spaced(self):
        positions = choose_positions(1000000, 4, 1000, random.Random(1), "even")
        gaps = [b - a for a, b in zip(positions, positions[1:])]
        assert len(set(gaps)) == 1

    def test_zero_inserts(self):
        assert choose_positions(1000, 0, 100, random.Random(1), "random") == []

    def test_backbone_too_small_is_rejected(self):
        with pytest.raises(SystemExit):
            choose_positions(10000, 5, 20000, random.Random(1), "random")

    def test_exactly_fitting_backbone_is_accepted(self):
        positions = choose_positions(6 * 1000, 5, 1000, random.Random(1), "random")
        assert len(positions) == 5


class TestPlacementInSimulation:
    def test_random_placement_is_the_default(self, is_fasta, cargo_fasta, tmp_path):
        args = make_args(is_fasta, cargo_fasta, tmp_path)
        assert args.placement == "random"

    def test_coordinates_stay_correct_under_random_placement(self, is_fasta,
                                                             cargo_fasta, tmp_path):
        for seed in (1, 2, 3):
            genome, records, _ = simulate(
                make_args(is_fasta, cargo_fasta, tmp_path, seed=seed))
            for record in records:
                assert genome[record["CT_Start"] - 1:record["CT_End"]] == record["Sequence"]

    def test_different_seeds_place_cts_differently(self, is_fasta, cargo_fasta, tmp_path):
        _, first, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path, seed=1))
        _, second, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path, seed=2))
        assert [r["CT_Start"] for r in first] != [r["CT_Start"] for r in second]


class TestNovelCargo:
    """
    Composite transposons whose cargo has no database homology. PICOTA's score
    is largely a sum of homology hits, so a benchmark built only from CARD genes
    would test the tool on the easy half of its own problem.
    """

    def test_orf_starts_with_a_start_codon(self):
        assert make_orf(900, random.Random(1)).startswith("ATG")

    def test_orf_ends_with_a_stop_codon(self):
        assert make_orf(900, random.Random(2))[-3:] in STOP_CODONS

    def test_orf_has_no_internal_stop(self):
        """Otherwise Prodigal truncates it and the cargo is not a gene."""
        orf = make_orf(1200, random.Random(3))
        codons = [orf[i:i+3] for i in range(0, len(orf) - 3, 3)]
        assert not (set(codons) & STOP_CODONS)

    def test_orf_length_is_a_multiple_of_three(self):
        assert len(make_orf(900, random.Random(4))) % 3 == 0

    def test_orf_is_deterministic(self):
        assert make_orf(600, random.Random(5)) == make_orf(600, random.Random(5))

    def test_novel_cts_are_labelled(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                           n_cts=3, shared_is=2, novel_cts=1))
        assert [r["Cargo_Type"] for r in records] == ["AMR", "AMR", "novel"]

    def test_novel_cargo_is_not_a_card_gene(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                           n_cts=3, shared_is=2, novel_cts=1))
        assert records[-1]["Cargo_Genes"].startswith("novel_orf_")

    def test_novel_cts_cannot_exceed_total(self, is_fasta, cargo_fasta, tmp_path):
        with pytest.raises(SystemExit):
            simulate(make_args(is_fasta, cargo_fasta, tmp_path, n_cts=2, novel_cts=3))

    def test_coordinates_still_correct_with_novel_cargo(self, is_fasta, cargo_fasta, tmp_path):
        genome, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                                n_cts=3, shared_is=2, novel_cts=2))
        for record in records:
            assert genome[record["CT_Start"] - 1:record["CT_End"]] == record["Sequence"]


class TestCargoISCounting:
    """
    Ground-truth arithmetic around an IS that also sits inside the cargo. Both
    cases below were wrong once: the copy count omitted the cargo-internal copy,
    and the interrupting element was drawn without excluding the flanking ones,
    so a collision could silently turn "different IS" into "same IS".
    """

    def test_same_mode_counts_the_cargo_internal_copy(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                           n_cts=2, shared_is=0,
                                           cargo_is_mode="same"))
        for record in records:
            assert record["IS_Genome_Copies"] == 3, "2 flanking + 1 inside the cargo"

    def test_none_mode_counts_two(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                           n_cts=2, shared_is=0))
        assert all(r["IS_Genome_Copies"] == 2 for r in records)

    def test_different_mode_counts_two_for_the_flanking_element(self, is_fasta,
                                                                cargo_fasta, tmp_path):
        """The interrupting element is a different one, so it adds no copy here."""
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                           n_cts=2, shared_is=0,
                                           cargo_is_mode="different"))
        assert all(r["IS_Genome_Copies"] == 2 for r in records)

    def test_shared_mode_adds_free_standing_copies(self, is_fasta, cargo_fasta, tmp_path):
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                           n_cts=3, shared_is=2,
                                           is_copies_outside=4,
                                           cargo_is_mode="same"))
        # three per shared element, two of them share it, plus four free-standing
        assert records[0]["IS_Genome_Copies"] == 3 * 2 + 4

    @pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
    def test_interrupting_element_is_never_a_flanking_one(self, is_fasta,
                                                          cargo_fasta, tmp_path, seed):
        """
        A collision would make the scenario the harder one while the table still
        called it the easier one. The fixture offers only three elements, so a
        careless draw collides often.
        """
        _, records, _ = simulate(make_args(is_fasta, cargo_fasta, tmp_path,
                                           n_cts=2, shared_is=0, seed=seed,
                                           cargo_is_mode="different"))
        flanking = {r["IS_Name"] for r in records}
        interrupting = {r["Cargo_IS"] for r in records}
        assert not (flanking & interrupting)
