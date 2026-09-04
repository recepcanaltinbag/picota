"""
Unit tests for scoringv4ProtBlast.py

Tested fonksiyonlar:
  - CodingRegion class
  - GeneticInfo class
  - calculate_total_score()
  - merge_intervals()
  - parsing_blast_file()
  - parsing_blast_file_grouped() / parsing_blast_file_merged_grouped()
  - concat_fasta()
  - xenobiotic_names()
  - src.stale_output
"""

import sys
import os
import tempfile
import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.scoringv4ProtBlast import (
    CodingRegion,
    GeneticInfo,
    calculate_total_score,
    merge_intervals,
    parsing_blast_file,
    parsing_blast_file_grouped,
    parsing_blast_file_merged_grouped,
    concat_fasta,
    xenobiotic_names,
)
from src import stale_output


# ─────────────────────────────────────────────
# CodingRegion
# ─────────────────────────────────────────────

class TestCodingRegion:
    def test_basic_creation(self):
        cr = CodingRegion(10, 200, 1, "gene|product|name", "Antibiotics", 85.0)
        assert cr.start == 10
        assert cr.end == 200
        assert cr.strand == 1
        assert cr.r_type == "Antibiotics"
        assert cr.score == 85.0
        assert cr.product == ''
        assert cr.gene == ''

    def test_reverse_strand(self):
        cr = CodingRegion(500, 100, -1, "gene", "InsertionSequences", 90.0)
        assert cr.strand == -1

    def test_fullname_stored(self):
        cr = CodingRegion(0, 0, 1, "db|product|gene_name", "Xenobiotics", 50.0)
        assert cr.fullname == "db|product|gene_name"


# ─────────────────────────────────────────────
# GeneticInfo
# ─────────────────────────────────────────────

class TestGeneticInfo:
    def test_basic_creation(self):
        gi = GeneticInfo("acc1", "qseq1", "description", [], "ATCG", 10.0, 20.0, 30.0)
        assert gi.seq_acc == "acc1"
        assert gi.seq_id == "qseq1"
        assert gi.feature_list == []
        assert gi.score0 == 10.0
        assert gi.score1 == 20.0
        assert gi.score2 == 30.0

    def test_feature_list_stored(self):
        cr = CodingRegion(1, 100, 1, "x", "Antibiotics", 80.0)
        gi = GeneticInfo("acc", "q", "desc", [cr], "ATCG", 0, 0, 0)
        assert len(gi.feature_list) == 1


# ─────────────────────────────────────────────
# calculate_total_score
# ─────────────────────────────────────────────

class TestCalculateTotalScore:
    """
    Parametreler:
      total_score_type: 0, 1 veya 2
      dist_type: 1 veya diğer
      max_z: float (ör. 20)
      mean_of_CompTns, std_of_CompTns: float
      len_of_cycle: int
      lst_ant, lst_is, lst_xe: list of floats
      comp_number: int
    """

    BASE = dict(
        dist_type=0,
        max_z=20.0,
        mean_of_CompTns=5000.0,
        std_of_CompTns=1000.0,
        len_of_cycle=5000,
        lst_ant=[80.0],
        lst_is=[90.0],
        lst_xe=[],
        comp_number=2,
    )

    def test_type0_returns_float(self):
        score = calculate_total_score(total_score_type=0, **self.BASE)
        assert isinstance(score, float)

    def test_type1_returns_float(self):
        score = calculate_total_score(total_score_type=1, **self.BASE)
        assert isinstance(score, float)

    def test_type2_returns_float(self):
        score = calculate_total_score(total_score_type=2, **self.BASE)
        assert isinstance(score, float)

    def test_invalid_type_raises(self):
        with pytest.raises(Exception):
            calculate_total_score(total_score_type=99, **self.BASE)

    def test_empty_lists_score_type0(self):
        params = dict(self.BASE)
        params.update(lst_ant=[], lst_is=[], lst_xe=[])
        score = calculate_total_score(total_score_type=0, **params)
        # 0^z_c_l = 0
        assert score == 0.0

    def test_empty_lists_score_type1(self):
        params = dict(self.BASE)
        params.update(lst_ant=[], lst_is=[], lst_xe=[])
        score = calculate_total_score(total_score_type=1, **params)
        assert score == 0.0

    def test_type2_with_is_only(self):
        # Sadece IS var, ant/xe yok → antcxc=0, isc=1 → (0*1) + 10^z_c_l
        params = dict(self.BASE)
        params.update(lst_ant=[], lst_is=[90.0], lst_xe=[])
        score = calculate_total_score(total_score_type=2, **params)
        assert score > 0

    def test_dist_type1_short_cycle_zero_z(self):
        # len_of_cycle < mean → dist_type=1 → z=0
        score_dist1 = calculate_total_score(
            total_score_type=0, dist_type=1, max_z=20.0,
            mean_of_CompTns=5000.0, std_of_CompTns=1000.0,
            len_of_cycle=3000,
            lst_ant=[80.0], lst_is=[90.0], lst_xe=[],
            comp_number=2
        )
        score_dist0 = calculate_total_score(
            total_score_type=0, dist_type=0, max_z=20.0,
            mean_of_CompTns=5000.0, std_of_CompTns=1000.0,
            len_of_cycle=3000,
            lst_ant=[80.0], lst_is=[90.0], lst_xe=[],
            comp_number=2
        )
        # dist_type=1 kısa döngüde z=0 demek → daha yüksek z_c_l → yüksek skor
        assert score_dist1 >= score_dist0

    def test_z_capped_at_max_z(self):
        # Çok uzun döngü → z büyük ama max_z'de kesiyor olmalı
        score_huge = calculate_total_score(
            total_score_type=0, dist_type=0, max_z=20.0,
            mean_of_CompTns=5000.0, std_of_CompTns=1000.0,
            len_of_cycle=100000,
            lst_ant=[80.0], lst_is=[90.0], lst_xe=[],
            comp_number=2
        )
        score_normal = calculate_total_score(
            total_score_type=0, dist_type=0, max_z=20.0,
            mean_of_CompTns=5000.0, std_of_CompTns=1000.0,
            len_of_cycle=25000,
            lst_ant=[80.0], lst_is=[90.0], lst_xe=[],
            comp_number=2
        )
        # Her ikisi de max_z'e takılmalı → aynı skor
        assert score_huge == score_normal

    def test_score_positive_with_hits(self):
        score = calculate_total_score(total_score_type=0, **self.BASE)
        assert score > 0


# ─────────────────────────────────────────────
# merge_intervals
# ─────────────────────────────────────────────

class TestMergeIntervals:
    def test_empty_list(self):
        assert merge_intervals([]) == []

    def test_single_interval(self):
        assert merge_intervals([(1, 5)]) == [(1, 5)]

    def test_non_overlapping(self):
        result = merge_intervals([(1, 3), (5, 8)])
        assert result == [(1, 3), (5, 8)]

    def test_overlapping(self):
        result = merge_intervals([(1, 5), (3, 8)])
        assert result == [(1, 8)]

    def test_touching_intervals(self):
        # (1,5) ve (5,10) → merge edilmeli mi? ≤ kontrolüne bağlı
        result = merge_intervals([(1, 5), (5, 10)])
        # Touching (end == start) → merge
        assert result == [(1, 10)]

    def test_fully_contained(self):
        result = merge_intervals([(1, 10), (3, 7)])
        assert result == [(1, 10)]

    def test_multiple_merges(self):
        result = merge_intervals([(1, 4), (2, 6), (5, 9), (15, 20)])
        assert result == [(1, 9), (15, 20)]

    def test_unsorted_input(self):
        result = merge_intervals([(10, 15), (1, 5), (3, 8)])
        assert result == [(1, 8), (10, 15)]

    def test_negative_coordinates(self):
        result = merge_intervals([(-10, -5), (-7, -2)])
        assert result == [(-10, -2)]

    def test_all_same(self):
        result = merge_intervals([(3, 7), (3, 7), (3, 7)])
        assert result == [(3, 7)]


# ─────────────────────────────────────────────
# parsing_blast_file
# ─────────────────────────────────────────────

def _write_blast_file(rows):
    """
    rows: list of tuples (qseqid, sseqid, pident, length, mismatch, gapopen,
                          qstart, qend, sstart, send, evalue, bitscore, slen, qlen)
    """
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    for row in rows:
        f.write('\t'.join(str(x) for x in row) + '\n')
    f.flush()
    return f.name


class TestParsingBlastFile:
    def test_empty_file_returns_empty_list(self):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
        f.close()
        result = parsing_blast_file(f.name, 'Antibiotics', 50.0, {})
        os.unlink(f.name)
        assert result == []

    def test_basic_antibiotic_hit(self):
        # pident=90, evalue=1e-5, slen=300, length=270 → score = (270/300)*90 = 81
        rows = [("q1_1", "db|product|gene", 90.0, 270, 0, 0, 1, 90, 1, 90, 1e-5, 200, 300, 270)]
        blast_file = _write_blast_file(rows)
        info_prod_dict = {"q1_1": (0, 1)}   # offset=0
        result = parsing_blast_file(blast_file, 'Antibiotics', 50.0, info_prod_dict)
        os.unlink(blast_file)
        assert len(result) == 1
        assert result[0].r_type == 'Antibiotics'

    def test_low_pident_filtered_out(self):
        # pident=50 < 80 → filtrelenmeli
        rows = [("q1_1", "db|product|gene", 50.0, 270, 0, 0, 1, 90, 1, 90, 1e-5, 200, 300, 270)]
        blast_file = _write_blast_file(rows)
        result = parsing_blast_file(blast_file, 'Antibiotics', 50.0, {"q1_1": (0, 1)})
        os.unlink(blast_file)
        assert result == []

    def test_high_evalue_filtered_out(self):
        # evalue=1e11 > 1e10 → filtrelenmeli
        rows = [("q1_1", "db|product|gene", 90.0, 270, 0, 0, 1, 90, 1, 90, 1e11, 200, 300, 270)]
        blast_file = _write_blast_file(rows)
        result = parsing_blast_file(blast_file, 'Antibiotics', 50.0, {"q1_1": (0, 1)})
        os.unlink(blast_file)
        assert result == []

    def test_below_threshold_filtered_out(self):
        # score = (10/300)*90 = 3 < threshold=50 → filtrelenmeli
        rows = [("q1_1", "db|product|gene", 90.0, 10, 0, 0, 1, 10, 1, 10, 1e-5, 50, 300, 10)]
        blast_file = _write_blast_file(rows)
        result = parsing_blast_file(blast_file, 'Antibiotics', 50.0, {"q1_1": (0, 1)})
        os.unlink(blast_file)
        assert result == []

    def test_insertion_sequence_no_coord_conversion(self):
        # IS için koordinat dönüşümü yapılmamalı (protein→nucleotide)
        rows = [("q1", "IS26", 95.0, 800, 0, 0, 100, 900, 100, 900, 1e-50, 500, 820, 900)]
        blast_file = _write_blast_file(rows)
        result = parsing_blast_file(blast_file, 'InsertionSequences', 10.0, {})
        os.unlink(blast_file)
        assert len(result) == 1
        assert result[0].start == 100
        assert result[0].end == 900

    def test_reverse_strand_detected(self):
        # qstart > qend → strand = -1, koordinatlar swap edilmeli
        rows = [("q1_1", "db|product|gene", 95.0, 270, 0, 0, 90, 1, 90, 1, 1e-5, 200, 300, 90)]
        blast_file = _write_blast_file(rows)
        info_prod_dict = {"q1_1": (0, 1)}
        result = parsing_blast_file(blast_file, 'Antibiotics', 50.0, info_prod_dict)
        os.unlink(blast_file)
        assert len(result) == 1
        assert result[0].strand == -1
        assert result[0].start < result[0].end

    def test_best_hit_selected_per_query(self):
        # Aynı qseqid için iki hit → yüksek skorlu seçilmeli
        rows = [
            ("q1_1", "gene_low",  90.0, 100, 0, 0, 1, 100, 1, 100, 1e-5, 100, 300, 100),
            ("q1_1", "gene_high", 90.0, 270, 0, 0, 1, 270, 1, 270, 1e-5, 200, 300, 270),
        ]
        blast_file = _write_blast_file(rows)
        result = parsing_blast_file(blast_file, 'Antibiotics', 30.0, {"q1_1": (0, 1)})
        os.unlink(blast_file)
        assert len(result) == 1
        assert result[0].fullname == "gene_high"

    def test_fullname_stored_correctly(self):
        # parsing_blast_file fullname'i saklar; product/gene sonraki adımda set edilir
        rows = [("q1_1", "db|MyProduct|MyGene", 90.0, 270, 0, 0, 1, 90, 1, 90, 1e-5, 200, 300, 270)]
        blast_file = _write_blast_file(rows)
        result = parsing_blast_file(blast_file, 'Antibiotics', 50.0, {"q1_1": (0, 1)})
        os.unlink(blast_file)
        assert len(result) == 1
        assert result[0].fullname == "db|MyProduct|MyGene"
        # product/gene bu aşamada boş string (sonraki scoring adımında doldurulur)
        assert result[0].product == ''
        assert result[0].gene == ''

    def test_fullname_no_pipe_stored(self):
        rows = [("q1_1", "NoPipeFullName", 90.0, 270, 0, 0, 1, 90, 1, 90, 1e-5, 200, 300, 270)]
        blast_file = _write_blast_file(rows)
        result = parsing_blast_file(blast_file, 'Antibiotics', 50.0, {"q1_1": (0, 1)})
        os.unlink(blast_file)
        assert len(result) == 1
        assert result[0].fullname == "NoPipeFullName"


# ─────────────────────────────────────────────
# Batched BLAST: grouping and concatenation
# ─────────────────────────────────────────────

class TestGroupedParsing:
    """
    Batched searching puts every cycle's queries in one result file, so the
    parsers have to attribute each hit back to its cycle. The unbatched
    functions are thin wrappers over these, which is what keeps the two paths
    from drifting apart.
    """

    def test_hits_go_to_the_cycle_their_query_came_from(self):
        rows = [
            ("cycA_1", "db|product|geneA", 90.0, 270, 0, 0, 1, 90, 1, 90, 1e-5, 200, 300, 270),
            ("cycB_1", "db|product|geneB", 95.0, 280, 0, 0, 1, 93, 1, 93, 1e-5, 200, 300, 280),
        ]
        blast_file = _write_blast_file(rows)
        info_prod_dict = {"cycA_1": (0, 1), "cycB_1": (0, 1)}
        owner = {"cycA_1": "cycA", "cycB_1": "cycB"}
        grouped = parsing_blast_file_grouped(
            blast_file, 'Antibiotics', 50.0, info_prod_dict, owner)
        os.unlink(blast_file)

        assert set(grouped) == {"cycA", "cycB"}
        assert grouped["cycA"][0].fullname == "db|product|geneA"
        assert grouped["cycB"][0].fullname == "db|product|geneB"

    def test_grouped_without_owner_matches_the_unbatched_parser(self):
        rows = [
            ("q1_1", "db|product|gene", 90.0, 270, 0, 0, 1, 90, 1, 90, 1e-5, 200, 300, 270),
            ("q1_2", "db|other|gene2", 85.0, 250, 0, 0, 1, 83, 1, 83, 1e-5, 180, 280, 250),
        ]
        blast_file = _write_blast_file(rows)
        info_prod_dict = {"q1_1": (0, 1), "q1_2": (0, 1)}
        flat = parsing_blast_file(blast_file, 'Antibiotics', 50.0, info_prod_dict)
        grouped = parsing_blast_file_grouped(
            blast_file, 'Antibiotics', 50.0, info_prod_dict)
        os.unlink(blast_file)

        assert [r.fullname for r in grouped[None]] == [r.fullname for r in flat]

    def test_best_transposon_is_chosen_per_cycle_not_per_file(self):
        """
        The merged parser keeps one hit only. Resolved file-wide, as it was when
        a file held one cycle, the strongest cycle's transposon call would be
        handed to every cycle in a batched run and the rest would lose theirs.
        """
        rows = [
            # cycA: a strong Tn5 hit
            ("cycA", "Tn5", 99.0, 900, 0, 0, 1, 900, 1, 900, 1e-50, 900, 900, 5000),
            # cycB: a weaker but still reportable Tn3 hit -- (850/900)*95 = 89.7
            ("cycB", "Tn3", 95.0, 850, 0, 0, 1, 850, 1, 850, 1e-40, 800, 900, 5000),
        ]
        blast_file = _write_blast_file(rows)
        grouped = parsing_blast_file_merged_grouped(
            blast_file, 'CompTNs', 80.0, {}, {"cycA": "cycA", "cycB": "cycB"})
        os.unlink(blast_file)

        assert set(grouped) == {"cycA", "cycB"}
        assert grouped["cycA"][0].fullname == "Tn5"
        assert grouped["cycB"][0].fullname == "Tn3"

    def test_merged_keeps_the_single_best_subject_within_one_cycle(self):
        rows = [
            ("cycA", "Tn5", 99.0, 900, 0, 0, 1, 900, 1, 900, 1e-50, 900, 900, 5000),
            ("cycA", "Tn3", 95.0, 850, 0, 0, 1, 850, 1, 850, 1e-40, 800, 900, 5000),
        ]
        blast_file = _write_blast_file(rows)
        grouped = parsing_blast_file_merged_grouped(
            blast_file, 'CompTNs', 80.0, {}, {"cycA": "cycA"})
        os.unlink(blast_file)

        assert len(grouped["cycA"]) == 1
        assert grouped["cycA"][0].fullname == "Tn5"


class TestConcatFasta:
    def test_missing_trailing_newline_does_not_fuse_records(self, tmp_path):
        """split_fasta writes its records without a trailing newline."""
        a = tmp_path / "a.fa"
        b = tmp_path / "b.fa"
        a.write_text(">cycA\nACGTACGT")
        b.write_text(">cycB\nTTTTGGGG")
        merged = tmp_path / "merged" / "all.fa"

        written = concat_fasta([str(a), str(b)], str(merged))

        assert written == 2
        assert merged.read_text() == ">cycA\nACGTACGT\n>cycB\nTTTTGGGG\n"

    def test_empty_and_missing_files_are_skipped(self, tmp_path):
        present = tmp_path / "a.fa"
        present.write_text(">cycA\nACGT\n")
        empty = tmp_path / "empty.fa"
        empty.write_text("")
        merged = tmp_path / "all.fa"

        written = concat_fasta(
            [str(present), str(empty), str(tmp_path / "gone.fa")], str(merged))

        assert written == 1
        assert merged.read_text() == ">cycA\nACGT\n"


# ─────────────────────────────────────────────
# Xenobiotics hit names
# ─────────────────────────────────────────────

class TestXenobioticNames:
    """
    The KEGG-built reference set names a sequence with its KO, EC and pathway.
    The classified set that preceded it has no structure at all, and both have
    to keep working, so the format is read off the name rather than configured.
    """

    KEGG = ("P9WQP7|KO:K16045|EC:1.1.1.145|PATH:map00984|"
            "3_beta-hydroxysteroid_dehydrogenase|Mycobacterium_tuberculosis")
    CLASSIFIED = ("QIH10984.1_aromatic_ring-hydroxylating_dioxygenase_"
                  "subunit_alpha_[Pseudomonas_sp._BIOMIG1BAC]")

    def test_kegg_product_carries_the_function_not_the_accession(self):
        """
        Regression: the accession alone was all that reached the result table.
        BLAST cuts sseqid at the first space, so the first build's descriptions
        -- written with spaces -- were truncated mid-word in all 16,539 headers
        and 'P9WQP7' was reported as the product.
        """
        product, gene = xenobiotic_names(self.KEGG)
        assert product.startswith("3_beta-hydroxysteroid_dehydrogenase")
        assert "KO:K16045" in product
        assert "EC:1.1.1.145" in product
        assert gene == "K16045"

    def test_classified_names_pass_through_unchanged(self):
        product, gene = xenobiotic_names(self.CLASSIFIED)
        assert product == self.CLASSIFIED
        assert gene == self.CLASSIFIED

    def test_kegg_entry_without_an_ec_still_names_the_function(self):
        name = "Q1234|KO:K00001|EC:|PATH:map00361|haloalkane_dehalogenase|Xanthobacter"
        product, gene = xenobiotic_names(name)
        assert product == "haloalkane_dehalogenase|KO:K00001"
        assert gene == "K00001"

    def test_a_name_with_no_fields_is_returned_as_is(self):
        product, gene = xenobiotic_names("bare_name")
        assert product == "bare_name"
        assert gene == "bare_name"

    def test_none_does_not_raise(self):
        assert xenobiotic_names(None) == (None, None)


# ─────────────────────────────────────────────
# Reusing a previous run's output
# ─────────────────────────────────────────────

class TestStaleOutput:
    """
    The pipeline skips work whose output exists, which is what makes a long run
    resumable and also how a result outlives its inputs. These fix the second
    without giving up the first.
    """

    def test_output_without_a_sidecar_is_never_current(self, tmp_path):
        """An interrupted run leaves a truncated table and no sidecar."""
        out = tmp_path / "blast.out"
        out.write_text("half a table\n")
        src = tmp_path / "query.faa"
        src.write_text(">q\nMKV\n")

        sig = stale_output.signature(inputs=[str(src)])

        assert stale_output.is_current(str(out), sig) is False

    def test_recorded_output_is_reused(self, tmp_path):
        out = tmp_path / "blast.out"
        out.write_text("a table\n")
        src = tmp_path / "query.faa"
        src.write_text(">q\nMKV\n")

        sig = stale_output.signature(inputs=[str(src)])
        stale_output.record(str(out), sig)

        assert stale_output.is_current(str(out), sig) is True

    def test_a_changed_input_invalidates_the_output(self, tmp_path):
        out = tmp_path / "blast.out"
        out.write_text("a table\n")
        db = tmp_path / "refs.fasta"
        db.write_text(">ref\nMKV\n")
        stale_output.record(str(out), stale_output.signature(inputs=[str(db)]))

        db.write_text(">ref\nMKVQQ\n")   # database rebuilt

        assert stale_output.is_current(
            str(out), stale_output.signature(inputs=[str(db)])) is False

    def test_touching_a_file_without_changing_it_keeps_the_output(self, tmp_path):
        """
        A checkout, a copy or an rsync moves mtime without changing a byte.
        Redoing a whole phase for that is its own wrong answer.
        """
        out = tmp_path / "blast.out"
        out.write_text("a table\n")
        db = tmp_path / "refs.fasta"
        db.write_text(">ref\nMKV\n")
        stale_output.record(str(out), stale_output.signature(inputs=[str(db)]))

        os.utime(str(db), (0, 0))

        assert stale_output.is_current(
            str(out), stale_output.signature(inputs=[str(db)])) is True

    def test_a_changed_scalar_invalidates_the_output(self, tmp_path):
        out = tmp_path / "blast.out"
        out.write_text("a table\n")
        stale_output.record(str(out), stale_output.signature(extra={"outfmt": "6 qseqid"}))

        assert stale_output.is_current(
            str(out), stale_output.signature(extra={"outfmt": "6 qseqid stitle"})) is False

    def test_a_missing_input_is_part_of_the_signature(self, tmp_path):
        """A search run while a database was absent must not survive its arrival."""
        out = tmp_path / "blast.out"
        out.write_text("empty\n")
        db = tmp_path / "refs.fasta"
        stale_output.record(str(out), stale_output.signature(inputs=[str(db)]))

        db.write_text(">ref\nMKV\n")

        assert stale_output.is_current(
            str(out), stale_output.signature(inputs=[str(db)])) is False

    def test_discard_removes_the_sidecar_too(self, tmp_path):
        out = tmp_path / "blast.out"
        out.write_text("a table\n")
        stale_output.record(str(out), stale_output.signature())

        stale_output.discard(str(out))

        assert not out.exists()
        assert not os.path.exists(stale_output.signature_path(str(out)))
