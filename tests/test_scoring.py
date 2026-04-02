"""
Unit tests for scoringv4ProtBlast.py

Tested fonksiyonlar:
  - CodingRegion class
  - GeneticInfo class
  - calculate_total_score()
  - merge_intervals()
  - parsing_blast_file()
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
)


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
