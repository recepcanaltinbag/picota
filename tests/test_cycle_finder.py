"""
Unit tests for cycle_finderv2.py

Tested fonksiyonlar:
  - reverse_complement()
  - find_overlap_length()
  - is_similar()
  - is_similar_polynomial()
  - cycle_match_based_on_contig_id()
  - cycle_info()
  - cycle_info_optimized()
  - Graph class
  - GraphWork.parse_gfa()
  - GraphWork.dfs_iterative()
"""

import sys
import os
import tempfile
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.cycle_finderv2 import (
    reverse_complement,
    find_overlap_length,
    is_similar,
    is_similar_polynomial,
    cycle_match_based_on_contig_id,
    cycle_info,
    cycle_info_optimized,
    Graph,
    GraphWork,
    Cycle,
)


# ─────────────────────────────────────────────
# reverse_complement
# ─────────────────────────────────────────────

class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_all_adenine(self):
        assert reverse_complement("AAAA") == "TTTT"

    def test_palindrome(self):
        # ATAT reverse complement = ATAT
        assert reverse_complement("ATAT") == "ATAT"

    def test_lowercase(self):
        assert reverse_complement("atcg") == "cgat"

    def test_mixed_case(self):
        assert reverse_complement("AtCg") == "cGaT"

    def test_single_nucleotide(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("G") == "C"
        assert reverse_complement("C") == "G"

    def test_empty_string(self):
        assert reverse_complement("") == ""

    def test_longer_sequence(self):
        seq = "ATCGATCGATCG"
        rc = reverse_complement(seq)
        assert len(rc) == len(seq)
        # rc'nin rc'si orijinal olmalı
        assert reverse_complement(rc) == seq

    def test_invalid_character_raises(self):
        # 'N' complement_dict içinde yok → KeyError beklenir
        with pytest.raises(KeyError):
            reverse_complement("ATCGN")


# ─────────────────────────────────────────────
# find_overlap_length
# ─────────────────────────────────────────────

class TestFindOverlapLength:
    def test_zero_overlap(self):
        assert find_overlap_length("0M") == 0

    def test_small_overlap(self):
        assert find_overlap_length("5M") == 5

    def test_large_overlap(self):
        assert find_overlap_length("127M") == 127

    def test_single_digit(self):
        assert find_overlap_length("1M") == 1

    def test_returns_int(self):
        result = find_overlap_length("12M")
        assert isinstance(result, int)


# ─────────────────────────────────────────────
# is_similar
# ─────────────────────────────────────────────

class TestIsSimilar:
    def test_identical_sequences(self):
        seq = "ATCGATCGATCGATCGATCGATCGATCG"
        assert is_similar(seq, seq, k_mer_sim=5, threshold_sim=95) is True

    def test_completely_different(self):
        # "AAAA..." RC'si "TTTT..." olduğu için is_similar True döner — bu biyolojik olarak doğru
        seq1 = "AAAAAAAAAAAAAAAAAAAAAAAAA"
        seq2 = "TTTTTTTTTTTTTTTTTTTTTTTTT"
        assert is_similar(seq1, seq2, k_mer_sim=5, threshold_sim=95) is True

    def test_reverse_complement_match(self):
        seq = "ATCGATCGATCGATCG"
        rc = reverse_complement(seq)
        assert is_similar(seq, rc, k_mer_sim=4, threshold_sim=90) is True

    def test_threshold_at_zero(self):
        # is_similar strict olarak > kullanıyor (>= değil)
        # Eşleşen k-mer yoksa (exact_count=0): 0 > 0 False döner
        seq1 = "ATCGATCG"
        seq2 = "TTTTTTTT"
        # "TTTT" k-mer'i "ATCGATCG" veya RC'sinde yok → False
        assert is_similar(seq1, seq2, k_mer_sim=4, threshold_sim=0) is False

    def test_partial_overlap(self):
        # İlk yarısı aynı, ikinci yarısı farklı
        seq1 = "ATCGATCGATCGATCGATCG"
        seq2 = "ATCGATCGATCGATCGTTTT"
        # ~%75 similarity beklenir
        result_low = is_similar(seq1, seq2, k_mer_sim=4, threshold_sim=50)
        result_high = is_similar(seq1, seq2, k_mer_sim=4, threshold_sim=99)
        assert result_low is True
        assert result_high is False

    def test_seq2_shorter_than_kmer(self):
        # seq2 k_mer_sim'den kısa → C_t boş → ZeroDivisionError riski
        with pytest.raises((ZeroDivisionError, ValueError)):
            is_similar("ATCGATCG", "AT", k_mer_sim=5, threshold_sim=95)


# ─────────────────────────────────────────────
# is_similar_polynomial
# ─────────────────────────────────────────────

class TestIsSimilarPolynomial:
    def test_identical(self):
        seq = "ATCGATCGATCGATCG"
        assert is_similar_polynomial(seq, seq, k_mer_sim=4, threshold_sim=95) is True

    def test_no_match(self):
        # RC(AAAA...) = TTTT... → is_similar_polynomial True döner (RC kontrolü yapıyor)
        seq1 = "AAAAAAAAAAAAAAAA"
        seq2 = "TTTTTTTTTTTTTTTT"
        assert is_similar_polynomial(seq1, seq2, k_mer_sim=4, threshold_sim=50) is True

    def test_seq2_shorter_than_kmer(self):
        # total_kmers = 2-5+1 = -2 → range(-2) boş → exact_count=0, 0/-2=0 → False döner, hata vermez
        result = is_similar_polynomial("ATCGATCG", "AT", k_mer_sim=5, threshold_sim=95)
        assert result is False


# ─────────────────────────────────────────────
# cycle_match_based_on_contig_id
# ─────────────────────────────────────────────

class TestCycleMatchBasedOnContigId:
    def setup_method(self):
        self.nodes_len = {
            "10+": 100, "10-": 100,
            "20+": 200, "20-": 200,
            "30+": 150, "30-": 150,
        }

    def test_empty_new_parts_returns_true(self):
        path = ["10+", "20+"]
        assert cycle_match_based_on_contig_id(path, self.nodes_len, []) is True

    def test_identical_path_returns_false(self):
        path = ["10+", "20+"]
        new_parts = [["10+", "20+"]]
        assert cycle_match_based_on_contig_id(path, self.nodes_len, new_parts) is False

    def test_completely_different_path_returns_true(self):
        path = ["30+"]
        new_parts = [["10+", "20+"]]
        assert cycle_match_based_on_contig_id(path, self.nodes_len, new_parts) is True

    def test_partial_overlap_below_threshold(self):
        # 10+ paylaşılıyor (100bp), path total = 300bp → %33 < 70 threshold
        path = ["10+", "20+", "30+"]
        new_parts = [["10+"]]
        assert cycle_match_based_on_contig_id(path, self.nodes_len, new_parts) is True

    def test_partial_overlap_above_threshold(self):
        # path=["10+","20+"] (300bp), existing=["10+","20+","30+"] (450bp)
        # common=300, max=450 → identity=66.7% < 70 threshold → True (farklı sayılır)
        path = ["10+", "20+"]
        new_parts = [["10+", "20+", "30+"]]
        assert cycle_match_based_on_contig_id(path, self.nodes_len, new_parts) is True

    def test_custom_threshold(self):
        path = ["10+", "20+"]
        new_parts = [["10+"]]
        # 100/300 = %33 → threshold=20'nin üstünde → False
        assert cycle_match_based_on_contig_id(path, self.nodes_len, new_parts, threshold=20) is False
        # threshold=50 → False olmaz
        assert cycle_match_based_on_contig_id(path, self.nodes_len, new_parts, threshold=50) is True


# ─────────────────────────────────────────────
# Graph class
# ─────────────────────────────────────────────

class TestGraph:
    def test_basic_construction(self):
        edges = [("A+", "B+"), ("B+", "C+")]
        g = Graph(edges)
        assert "A+" in g.adj
        assert "B+" in g.adj["A+"]

    def test_empty_edges(self):
        g = Graph([])
        assert len(g.adj) == 0

    def test_self_loop(self):
        g = Graph([("A+", "A+")])
        assert "A+" in g.adj["A+"]

    def test_adjacency_list_correct(self):
        edges = [("1+", "2+"), ("1+", "3-"), ("2+", "3-")]
        g = Graph(edges)
        assert set(g.adj["1+"]) == {"2+", "3-"}
        assert g.adj["2+"] == ["3-"]


# ─────────────────────────────────────────────
# cycle_info & cycle_info_optimized
# ─────────────────────────────────────────────

def make_nodes_and_edges(path_seqs, overlaps):
    """
    path_seqs: list of (node_id, sequence)
    overlaps:  list of overlap strings between consecutive nodes (len = len(path_seqs) - 1)
    """
    nodes = {}
    for node_id, seq in path_seqs:
        nodes[node_id] = {"Name": node_id, "Sequence": seq}

    edges = {}
    for i in range(len(overlaps)):
        n1 = path_seqs[i][0]
        n2 = path_seqs[i + 1][0]
        edges[(n1, n2)] = {"Overlap": overlaps[i]}

    path = [n for n, _ in path_seqs]
    return path, nodes, edges


class TestCycleInfo:
    def test_two_node_cycle_no_overlap(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCGATCG"), ("2+", "TTTTGGGG")],
            ["0M"]
        )
        result = cycle_info(path, nodes, edges, [])
        assert isinstance(result, Cycle)
        assert result.sequence == "ATCGATCGTTTTGGGG"
        assert result.component_number == 2

    def test_two_node_cycle_with_overlap(self):
        # 4bp overlap: "ATCGATCG" + "ATCGTTTT" overlap=4M → "ATCG" cut from first
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCGATCG"), ("2+", "ATCGTTTT")],
            ["4M"]
        )
        result = cycle_info(path, nodes, edges, [])
        assert isinstance(result, Cycle)
        # İlk düğümden son 4bp kesilir: "ATCG" + ikinci düğüm tam = "ATCGATCGTTTT"
        assert result.sequence == "ATCGATCGTTTT"

    def test_star_overlap_returns_none(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCG"), ("2+", "TTTT")],
            ["*"]
        )
        result = cycle_info(path, nodes, edges, [])
        assert result is None

    def test_duplicate_returns_pass(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCGATCG"), ("2+", "TTTTGGGG")],
            ["0M"]
        )
        existing = Cycle("c1", "ATCGATCGTTTTGGGG", 16, 2, path)
        result = cycle_info(path, nodes, edges, [existing])
        assert result == 'Pass'

    def test_reverse_complement_duplicate_returns_pass(self):
        seq = "ATCGATCGTTTTGGGG"
        rc = reverse_complement(seq)
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCGATCG"), ("2+", "TTTTGGGG")],
            ["0M"]
        )
        existing = Cycle("c1", rc, len(rc), 2, path)
        result = cycle_info(path, nodes, edges, [existing])
        assert result == 'Pass'

    def test_three_node_cycle(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "AAAA"), ("2+", "CCCC"), ("3+", "GGGG")],
            ["0M", "0M"]
        )
        result = cycle_info(path, nodes, edges, [])
        assert isinstance(result, Cycle)
        assert result.component_number == 3
        assert result.sequence == "AAAACCCCGGGG"

    def test_empty_cycle_list_allowed(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCG"), ("2+", "GCTA")],
            ["0M"]
        )
        result = cycle_info(path, nodes, edges, [])
        assert isinstance(result, Cycle)


class TestCycleInfoOptimized:
    """cycle_info_optimized aynı davranışı göstermeli."""

    def test_two_node_no_overlap(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCGATCG"), ("2+", "TTTTGGGG")],
            ["0M"]
        )
        result = cycle_info_optimized(path, nodes, edges, [])
        assert isinstance(result, Cycle)
        assert result.sequence == "ATCGATCGTTTTGGGG"

    def test_star_overlap_returns_none(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCG"), ("2+", "TTTT")],
            ["*"]
        )
        result = cycle_info_optimized(path, nodes, edges, [])
        assert result is None

    def test_duplicate_returns_pass(self):
        path, nodes, edges = make_nodes_and_edges(
            [("1+", "ATCGATCG"), ("2+", "TTTTGGGG")],
            ["0M"]
        )
        existing = Cycle("c1", "ATCGATCGTTTTGGGG", 16, 2, path)
        result = cycle_info_optimized(path, nodes, edges, [existing])
        assert result == 'Pass'

    def test_consistent_with_cycle_info(self):
        """İki fonksiyon aynı sequence üretmeli."""
        path, nodes, edges = make_nodes_and_edges(
            [("A+", "ATCGATCGATCG"), ("B+", "GCTAGCTAGCTA"), ("C+", "TTTTAAAACCCC")],
            ["3M", "3M"]
        )
        r1 = cycle_info(path, nodes, edges, [])
        r2 = cycle_info_optimized(path, nodes, edges, [])
        assert isinstance(r1, Cycle)
        assert isinstance(r2, Cycle)
        assert r1.sequence == r2.sequence


# ─────────────────────────────────────────────
# GraphWork.parse_gfa
# ─────────────────────────────────────────────

MINIMAL_GFA = """\
H\tVN:Z:1.0
S\t1\tATCGATCG
S\t2\tTTTTGGGG
L\t1\t+\t2\t-\t4M
"""

class TestParseGfa:
    def _write_gfa(self, content):
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.gfa', delete=False)
        f.write(content)
        f.flush()
        return f.name

    def test_nodes_created(self):
        gfa_path = self._write_gfa(MINIMAL_GFA)
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(gfa_path)
        os.unlink(gfa_path)
        assert "1+" in nodes
        assert "1-" in nodes
        assert "2+" in nodes
        assert "2-" in nodes

    def test_reverse_complement_node(self):
        gfa_path = self._write_gfa(MINIMAL_GFA)
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(gfa_path)
        os.unlink(gfa_path)
        assert nodes["1-"]["Sequence"] == reverse_complement("ATCGATCG")

    def test_edges_created(self):
        gfa_path = self._write_gfa(MINIMAL_GFA)
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(gfa_path)
        os.unlink(gfa_path)
        # Forward ve reverse edge oluşturulmalı
        assert ("1+", "2-") in edges
        assert ("2+", "1-") in edges  # reverse link

    def test_edge_overlap(self):
        gfa_path = self._write_gfa(MINIMAL_GFA)
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(gfa_path)
        os.unlink(gfa_path)
        assert edges[("1+", "2-")]["Overlap"] == "4M"

    def test_empty_gfa(self):
        gfa_path = self._write_gfa("H\tVN:Z:1.0\n")
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(gfa_path)
        os.unlink(gfa_path)
        assert len(nodes) == 0
        assert len(edges) == 0

    def test_multiple_nodes(self):
        content = "H\tVN:Z:1.0\n"
        for i in range(1, 6):
            content += f"S\t{i}\tATCGATCG\n"
        gfa_path = self._write_gfa(content)
        gw = GraphWork()
        nodes, edges = gw.parse_gfa(gfa_path)
        os.unlink(gfa_path)
        assert len(nodes) == 10  # 5 forward + 5 reverse


# ─────────────────────────────────────────────
# GraphWork.dfs_iterative — cycle detection
# ─────────────────────────────────────────────

class TestDfsIterative:
    def test_simple_cycle_detected(self):
        # A+ → B+ → A+ döngüsü
        edges = [("A+", "B+"), ("B+", "A+")]
        g = Graph(edges)
        gw = GraphWork()
        gw.dfs_iterative(g)
        assert len(gw.cycles) > 0

    def test_no_cycle_in_dag(self):
        # Yönlü asiklik graf
        edges = [("A+", "B+"), ("B+", "C+")]
        g = Graph(edges)
        gw = GraphWork()
        gw.dfs_iterative(g)
        assert gw.cycles == []
        assert gw.reverse_or_cycles == []

    def test_self_loop_detected(self):
        edges = [("A+", "A+")]
        g = Graph(edges)
        gw = GraphWork()
        gw.dfs_iterative(g)
        assert len(gw.cycles) > 0

    def test_reverse_oriented_cycle_detected(self):
        # A+ → B- → A- (ters yönelimli döngü)
        edges = [("A+", "B-"), ("B-", "A-")]
        g = Graph(edges)
        gw = GraphWork()
        gw.dfs_iterative(g)
        assert len(gw.reverse_or_cycles) > 0

    def test_multiple_cycles(self):
        edges = [
            ("A+", "B+"), ("B+", "A+"),   # döngü 1
            ("C+", "D+"), ("D+", "C+"),   # döngü 2
        ]
        g = Graph(edges)
        gw = GraphWork()
        gw.dfs_iterative(g)
        assert len(gw.cycles) >= 2

    def test_empty_graph(self):
        g = Graph([])
        gw = GraphWork()
        gw.dfs_iterative(g)
        assert gw.cycles == []
