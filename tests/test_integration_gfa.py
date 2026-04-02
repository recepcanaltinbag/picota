"""
Integration test: gerçek testNitro.gfa dosyasından cycle tespiti.

Assembly adımını atlar — GFA doğrudan verilir.
Beklenen çıktı cyclesOut.fasta ile karşılaştırılır.
"""

import sys
import os
import time
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'picota'))

from src.cycle_finderv2 import (
    GraphWork,
    cycle_info_optimized,
    cycle_match_based_on_contig_id,
    Cycle,
)
from src.cycle_kmer_hash import filter_cycles_with_kmer

# ─── Dosya yolları ────────────────────────────────────────────────────────────

TEST_DATA = os.path.join(os.path.dirname(__file__), '..', 'picota', 'test_data')
GFA_FILE        = os.path.join(TEST_DATA, 'testNitro.gfa')
EXPECTED_FASTA  = os.path.join(TEST_DATA, 'cyclesOut.fasta')

# Dosyalar yoksa testin tamamını atla
pytestmark = pytest.mark.skipif(
    not os.path.exists(GFA_FILE),
    reason="testNitro.gfa bulunamadı"
)


# ─── Yardımcı ────────────────────────────────────────────────────────────────

def read_fasta_lengths(fasta_path):
    """FASTA'daki tüm sekansların uzunluklarını döner."""
    lengths = []
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    lengths.append(len(''.join(current_seq)))
                    current_seq = []
            else:
                current_seq.append(line)
    if current_seq:
        lengths.append(len(''.join(current_seq)))
    return sorted(lengths)


def run_cycle_pipeline(gfa_path,
                       min_size=3000, max_size=100000,
                       min_comp=1, max_comp=25,
                       k_mer_sim=200, threshold_sim=99,
                       path_limit=15):
    """
    GFA'dan cycle listesi üretir (assembly adımı olmadan).
    Döner: (filtered_cycles, timing_dict)
    """
    timings = {}

    GW = GraphWork()
    GW.find_all_path = False
    GW.path_limit = path_limit

    t0 = time.time()
    node_dict, edge_dict = GW.parse_gfa(gfa_path)
    timings['parse_gfa'] = time.time() - t0

    t1 = time.time()
    genome_graph = GW.generate_genome_graph(node_dict, edge_dict)
    GW.dfs_iterative(genome_graph)
    timings['dfs'] = time.time() - t1

    # Path deduplication
    node_lengths = {k: len(v['Sequence']) for k, v in node_dict.items()}
    unique_paths = []
    for path in GW.paths:
        if cycle_match_based_on_contig_id(path, node_lengths, unique_paths):
            unique_paths.append(path)

    # Cycle sekansı oluşturma + boyut/component filtresi
    t2 = time.time()
    cycle_info_list = []
    for path in unique_paths:
        obj = cycle_info_optimized(path, node_dict, edge_dict, cycle_info_list)
        if obj is None or obj == 'Pass':
            continue
        if (min_size <= obj.length <= max_size and
                min_comp <= obj.component_number <= max_comp):
            obj.name = f'Cycle_{len(cycle_info_list) + 1}'
            cycle_info_list.append(obj)
    timings['cycle_build'] = time.time() - t2

    # K-mer benzerlik filtresi
    t3 = time.time()
    filtered = filter_cycles_with_kmer(cycle_info_list, k_mer_sim, threshold_sim, 'Cycle')
    timings['kmer_filter'] = time.time() - t3

    return filtered, timings


# ─── Fixture ─────────────────────────────────────────────────────────────────

@pytest.fixture(scope='module')
def pipeline_result():
    """Tüm modül boyunca tek seferinde çalışır."""
    cycles, timings = run_cycle_pipeline(GFA_FILE)
    return cycles, timings


# ─── Testler ─────────────────────────────────────────────────────────────────

class TestGfaParsing:
    def test_node_count(self):
        GW = GraphWork()
        nodes, edges = GW.parse_gfa(GFA_FILE)
        # testNitro.gfa: 272 S satırı → 544 node (forward + reverse)
        assert len(nodes) == 544

    def test_edge_count(self):
        GW = GraphWork()
        nodes, edges = GW.parse_gfa(GFA_FILE)
        # 311 L satırı → 622 edge (forward + reverse)
        assert len(edges) == 622

    def test_nodes_have_both_strands(self):
        GW = GraphWork()
        nodes, _ = GW.parse_gfa(GFA_FILE)
        for key in list(nodes.keys())[:10]:
            node_id = key[:-1]
            assert node_id + '+' in nodes
            assert node_id + '-' in nodes

    def test_reverse_complement_node(self):
        from src.cycle_finderv2 import reverse_complement
        GW = GraphWork()
        nodes, _ = GW.parse_gfa(GFA_FILE)
        for key in list(nodes.keys())[:5]:
            if key.endswith('+'):
                node_id = key[:-1]
                fwd_seq = nodes[node_id + '+']['Sequence']
                rev_seq = nodes[node_id + '-']['Sequence']
                assert rev_seq == reverse_complement(fwd_seq)

    def test_edges_have_overlap(self):
        GW = GraphWork()
        _, edges = GW.parse_gfa(GFA_FILE)
        for edge_data in list(edges.values())[:20]:
            assert 'Overlap' in edge_data
            assert edge_data['Overlap'] != ''


class TestCycleDetection:
    def test_cycles_found(self, pipeline_result):
        cycles, _ = pipeline_result
        assert len(cycles) > 0

    def test_cycle_count_reasonable(self, pipeline_result):
        # testNitro.gfa'dan 5-15 arası cycle beklenir
        cycles, _ = pipeline_result
        assert 5 <= len(cycles) <= 15

    def test_all_cycles_within_size_bounds(self, pipeline_result):
        cycles, _ = pipeline_result
        for c in cycles:
            assert 3000 <= c.length <= 100000, \
                f"{c.name}: {c.length}bp boyut sınırları dışında"

    def test_all_cycles_within_component_bounds(self, pipeline_result):
        cycles, _ = pipeline_result
        for c in cycles:
            assert 1 <= c.component_number <= 25, \
                f"{c.name}: {c.component_number} component sınırları dışında"

    def test_no_duplicate_names(self, pipeline_result):
        cycles, _ = pipeline_result
        names = [c.name for c in cycles]
        assert len(names) == len(set(names))

    def test_sequences_are_dna(self, pipeline_result):
        cycles, _ = pipeline_result
        valid = set('ATCGatcg')
        for c in cycles:
            invalid = set(c.sequence) - valid
            assert not invalid, \
                f"{c.name}: Geçersiz karakterler: {invalid}"

    def test_sequence_length_matches_attribute(self, pipeline_result):
        cycles, _ = pipeline_result
        for c in cycles:
            assert len(c.sequence) == c.length, \
                f"{c.name}: sequence uzunluğu {len(c.sequence)} ama length={c.length}"


class TestExpectedOutput:
    """cyclesOut.fasta ile karşılaştırma."""

    @pytest.mark.skipif(
        not os.path.exists(EXPECTED_FASTA),
        reason="cyclesOut.fasta bulunamadı"
    )
    def test_known_cycle_lengths_found(self, pipeline_result):
        """
        cyclesOut.fasta'daki bilinen uzunluklar bulunan cycle'larda mevcut olmalı.
        Birebir eşleşme değil — aynı uzunlukta cycle var mı kontrolü.
        """
        cycles, _ = pipeline_result
        expected_lengths = set(read_fasta_lengths(EXPECTED_FASTA))
        found_lengths = set(c.length for c in cycles)

        # En az bir beklenen uzunluk bulunmalı
        matches = expected_lengths & found_lengths
        assert len(matches) > 0, \
            f"Beklenen uzunluklardan hiçbiri bulunamadı.\n" \
            f"Beklenen: {sorted(expected_lengths)}\n" \
            f"Bulunan:  {sorted(found_lengths)}"

    @pytest.mark.skipif(
        not os.path.exists(EXPECTED_FASTA),
        reason="cyclesOut.fasta bulunamadı"
    )
    def test_cycle_count_vs_expected(self, pipeline_result):
        """
        Bulunan cycle sayısı beklenen sayının ±50% içinde olmalı.
        Tam eşleşme değil — kmer threshold ve path_limit'e bağlı.
        """
        cycles, _ = pipeline_result
        expected_count = len(read_fasta_lengths(EXPECTED_FASTA))
        found_count = len(cycles)
        assert abs(found_count - expected_count) <= max(expected_count // 2, 3), \
            f"Cycle sayısı çok farklı: beklenen {expected_count}, bulunan {found_count}"


class TestPerformance:
    """Pipeline'ın makul sürede çalışmasını kontrol eder."""

    def test_parse_gfa_fast(self, pipeline_result):
        _, timings = pipeline_result
        # 272 node'lu GFA 1 saniyeden kısa parse edilmeli
        assert timings['parse_gfa'] < 1.0, \
            f"parse_gfa çok yavaş: {timings['parse_gfa']:.2f}s"

    def test_dfs_fast(self, pipeline_result):
        _, timings = pipeline_result
        # DFS 1 saniyeden kısa bitmeli
        assert timings['dfs'] < 1.0, \
            f"DFS çok yavaş: {timings['dfs']:.2f}s"

    def test_cycle_build_fast(self, pipeline_result):
        _, timings = pipeline_result
        # Cycle sekans oluşturma 1 saniyeden kısa
        assert timings['cycle_build'] < 1.0, \
            f"cycle_build çok yavaş: {timings['cycle_build']:.2f}s"

    def test_total_pipeline_under_5s(self):
        """Assembly dahil değil, sadece GFA→Cycles."""
        start = time.time()
        run_cycle_pipeline(GFA_FILE)
        elapsed = time.time() - start
        assert elapsed < 5.0, \
            f"Pipeline çok yavaş: {elapsed:.2f}s (5s limit)"

    def test_kmer_filter_bottleneck_documented(self, pipeline_result):
        """
        k-mer filtresi şu an O(n²) — bu test bottleneck'i belgeler.
        10 cycle için 1s'den az olmalı; daha fazla cycle'da yavaşlayacak.
        """
        cycles, timings = pipeline_result
        n = len(cycles)
        assert timings['kmer_filter'] < 2.0, \
            f"kmer_filter {n} cycle için {timings['kmer_filter']:.2f}s aldı"
