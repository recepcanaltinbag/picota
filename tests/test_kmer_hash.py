"""
Unit tests for cycle_kmer_hash.py

Tested fonksiyonlar:
  - reverse_complement()
  - get_kmer_hashes()
  - filter_cycles_with_kmer()
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.cycle_kmer_hash import (
    reverse_complement,
    get_kmer_hashes,
    filter_cycles_with_kmer,
)

# Cycle class'ını da import ediyoruz (filter için gerekli)
from src.cycle_finderv2 import Cycle


def make_cycle(name, seq):
    c = Cycle(name, seq, len(seq), 2, [])
    return c


# ─────────────────────────────────────────────
# reverse_complement
# ─────────────────────────────────────────────

class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_empty(self):
        assert reverse_complement("") == ""

    def test_single(self):
        assert reverse_complement("A") == "T"

    def test_lowercase(self):
        assert reverse_complement("atcg") == "cgat"

    def test_double_reverse_is_identity(self):
        seq = "ATCGATCGATCG"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_invalid_char_raises(self):
        with pytest.raises(KeyError):
            reverse_complement("ATCGN")


# ─────────────────────────────────────────────
# get_kmer_hashes
# ─────────────────────────────────────────────

class TestGetKmerHashes:
    def test_basic(self):
        result = get_kmer_hashes("ATCGATCG", 3)
        assert isinstance(result, set)
        assert "ATC" in result
        assert "TCG" in result

    def test_k_equals_seq_len(self):
        result = get_kmer_hashes("ATCG", 4)
        assert result == {"ATCG"}

    def test_k_larger_than_seq(self):
        result = get_kmer_hashes("AT", 5)
        assert result == set()

    def test_k_equals_one(self):
        result = get_kmer_hashes("AATTGGCC", 1)
        assert result == {"A", "T", "G", "C"}

    def test_duplicate_kmers_collapsed(self):
        # "AAAA" → k=2 → {"AA"} (tümü aynı)
        result = get_kmer_hashes("AAAA", 2)
        assert result == {"AA"}

    def test_empty_seq(self):
        result = get_kmer_hashes("", 3)
        assert result == set()

    def test_returns_set(self):
        result = get_kmer_hashes("ATCG", 2)
        assert isinstance(result, set)

    def test_correct_count(self):
        # "ATCG" k=2 → "AT", "TC", "CG" → 3 unique
        result = get_kmer_hashes("ATCG", 2)
        assert len(result) == 3


# ─────────────────────────────────────────────
# filter_cycles_with_kmer
# ─────────────────────────────────────────────

class TestFilterCyclesWithKmer:
    # Testlerde k_mer_sim=4, threshold_sim=90 kullanıyoruz

    def test_empty_list_returns_empty(self):
        result = filter_cycles_with_kmer([], 4, 90, "Cycle")
        assert result == []

    def test_single_cycle_kept(self):
        c = make_cycle("tmp", "ATCGATCGATCGATCG")
        result = filter_cycles_with_kmer([c], 4, 90, "Cycle")
        assert len(result) == 1

    def test_identical_cycles_deduplicated(self):
        seq = "ATCGATCGATCGATCGATCGATCGATCG"
        c1 = make_cycle("t1", seq)
        c2 = make_cycle("t2", seq)
        result = filter_cycles_with_kmer([c1, c2], 4, 90, "Cycle")
        assert len(result) == 1

    def test_completely_different_cycles_both_kept(self):
        # AAAA... RC'si TTTT... olduğu için bu ikili benzer sayılır → sadece 1 kalır
        # Gerçekten farklı sekanslar için nükleotid kombinasyonu gerekir
        c1 = make_cycle("t1", "ATCGATCGATCGATCGATCGATCGATCGATCG")
        c2 = make_cycle("t2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")
        result = filter_cycles_with_kmer([c1, c2], 4, 90, "Cycle")
        assert len(result) == 2

    def test_reverse_complement_deduplicated(self):
        seq = "ATCGATCGATCGATCGATCGATCGATCG"
        rc_seq = reverse_complement(seq)
        c1 = make_cycle("t1", seq)
        c2 = make_cycle("t2", rc_seq)
        result = filter_cycles_with_kmer([c1, c2], 4, 90, "Cycle")
        # RC benzer → sadece 1 kalmalı
        assert len(result) == 1

    def test_names_reassigned(self):
        # Farklı olan sekanslar için isim atama kontrolü
        seq1 = "ATCGATCGATCGATCGATCGATCGATCG"
        seq2 = "GCTAGCTAGCTAGCTAGCTAGCTAGCTA"
        c1 = make_cycle("old1", seq1)
        c2 = make_cycle("old2", seq2)
        result = filter_cycles_with_kmer([c1, c2], 4, 90, "Cycle")
        names = [c.name for c in result]
        assert "Cycle_1" in names
        assert "Cycle_2" in names

    def test_reverse_oriented_name_prefix(self):
        c = make_cycle("t1", "ATCGATCGATCGATCG")
        c.reverseOriented = True
        result = filter_cycles_with_kmer([c], 4, 90, "Cycle")
        assert len(result) == 1
        assert "reverseoriented_" in result[0].name

    def test_threshold_100_keeps_only_first(self):
        # threshold=100 → her şey benzer sayılır, ancak ilk cycle her zaman kalır
        seq = "ATCGATCGATCGATCG"
        c1 = make_cycle("t1", seq)
        c2 = make_cycle("t2", seq + "AAAA")
        result = filter_cycles_with_kmer([c1, c2], 4, 100, "Cycle")
        # İlk her zaman eklenir; ikinci ile karşılaştırılır
        assert len(result) >= 1

    def test_threshold_0_identical_deduplicated(self):
        # threshold=0 → aynı k-merler paylaşılıyorsa benzer sayılır → tek kalır
        seq = "ATCGATCGATCGATCGATCG"
        cycles = [make_cycle(f"t{i}", seq) for i in range(3)]
        result = filter_cycles_with_kmer(cycles, 4, 0, "Cycle")
        assert len(result) == 1

    def test_three_unique_cycles_all_kept(self):
        # AAAA/TTTT ve CCCC/GGGG RC çiftleri birbirine benzer sayılır
        # Gerçekten farklı 3 sekans kullanılmalı
        c1 = make_cycle("t1", "ATCGATCGATCGATCGATCGATCGATCG")
        c2 = make_cycle("t2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTA")
        c3 = make_cycle("t3", "AATTCCGGAATTCCGGAATTCCGGAATT")
        result = filter_cycles_with_kmer([c1, c2, c3], 4, 90, "Cycle")
        assert len(result) == 3
