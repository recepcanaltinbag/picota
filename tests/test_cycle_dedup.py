"""
Unit tests for src/cycle_dedup.py (roadmap phase 2).

The end-to-end effect of this module is covered in test_known_defects.py; here
the individual measures are pinned, because the whole point of the rewrite is
that they have properties the legacy ones lacked -- symmetry, multiplicity
awareness, and invariance under how a cycle happens to be written down.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "picota"))
sys.path.insert(0, os.path.dirname(__file__))

from src.cycle_dedup import (  # noqa: E402
    canonical_kmers,
    canonical_path_key,
    dedup_paths,
    filter_cycles_multiset,
    flip_orientation,
    multiset_jaccard,
    paths_are_duplicates,
    reverse_complement,
    reverse_traversal,
)
from src.cycle_finderv2 import Cycle  # noqa: E402
from synthetic_gfa import make_dna  # noqa: E402


class TestOrientation:
    def test_flip_plus_to_minus(self):
        assert flip_orientation("14349+") == "14349-"

    def test_flip_minus_to_plus(self):
        assert flip_orientation("14349-") == "14349+"

    def test_unsigned_node_unchanged(self):
        assert flip_orientation("14349") == "14349"

    def test_reverse_traversal_reverses_and_flips(self):
        assert reverse_traversal(["a+", "b+", "c-"]) == ["c+", "b-", "a-"]

    def test_reverse_traversal_is_an_involution(self):
        path = ["a+", "b-", "c+"]
        assert reverse_traversal(reverse_traversal(path)) == path


class TestCanonicalPathKey:
    def test_rotations_share_a_key(self):
        """A cycle has no start, so where the path was cut must not matter."""
        base = ["a+", "b+", "c+"]
        for offset in range(3):
            rotated = base[offset:] + base[:offset]
            assert canonical_path_key(rotated) == canonical_path_key(base)

    def test_reverse_traversal_shares_a_key(self):
        """
        The same physical cycle read the other way round is not a new candidate.
        This is what legacy got wrong: its similarity test compared
        orientation-tagged names, so a-/b- looked unrelated to a+/b+.
        """
        path = ["repeat+", "cargo+"]
        assert canonical_path_key(reverse_traversal(path)) == canonical_path_key(path)

    def test_bubbles_sharing_a_repeat_have_different_keys(self):
        """The defining case: shared repeat, different cargo, two real CTs."""
        a = canonical_path_key(["repeat+", "cargoA+"])
        b = canonical_path_key(["repeat+", "cargoB+"])
        assert a != b

    def test_same_nodes_in_different_cyclic_order_differ(self):
        assert canonical_path_key(["a+", "b+", "c+"]) != canonical_path_key(["a+", "c+", "b+"])

    def test_empty_path(self):
        assert canonical_path_key([]) == ()

    def test_paths_are_duplicates_wraps_the_key(self):
        assert paths_are_duplicates(["a+", "b+"], ["b+", "a+"]) is True
        assert paths_are_duplicates(["a+", "b+"], ["a+", "c+"]) is False


class TestDedupPaths:
    def test_keeps_every_bubble_through_a_shared_repeat(self):
        paths = [["repeat+", f"cargo{i}+"] for i in range(5)]
        assert len(dedup_paths(paths)) == 5

    def test_collapses_rotations_and_reverse_traversals(self):
        paths = [["a+", "b+", "c+"],
                 ["b+", "c+", "a+"],          # rotation
                 ["c-", "b-", "a-"]]          # reverse traversal
        assert dedup_paths(paths) == [["a+", "b+", "c+"]]

    def test_input_order_is_preserved(self):
        paths = [["b+", "x+"], ["a+", "y+"]]
        assert dedup_paths(paths) == paths

    def test_empty_input(self):
        assert dedup_paths([]) == []


class TestCanonicalKmers:
    def test_strand_folding_makes_counts_identical(self):
        seq = make_dna(200, 1)
        assert canonical_kmers(seq, 21) == canonical_kmers(reverse_complement(seq), 21)

    def test_rotation_of_a_cycle_gives_identical_counts(self):
        """Circular extraction is why two rotations of one cycle compare equal."""
        seq = make_dna(200, 2)
        rotated = seq[50:] + seq[:50]
        assert canonical_kmers(seq, 21) == canonical_kmers(rotated, 21)

    def test_repeat_is_counted_twice(self):
        """The multiset property D3 was missing."""
        unit = make_dna(150, 3)
        once = canonical_kmers(unit, 21)
        twice = canonical_kmers(unit + unit, 21)
        assert sum(twice.values()) == pytest.approx(2 * sum(once.values()), rel=0.2)
        assert max(twice.values()) >= 2

    def test_sequence_shorter_than_k_yields_nothing(self):
        assert canonical_kmers("ACGT", 21) == {}

    def test_non_positive_k_yields_nothing(self):
        assert canonical_kmers(make_dna(100, 4), 0) == {}


class TestMultisetJaccard:
    def test_identical_sequences_score_one(self):
        counts = canonical_kmers(make_dna(500, 5), 21)
        assert multiset_jaccard(counts, counts) == pytest.approx(1.0)

    def test_symmetric(self):
        """D4: the outcome must not depend on which came first."""
        a = canonical_kmers(make_dna(500, 6), 21)
        b = canonical_kmers(make_dna(700, 7), 21)
        assert multiset_jaccard(a, b) == multiset_jaccard(b, a)

    def test_containment_does_not_score_one(self):
        """
        A short cycle fully inside a long one scored ~100% under the legacy
        containment measure and was deleted. Jaccard reflects the size gap.
        """
        shared = make_dna(1500, 8)
        small = canonical_kmers(shared, 21)
        large = canonical_kmers(shared + make_dna(1500, 9), 21)
        assert multiset_jaccard(small, large) == pytest.approx(0.5, abs=0.05)

    def test_repeat_count_lowers_the_score(self):
        unit = make_dna(1500, 10)
        once = canonical_kmers(unit, 21)
        twice = canonical_kmers(unit + unit, 21)
        assert multiset_jaccard(once, twice) < 0.75

    def test_disjoint_sequences_score_zero(self):
        a = canonical_kmers(make_dna(500, 11), 31)
        b = canonical_kmers(make_dna(500, 12), 31)
        assert multiset_jaccard(a, b) == pytest.approx(0.0, abs=0.01)

    def test_empty_input_scores_zero(self):
        assert multiset_jaccard({}, canonical_kmers(make_dna(200, 13), 21)) == 0.0


class TestFilterCyclesMultiset:
    def _cycle(self, name, seq, reverse=False):
        cycle = Cycle(name, seq, len(seq), 2, [])
        cycle.reverseOriented = reverse
        return cycle

    def test_empty_list_returned_unchanged(self):
        assert filter_cycles_multiset([], 80, 80, "Cycle") == []

    def test_survivors_are_renamed_sequentially(self):
        cycles = [self._cycle("x", make_dna(4000, 20)),
                  self._cycle("y", make_dna(4000, 21))]
        kept = filter_cycles_multiset(cycles, 80, 80, "Cycle")
        assert [c.name for c in kept] == ["Cycle_1", "Cycle_2"]

    def test_reverse_oriented_cycles_keep_their_name_marker(self):
        """Downstream parsing relies on this prefix, so it must not change."""
        cycles = [self._cycle("x", make_dna(4000, 22), reverse=True)]
        assert filter_cycles_multiset(cycles, 80, 80, "Cycle")[0].name == \
            "Cycle_reverseoriented_1"

    def test_cycles_shorter_than_k_are_kept(self):
        """No k-mers means no evidence of duplication, so keep the candidate."""
        cycles = [self._cycle("x", "ACGT"), self._cycle("y", "ACGT")]
        assert len(filter_cycles_multiset(cycles, 80, 80, "Cycle")) == 2

    def test_threshold_is_a_percentage(self):
        """Lowering the threshold merges candidates that a higher one separates."""
        shared = make_dna(3000, 23)
        cycles = [self._cycle("x", shared),
                  self._cycle("y", shared + make_dna(1000, 24))]
        assert len(filter_cycles_multiset(list(cycles), 21, 60, "Cycle")) == 1
        assert len(filter_cycles_multiset(list(cycles), 21, 90, "Cycle")) == 2

    def test_candidates_sharing_no_kmers_are_never_compared(self):
        """
        The inverted index only surfaces candidates sharing a k-mer, so two
        unrelated sequences stay separate even at a threshold of zero. Same
        property the legacy filter had; pinned so the optimisation is not
        mistaken for a scoring bug later.
        """
        cycles = [self._cycle("x", make_dna(4000, 25)),
                  self._cycle("y", make_dna(4000, 26))]
        assert len(filter_cycles_multiset(cycles, 21, 0, "Cycle")) == 2
