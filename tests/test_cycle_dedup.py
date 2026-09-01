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
    estimated_ani,
    canonical_path_key,
    dedup_paths,
    filter_cycles_multiset,
    flip_orientation,
    length_ratio,
    multiset_jaccard,
    paths_are_duplicates,
    reverse_complement,
    reverse_traversal,
)
from src.cycle_finderv2 import Cycle  # noqa: E402
from synthetic_gfa import make_dna, mutate  # noqa: E402


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


class TestEstimatedAni:
    """
    D5: raw k-mer overlap is a cliff whose position depends on k. The Mash
    transform turns it into a quantity that tracks true identity, which is what
    makes a threshold meaningful rather than an artefact of the parameter.
    """

    @pytest.mark.parametrize("divergence", [0.001, 0.005, 0.01, 0.02, 0.05])
    def test_recovers_true_identity(self, divergence):
        seq = make_dna(4000, 20)
        mutated = mutate(seq, divergence, 7)
        jaccard = multiset_jaccard(canonical_kmers(seq, 21),
                                   canonical_kmers(mutated, 21))
        assert estimated_ani(jaccard, 21) == pytest.approx(1 - divergence, abs=0.005)

    @pytest.mark.parametrize("k", [21, 31, 80])
    def test_estimate_is_stable_across_k(self, k):
        """
        The raw Jaccard for this pair is 0.81 at k=21 and 0.48 at k=80. The
        identity estimate has to stay put where the Jaccard does not.
        """
        seq = make_dna(4000, 20)
        mutated = mutate(seq, 0.005, 7)
        jaccard = multiset_jaccard(canonical_kmers(seq, k), canonical_kmers(mutated, k))
        assert estimated_ani(jaccard, k) == pytest.approx(0.995, abs=0.002)

    def test_identical_sequences_estimate_one(self):
        assert estimated_ani(1.0, 21) == 1.0

    def test_no_overlap_estimates_zero(self):
        assert estimated_ani(0.0, 21) == 0.0

    def test_negative_jaccard_is_clamped(self):
        assert estimated_ani(-0.5, 21) == 0.0


class TestLengthRatio:
    def test_equal_lengths(self):
        assert length_ratio(4000, 4000) == 1.0

    def test_shorter_over_longer_regardless_of_order(self):
        assert length_ratio(2400, 3900) == length_ratio(3900, 2400)
        assert length_ratio(2400, 3900) == pytest.approx(0.615, abs=0.001)

    def test_zero_length(self):
        assert length_ratio(0, 0) == 0.0


class TestAniDeduplication:
    """The two criteria together, which is how strict mode actually runs."""

    def _cycles(self, seq_a, seq_b):
        return [Cycle("a", seq_a, len(seq_a), 2, []),
                Cycle("b", seq_b, len(seq_b), 2, [])]

    @pytest.mark.parametrize("k", [21, 80])
    def test_same_ct_at_half_percent_divergence_is_merged(self, k):
        """D5 fixed: assembly-level noise no longer looks like a distinct CT."""
        seq = make_dna(4000, 20)
        cycles = self._cycles(seq, mutate(seq, 0.005, 7))
        assert len(filter_cycles_multiset(cycles, k, 80, "Cycle", min_ani=0.99)) == 1

    @pytest.mark.parametrize("k", [21, 80])
    def test_complete_ct_not_merged_into_its_fragment(self, k):
        """
        D3 must stay fixed. At k=80 these estimate 99.65% identity, above the
        threshold -- only the length criterion keeps them apart.
        """
        is_element, cargo = make_dna(1500, 10), make_dna(900, 11)
        cycles = self._cycles(is_element + cargo, is_element + cargo + is_element)
        assert len(filter_cycles_multiset(cycles, k, 80, "Cycle", min_ani=0.99)) == 2

    @pytest.mark.parametrize("k", [21, 80])
    def test_contained_cycle_not_merged(self, k):
        shared = make_dna(1500, 10) + make_dna(900, 11)
        cycles = self._cycles(shared, shared + make_dna(2000, 12))
        assert len(filter_cycles_multiset(cycles, k, 80, "Cycle", min_ani=0.99)) == 2

    def test_genuine_duplicates_still_merged(self):
        seq = make_dna(4000, 30)
        assert len(filter_cycles_multiset(self._cycles(seq, seq), 21, 80,
                                          "Cycle", min_ani=0.99)) == 1

    def test_ani_path_is_order_independent(self):
        seq = make_dna(4000, 20)
        mutated = mutate(seq, 0.005, 7)
        forward = filter_cycles_multiset(self._cycles(seq, mutated), 21, 80,
                                         "Cycle", min_ani=0.99)
        backward = filter_cycles_multiset(self._cycles(mutated, seq), 21, 80,
                                          "Cycle", min_ani=0.99)
        assert len(forward) == len(backward) == 1

    def test_length_guard_can_be_relaxed(self):
        """min_length_ratio=0 reduces the rule to identity alone."""
        is_element, cargo = make_dna(1500, 10), make_dna(900, 11)
        cycles = self._cycles(is_element + cargo, is_element + cargo + is_element)
        assert len(filter_cycles_multiset(cycles, 80, 80, "Cycle", min_ani=0.99,
                                          min_length_ratio=0.0)) == 1
