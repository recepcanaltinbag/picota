"""
Tests for total_score_type 3.

score2, the mode used for published analyses, puts 90 of its 100 points on a
single yes/no question -- does the cargo hit a database -- and compresses
everything structural into the remaining 10. Two consequences follow and both
are tested here: a candidate whose cargo is absent from every database can never
exceed 10 however good its structure, and a cycle of eighteen components scores
within four points of a clean two-component one, so the score labels rather than
ranks.

score3 keeps the same 0-100 range and inverts that balance.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "picota"))

from src.scoringv4ProtBlast import (  # noqa: E402
    _best,
    _component_fit,
    _length_fit,
    _multicopy_fit,
    calculate_total_score,
)

MEAN, STD, MAX_Z, DIST = 5850, 2586, 20, 1


def score(ant=(), is_=(100,), xe=(), length=5800, comp=2, depth_ratio=3.0,
          score_type=3):
    return calculate_total_score(score_type, DIST, MAX_Z, MEAN, STD, length,
                                 list(ant), list(is_), list(xe), comp,
                                 depth_ratio)


class TestComponentTerms:
    def test_length_fit_is_one_at_the_distribution_mean(self):
        assert _length_fit(MEAN, MEAN, STD, MAX_Z, DIST) == pytest.approx(1.0)

    def test_length_fit_falls_for_over_long_candidates(self):
        assert _length_fit(MEAN + 4 * STD, MEAN, STD, MAX_Z, DIST) < 0.85

    def test_short_candidates_are_not_penalised_under_dist_type_1(self):
        """Preserves the existing convention rather than changing it here."""
        assert _length_fit(1000, MEAN, STD, MAX_Z, 1) == pytest.approx(1.0)

    def test_component_fit_peaks_at_two(self):
        assert _component_fit(2) == pytest.approx(1.0)

    def test_component_fit_spreads_across_the_real_range(self):
        """
        The failure in score2: sqrt(|comp-2|)/20 moved the total under four
        points from two components to eighteen, so structure could not rank.
        """
        assert _component_fit(2) - _component_fit(18) > 0.6

    def test_component_fit_is_monotone(self):
        values = [_component_fit(c) for c in range(2, 25)]
        assert all(a >= b for a, b in zip(values, values[1:]))


class TestMulticopyTerm:
    def test_single_copy_scores_zero(self):
        """A composite transposon needs two copies of its flanking IS."""
        assert _multicopy_fit(1.0) == 0.0

    def test_two_and_a_half_copies_saturate(self):
        assert _multicopy_fit(2.5) == pytest.approx(1.0)

    def test_unknown_depth_is_neutral_not_zero(self):
        """An assembly reporting no coverage must not be penalised for it."""
        assert _multicopy_fit(None) == 0.5

    def test_below_one_is_clamped(self):
        assert _multicopy_fit(0.4) == 0.0


class TestBestHitQuality:
    def test_uses_the_best_hit_not_the_count(self):
        """
        score0 sums hit scores, so a cargo of three genes outscores one of a
        single gene for carrying more genes rather than for being more likely a
        composite transposon.
        """
        assert _best([90, 90, 90]) == _best([90])

    def test_no_hits_is_zero(self):
        assert _best([]) == 0.0

    def test_capped_at_one(self):
        assert _best([250]) == 1.0


class TestScore3:
    def test_ideal_candidate_reaches_one_hundred(self):
        assert score(ant=(100,)) == pytest.approx(100.0)

    def test_no_is_scores_zero(self):
        """An IS is a necessary condition, so it gates rather than contributes."""
        assert score(ant=(100,), is_=()) == 0.0

    def test_weak_is_lowers_the_score_without_halving_it(self):
        """
        The gate is presence-with-quality, not quality alone. A hit at 40%
        identity-coverage is still unambiguously an insertion sequence; scaling
        the whole score by 0.4 for it compounded with the component penalty and
        put genuine shared-IS candidates 0.1 points under the threshold.
        """
        weak, strong = score(ant=(100,), is_=(40,)), score(ant=(100,), is_=(100,))
        assert weak < strong
        assert weak > strong * 0.6

    def test_novel_cargo_stays_reportable(self):
        """
        The point of the exercise. Under score2 this candidate scores 10 and is
        unreportable at any sensible threshold, however good its structure.
        """
        novel = score(ant=())
        assert novel > 50
        assert novel < score(ant=(100,))

    def test_single_copy_is_penalised(self):
        """Homology cannot rescue a structure that is not multi-copy."""
        assert score(ant=(100,), depth_ratio=1.0) < score(ant=(100,), depth_ratio=3.0)

    def test_components_move_the_score_for_a_single_copy_element(self):
        """
        Where many components are genuinely suspicious -- a two-copy element
        should assemble into a clean bubble -- the penalty still bites.
        """
        spread = (score(ant=(100,), comp=2, depth_ratio=1.0)
                  - score(ant=(100,), comp=18, depth_ratio=1.0))
        assert spread > 10, "structure must rank, not merely label"

    def test_component_tolerance_widens_with_multi_copy_evidence(self):
        """
        The two signals are not independent. A cycle threading eighteen nodes is
        evidence against a single-copy element and expected of one whose IS sits
        in dozens of copies; charging both the same penalises the very structure
        being detected.
        """
        single = score(ant=(100,), comp=18, depth_ratio=1.0)
        multi = score(ant=(100,), comp=18, depth_ratio=54.0)
        assert multi > single + 15

    def test_multi_copy_evidence_outranks_an_unexplained_complex_cycle(self):
        """
        The discrimination the depth term exists for, and the reason it is not
        optional: without it this ordering inverts, and a single-copy element
        with a complex cycle outscores a genuine shared-IS composite transposon.
        """
        real_shared = score(ant=(100,), is_=(59,), length=3560, comp=18,
                            depth_ratio=54.0)
        unexplained = score(ant=(100,), is_=(100,), comp=18, depth_ratio=1.0)
        assert real_shared > unexplained

    def test_single_copy_clean_cycle_is_not_a_perfect_score(self):
        """A composite transposon needs two copies; one copy is not one."""
        assert score(ant=(100,), comp=2, depth_ratio=1.0) < 80

    def test_score2_barely_moves_over_the_same_range(self):
        """The comparison that motivates score3, pinned so it cannot be forgotten."""
        high = score(ant=(100,), comp=2, score_type=2)
        low = score(ant=(100,), comp=18, score_type=2)
        assert high - low < 5

    def test_stays_within_zero_and_one_hundred(self):
        for comp in (1, 2, 10, 25):
            for ratio in (None, 0.5, 1.0, 5.0):
                for ant in ((), (50,), (100,)):
                    value = score(ant=ant, comp=comp, depth_ratio=ratio)
                    assert 0.0 <= value <= 100.0

    def test_unknown_depth_does_not_break_ranking(self):
        assert score(ant=(100,), depth_ratio=None) > score(ant=(100,), depth_ratio=1.0)

    def test_invalid_score_type_still_raises(self):
        with pytest.raises(Exception):
            calculate_total_score(4, DIST, MAX_Z, MEAN, STD, 5800, [], [100], [], 2)
