"""
score3's terms, and the two defects that made them silent.

Every term here was at some point in the code doing nothing, or doing the
opposite of what its name said, without any test noticing:

  * the depth sidecar was never copied next to cycles.fasta, so the multi-copy
    term took its unknown-depth default for every candidate in every case, and
    that term is now gone: the ratio it used is bimodal, sitting below spurious
    cycles at two IS copies and above them at sixteen
  * the component tolerance was documented as scaling with a copy-number
    estimate that never arrived, so it was a constant in practice
  * the tool comparison unioned alignment fragments and credited matches that a
    single alignment never made

These tests pin the behaviour that replaced them.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "picota"))

from src.scoringv4ProtBlast import (  # noqa: E402
    _component_fit,
    calculate_total_score,
)

# (dist_type, max_z, mean, std) fixed across the score3 calls below so only the
# term under test varies.
DIST = (1, 3, 5000.0, 2000.0)


def score3(len_of_cycle=4800, ant=(95.0,), ises=(98.0,), xeno=(), comp=2):
    return calculate_total_score(3, *DIST, len_of_cycle, list(ant), list(ises),
                                 list(xeno), comp)


class TestComponentFit:
    def test_tolerance_is_constant(self):
        # Documented as scaling with multi-copy evidence, but that input never
        # reached the function and every published result used the constant.
        assert _component_fit(2) == pytest.approx(1.0)
        assert _component_fit(14) == pytest.approx(0.5)

    def test_more_components_score_lower(self):
        assert _component_fit(2) > _component_fit(11) > _component_fit(25)


class TestScore3:
    def test_no_insertion_sequence_scores_zero(self):
        # The IS is a necessary condition, not a bonus.
        assert score3(ises=()) == pytest.approx(0.0)

    def test_cargo_absent_from_databases_is_still_reportable(self):
        # The novel-cargo floor is what score2 lacks; a real element carrying
        # cargo in no database must still clear a threshold of 50.
        assert score3(ant=(), xeno=()) > 50.0

    def test_length_gates_rather_than_contributes(self):
        # A cycle far above the mean is rejected however good its homology,
        # which is what separates score3 from score2 on wild-type genomes.
        assert score3(len_of_cycle=36000) < 50.0

    def test_shorter_than_the_mean_pays_no_length_penalty(self):
        assert score3(len_of_cycle=3000) == pytest.approx(score3(len_of_cycle=4800))

    def test_weights_sum_to_one(self):
        # All three terms at their maximum, both gates open, must reach 100.
        perfect = calculate_total_score(3, *DIST, 4800, [100.0], [100.0], [], 2)
        assert perfect == pytest.approx(100.0)


class TestSidecarTravelsWithTheFasta:
    def test_run_scenarios_copies_the_depth_sidecar(self):
        # The sidecar sits next to cycles_<assembler>.fasta while scoring reads
        # cycles.fasta, so copying one without the other silences the term.
        source = open(os.path.join(os.path.dirname(__file__), "..", "scripts",
                                   "run_scenarios.py")).read()
        assert "cycles.depths.tsv" in source
        # copied, not merely mentioned
        assert source.count("shutil.copyfile") >= 2
