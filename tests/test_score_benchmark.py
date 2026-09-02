"""
Tests for scripts/score_picota_benchmark.py.

The precision definition is the delicate one and is tested hardest. An earlier
version asked only whether a candidate's own sequence was composite-transposon
derived, which a *fragment* of an element satisfies completely -- so partial
candidates covering 78% of an element were scored as true positives and
precision came out about four points too high. Being built from the right
sequence is not the same as being the right thing.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from score_picota_benchmark import (  # noqa: E402
    covered_positions,
    ct_id_from_subject,
    query_covered_fraction,
    read_ground_truth,
    score,
)


def blast_row(query, subject, sstart, send, pident=100.0, qlen=4000,
              qstart=1, qend=None):
    """One BLAST row in the tabular order the scorer expects."""
    qend = qend if qend is not None else qlen
    return [query, subject, str(pident), str(abs(send - sstart) + 1), str(qlen),
            "4000", str(qstart), str(qend), str(sstart), str(send)]


@pytest.fixture
def ground_truth(tmp_path):
    path = tmp_path / "ground_truth.tsv"
    path.write_text(
        "CT_ID\tCT_Length\tIS_Name\tCargo_Type\n"
        "CT001\t4000\tISfoo\tAMR\n"
        "CT002\t4000\tISfoo\tnovel\n")
    return read_ground_truth(str(path))


class TestSubjectParsing:
    def test_ct_id_taken_from_the_fasta_header(self):
        assert ct_id_from_subject("CT001_ISfoo_len4000") == "CT001"


class TestCoveredPositions:
    def test_overlapping_hsps_are_not_double_counted(self):
        """
        A cycle carrying two copies of one IS produces overlapping alignments;
        summing their lengths would push coverage past 100%.
        """
        rows = [blast_row("cyc1", "CT001_x", 1, 2000),
                blast_row("cyc1", "CT001_x", 1500, 3000)]
        covered = covered_positions(rows, 95.0)
        assert len(covered["CT001"]["cyc1"]) == 3000

    def test_low_identity_hits_are_ignored(self):
        rows = [blast_row("cyc1", "CT001_x", 1, 4000, pident=80.0)]
        assert covered_positions(rows, 95.0) == {}


class TestPrecisionRequiresACompleteElement:
    def test_candidate_covering_a_whole_element_counts(self, ground_truth):
        rows = [blast_row("cyc1", "CT001_x", 1, 4000)]
        result = score(ground_truth, rows, ["cyc1"], 95.0, 0.95)
        assert result["precision"] == (1, 1)

    def test_fragment_of_an_element_does_not_count(self, ground_truth):
        """
        The bug this test exists for. cyc1 is 2800 bp of pure CT001 sequence, so
        100% of the CANDIDATE is explained -- but it covers only 70% of the
        element, and reporting a fragment is not reporting the element.
        """
        rows = [blast_row("cyc1", "CT001_x", 1, 2800, qlen=2800)]
        result = score(ground_truth, rows, ["cyc1"], 95.0, 0.95)
        assert query_covered_fraction(rows, 95.0)["cyc1"] == pytest.approx(1.0)
        assert result["precision"] == (0, 1)

    def test_recall_and_precision_disagree_on_a_fragment(self, ground_truth):
        """A fragment is neither a recovery nor a true positive."""
        rows = [blast_row("cyc1", "CT001_x", 1, 2800, qlen=2800)]
        result = score(ground_truth, rows, ["cyc1"], 95.0, 0.95)
        assert result["ct_recall"] == (0, 2)
        assert result["precision"] == (0, 1)

    def test_two_candidates_for_one_element_both_count_as_precise(self, ground_truth):
        """
        Redundancy is not imprecision: both candidates really are the element.
        Recall still counts the element once.
        """
        rows = [blast_row("cyc1", "CT001_x", 1, 4000),
                blast_row("cyc2", "CT001_x", 1, 3900)]
        result = score(ground_truth, rows, ["cyc1", "cyc2"], 95.0, 0.95)
        assert result["precision"] == (2, 2)
        assert result["ct_recall"] == (1, 2)

    def test_candidate_matching_nothing_is_a_false_positive(self, ground_truth):
        result = score(ground_truth, [], ["cyc1"], 95.0, 0.95)
        assert result["precision"] == (0, 1)


class TestRecall:
    def test_every_element_recovered(self, ground_truth):
        rows = [blast_row("cyc1", "CT001_x", 1, 4000),
                blast_row("cyc2", "CT002_x", 1, 4000)]
        result = score(ground_truth, rows, ["cyc1", "cyc2"], 95.0, 0.95)
        assert result["ct_recall"] == (2, 2)

    def test_recall_split_by_cargo_type(self, ground_truth):
        """Novel cargo has to be readable separately; the score is homology-led."""
        rows = [blast_row("cyc1", "CT001_x", 1, 4000)]
        result = score(ground_truth, rows, ["cyc1"], 95.0, 0.95)
        assert result["by_cargo"]["AMR"] == [1, 1]
        assert result["by_cargo"]["novel"] == [0, 1]


class TestCopyDistinctness:
    def test_elements_sharing_an_is_reported_separately(self, ground_truth):
        rows = [blast_row("cyc1", "CT001_x", 1, 4000),
                blast_row("cyc2", "CT002_x", 1, 4000)]
        result = score(ground_truth, rows, ["cyc1", "cyc2"], 95.0, 0.95)
        assert result["copy_distinctness"] == (2, 2)

    def test_two_elements_collapsed_onto_one_candidate_count_once(self, ground_truth):
        """The failure the metric exists to catch."""
        rows = [blast_row("cyc1", "CT001_x", 1, 4000),
                blast_row("cyc1", "CT002_x", 1, 4000)]
        result = score(ground_truth, rows, ["cyc1"], 95.0, 0.95)
        assert result["copy_distinctness"] == (1, 2)
