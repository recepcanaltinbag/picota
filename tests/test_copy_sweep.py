"""
Tests for the copy-number axis.

Two things are pinned here. The first is that --is-copies-per-element does what
it claims, in the genome and not merely in the ground-truth column: an earlier
copy-count bug reported a number the sequence did not contain, and a copy count
that lies is worse than no copy count, because every depth-based conclusion is
compared against it.

The second is the benchmark harness threshold. run_scenarios.py hardcoded
min_size_of_cycle to 2000 while config.yaml moved to 1000, so the compact
scenario kept reporting the old threshold's result -- 8 of 12 elements -- long
after the change that was supposed to fix it. Nothing failed; the number simply
stayed wrong.
"""

import os
import subprocess
import sys

import pytest

SCRIPT_DIR = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, os.path.abspath(SCRIPT_DIR))

SIMULATE = os.path.join(SCRIPT_DIR, "simulate_ct_genome.py")


def read_ground_truth_copies(path):
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        column = header.index("IS_Genome_Copies")
        return [int(line.rstrip("\n").split("\t")[column]) for line in handle if line.strip()]


def simulate(tmp_path, **kwargs):
    out = str(tmp_path / "sim")
    command = [sys.executable, SIMULATE, "--out-dir", out,
               "--backbone-length", "600000", "--n-cts", "3",
               "--shared-is", "0", "--spacing", "5000", "--seed", "11"]
    for key, value in kwargs.items():
        command += ["--" + key.replace("_", "-"), str(value)]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    return out


class TestPerElementCopies:
    @pytest.mark.parametrize("extra,expected", [(0, 2), (2, 4), (6, 8)])
    def test_ground_truth_counts_two_plus_n(self, tmp_path, extra, expected):
        out = simulate(tmp_path, is_copies_per_element=extra)
        counts = read_ground_truth_copies(os.path.join(out, "ground_truth.tsv"))
        assert counts == [expected] * 3

    def test_the_genome_really_contains_them(self, tmp_path):
        """
        The column is only worth having if the sequence agrees with it: an
        earlier copy-count bug reported a number the sequence did not contain.

        Run at zero divergence so the assertion is exact. Every copy -- the two
        flanking ones and the six free-standing -- is otherwise mutated
        independently, so at the usual 0.5% no two are byte-identical and a
        count of exact matches would measure the mutation rate rather than the
        placement. Placement is what this test is for; divergence is covered by
        test_ground_truth_counts_two_plus_n through the recorded column.
        """
        out = simulate(tmp_path, is_copies_per_element=6, is_divergence=0)
        genome = "".join(
            line.strip() for line in open(os.path.join(out, "genome.fasta"))
            if not line.startswith(">"))

        with open(os.path.join(out, "ground_truth.tsv")) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            first = handle.readline().rstrip("\n").split("\t")
        is_length = int(first[header.index("IS_Length")])

        element = "".join(
            open(os.path.join(out, "ground_truth_cts.fasta")).read()
            .split("\n")[1:]).strip()
        flank = element[:is_length]

        assert genome.count(flank) == 8, \
            "2 flanking + 6 free-standing copies were requested; the genome " \
            "holds %d" % genome.count(flank)

    def test_shared_is_axis_is_untouched(self, tmp_path):
        """
        The new flag must not change what --shared-is did, or every published
        scenario number silently moves.
        """
        out = simulate(tmp_path, shared_is=3, is_copies_outside=4)
        counts = read_ground_truth_copies(os.path.join(out, "ground_truth.tsv"))
        assert counts == [10, 10, 10], "3 elements x 2 flanks + 4 free-standing"

    def test_default_is_off(self, tmp_path):
        out = simulate(tmp_path)
        assert read_ground_truth_copies(
            os.path.join(out, "ground_truth.tsv")) == [2, 2, 2]


class TestHarnessThreshold:
    def test_run_scenarios_does_not_hardcode_the_cycle_size(self):
        """
        Regression: the literal 2000 sat in the cycle_analysis call, so the
        compact scenario was still measured at the retired threshold after
        config.yaml moved to 1000.
        """
        source = open(os.path.join(SCRIPT_DIR, "run_scenarios.py")).read()
        call = source[source.index("cycle_analysis(gfa"):]
        call = call[:call.index(")")]
        assert "2000" not in call, "cycle size must come from --min-cycle-size"
        assert "args.min_cycle_size" in call

    def test_default_matches_the_shipped_config(self):
        import yaml

        sys.path.insert(0, os.path.abspath(SCRIPT_DIR))
        import run_scenarios

        config = os.path.join(os.path.dirname(__file__), "..", "picota",
                              "config.yaml")
        with open(config) as handle:
            shipped = yaml.safe_load(handle)

        parser_default = run_scenarios.build_parser().get_default("min_cycle_size")
        found = _find_min_size(shipped)
        assert found == parser_default, (
            "benchmark default %d and config.yaml %s disagree; the compact "
            "scenario measures whichever the harness happens to use"
            % (parser_default, found))


def _find_min_size(node):
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "min_size_of_cycle":
                return int(value)
            found = _find_min_size(value)
            if found is not None:
                return found
    return None
