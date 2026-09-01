"""
Tests for node depth parsing and the report-only depth summary (roadmap phase 1).

Depth is the signal that distinguishes a multi-copy IS node from a single-copy
region. Nothing in detection, scoring or filtering is allowed to depend on it
yet -- these tests pin the parsing and the sidecar report only.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "picota"))
sys.path.insert(0, os.path.dirname(__file__))

from src.cycle_finderv2 import (  # noqa: E402
    Cycle,
    GraphWork,
    cycle_analysis,
    parse_segment_depth,
    write_depth_report,
)
from synthetic_gfa import shared_repeat_n_cargos  # noqa: E402


class TestParseSegmentDepth:
    """Every depth encoding PICOTA can encounter in the wild."""

    def test_spades_name_encoding(self):
        assert parse_segment_depth("NODE_1_length_523_cov_9.5794_ID_1", "", 523) == pytest.approx(9.5794)

    def test_megahit_multi_in_tags(self):
        assert parse_segment_depth("k99_12", "flag=1 multi=12.3456 len=400", 400) == pytest.approx(12.3456)

    def test_megahit_multi_in_name(self):
        assert parse_segment_depth("k99_12_multi=7.5", "", 400) == pytest.approx(7.5)

    def test_gfa_dp_tag(self):
        assert parse_segment_depth("1", "LN:i:61\tdp:f:35.2", 61) == pytest.approx(35.2)

    def test_gfa_dp_tag_uppercase(self):
        assert parse_segment_depth("1", "LN:i:61\tDP:f:35.2", 61) == pytest.approx(35.2)

    def test_kmer_count_falls_back_to_kc_over_ln(self):
        assert parse_segment_depth("1", "LN:i:61\tKC:i:2440", 61) == pytest.approx(40.0)

    def test_kc_without_ln_uses_sequence_length(self):
        assert parse_segment_depth("1", "KC:i:2440", 61) == pytest.approx(40.0)

    def test_unknown_depth_is_none_not_zero(self):
        """None means 'assembler said nothing'. Zero would be a real measurement."""
        assert parse_segment_depth("1", "", 61) is None

    def test_spades_name_wins_over_tags(self):
        assert parse_segment_depth("NODE_1_length_523_cov_9.5_ID_1", "dp:f:99.0", 523) == pytest.approx(9.5)


class TestParseGfaDepth:
    def test_depth_attached_to_both_strands(self, tmp_path):
        gfa = shared_repeat_n_cargos(str(tmp_path / "g.gfa"), 3,
                                     repeat_depth=90.0, cargo_depth=30.0)
        node_dict, _ = GraphWork().parse_gfa(gfa)
        for name in ("repeat+", "repeat-"):
            assert node_dict[name]["Depth"] == pytest.approx(90.0, rel=0.05)
        assert node_dict["cargo0-"]["Depth"] == pytest.approx(30.0, rel=0.05)

    def test_untagged_gfa_yields_none(self, tmp_path):
        gfa = tmp_path / "plain.gfa"
        gfa.write_text("S\t1\tATCGATCG\nS\t2\tTTTTGGGG\nL\t1\t+\t2\t+\t0M\n")
        node_dict, _ = GraphWork().parse_gfa(str(gfa))
        assert node_dict["1+"]["Depth"] is None


class TestDepthRatio:
    def test_ratio_approximates_repeat_copy_number(self):
        cycle = Cycle("c", "ACGT", 4, 2, ["repeat+", "cargo+"], [90.0, 30.0])
        assert cycle.depth_ratio == pytest.approx(3.0)

    def test_single_copy_bubble_has_ratio_near_one(self):
        cycle = Cycle("c", "ACGT", 4, 2, ["a+", "b+"], [31.0, 30.0])
        assert cycle.depth_ratio == pytest.approx(1.033, rel=0.01)

    def test_none_when_depth_unknown(self):
        assert Cycle("c", "ACGT", 4, 2, ["a+", "b+"], [None, None]).depth_ratio is None

    def test_none_when_only_one_node_has_depth(self):
        assert Cycle("c", "ACGT", 4, 1, ["a+"], [30.0]).depth_ratio is None

    def test_zero_depth_is_ignored_not_divided_by(self):
        cycle = Cycle("c", "ACGT", 4, 3, ["a+", "b+", "c+"], [90.0, 0.0, 30.0])
        assert cycle.depth_ratio == pytest.approx(3.0)

    def test_default_node_depths_is_empty(self):
        """Cycle stays constructible with the pre-phase-1 signature."""
        assert Cycle("c", "ACGT", 4, 2, ["a+", "b+"]).depth_ratio is None


class TestDepthReport:
    def test_report_written_next_to_cycle_fasta(self, tmp_path):
        cycles = [Cycle("Cycle_1", "ACGT", 4, 2, ["repeat+", "cargo0+"], [90.0, 30.0])]
        out = write_depth_report(str(tmp_path / "sample_cycles.fasta"), cycles)
        assert os.path.basename(out) == "sample_cycles.depths.tsv"

    def test_report_contents(self, tmp_path):
        cycles = [
            Cycle("Cycle_1", "ACGT", 4000, 2, ["repeat+", "cargo0+"], [90.0, 30.0]),
            Cycle("Cycle_2", "ACGT", 3000, 2, ["a+", "b+"], [None, None]),
        ]
        out = write_depth_report(str(tmp_path / "s.fasta"), cycles)
        rows = [ln.rstrip("\n").split("\t") for ln in open(out)]
        assert rows[0][0] == "CycleID" and "DepthRatio" in rows[0]
        assert rows[1][0] == "Cycle_1-len4000-comp2-"
        assert rows[1][rows[0].index("DepthRatio")] == "3.000"
        assert rows[2][rows[0].index("DepthRatio")] == "NA"
        assert rows[2][rows[0].index("NodeDepths")] == "NA;NA"

    def test_cycle_analysis_emits_the_report(self, tmp_path):
        """The sidecar is produced by a real run, not only by direct calls."""
        gfa = shared_repeat_n_cargos(str(tmp_path / "g.gfa"), 3,
                                     repeat_len=2500, cargo_len=900,
                                     repeat_depth=90.0, cargo_depth=30.0)
        out_fasta = str(tmp_path / "out_cycles.fasta")
        cycle_analysis(gfa, out_fasta, True, 25, 2000, 40000, "Cycle", 1, 25, 80, 80)

        report = str(tmp_path / "out_cycles.depths.tsv")
        assert os.path.exists(report)
        rows = [ln.rstrip("\n").split("\t") for ln in open(report)]
        assert len(rows) > 1, "no cycles reported"
        ratios = [r[rows[0].index("DepthRatio")] for r in rows[1:]]
        assert all(r != "NA" for r in ratios)
        # repeat node at 90x against 30x cargo -> ratio 3
        assert all(float(r) == pytest.approx(3.0, rel=0.05) for r in ratios)


class TestDepthIsReportOnly:
    def test_cycle_count_unchanged_by_depth_tags(self, tmp_path):
        """Phase 1 must not alter which cycles are reported (see docs/ROADMAP.md)."""
        from test_known_defects import run_cycle_pipeline

        with_depth = shared_repeat_n_cargos(str(tmp_path / "d.gfa"), 3,
                                            repeat_depth=90.0, cargo_depth=30.0)
        # Same graph, depth tags stripped
        plain = tmp_path / "p.gfa"
        stripped = []
        for line in open(with_depth):
            fields = line.rstrip("\n").split("\t")
            kept = [f for f in fields if not f.startswith(("KC:i:", "dp:f:"))]
            stripped.append("\t".join(kept) + "\n")
        plain.write_text("".join(stripped))

        assert len(run_cycle_pipeline(with_depth)) == len(run_cycle_pipeline(str(plain)))


class TestNodeSpans:
    """
    Where each node sits inside the finished cycle. Without this the depth of a
    particular region -- the IS, say -- cannot be looked up at all.
    """

    def test_spans_cover_the_whole_cycle_without_gaps(self, tmp_path):
        gfa = shared_repeat_n_cargos(str(tmp_path / "g.gfa"), 2)
        from src.cycle_finderv2 import GraphWork, cycle_info_optimized
        graph_work = GraphWork()
        graph_work.find_all_path = True
        node_dict, edge_dict = graph_work.parse_gfa(gfa)
        graph_work.dfs_iterative(graph_work.generate_genome_graph(node_dict, edge_dict))
        cycle = cycle_info_optimized(graph_work.allPaths[0], node_dict, edge_dict, [])

        assert len(cycle.node_spans) == len(cycle.path)
        assert cycle.node_spans[0][0] == 1
        assert cycle.node_spans[-1][1] == cycle.length
        for (_, end), (next_start, _) in zip(cycle.node_spans, cycle.node_spans[1:]):
            assert next_start == end + 1

    def test_spans_appear_in_the_sidecar(self, tmp_path):
        cycles = [Cycle("Cycle_1", "ACGT", 4, 2, ["a+", "b+"], [90.0, 30.0],
                        [(1, 2), (3, 4)])]
        out = write_depth_report(str(tmp_path / "s.fasta"), cycles)
        rows = [ln.rstrip("\n").split("\t") for ln in open(out)]
        assert "NodeSpans" in rows[0]
        assert rows[1][rows[0].index("NodeSpans")] == "1-2;3-4"


class TestDepthAt:
    def _cycle(self):
        return Cycle("c", "N" * 4000, 4000, 3, ["is+", "cargo+", "tail+"],
                     [120.0, 30.0, 28.0], [(1, 1500), (1501, 3500), (3501, 4000)])

    def test_region_inside_one_node(self):
        assert self._cycle().depth_at(100, 900) == 120.0

    def test_region_inside_a_later_node(self):
        assert self._cycle().depth_at(2000, 3000) == 30.0

    def test_region_straddling_two_nodes_takes_the_median(self):
        assert self._cycle().depth_at(1400, 1600) == pytest.approx(75.0)

    def test_region_outside_every_node(self):
        assert self._cycle().depth_at(9000, 9500) is None

    def test_nodes_without_depth_are_skipped(self):
        cycle = Cycle("c", "N" * 100, 100, 2, ["a+", "b+"], [None, 40.0],
                      [(1, 50), (51, 100)])
        assert cycle.depth_at(1, 50) is None
        assert cycle.depth_at(60, 90) == 40.0


class TestJunctionDepthRatio:
    """
    The comparison the reviewer asked for: IS depth against cargo depth.

    The arithmetic is exact; what it means is not. Against known ground truth the
    ratio understates copy number severely -- 16 genomic copies read as 1.8 to
    3.8 -- so these tests pin the calculation, not an interpretation of it.
    """

    def _cycle(self, is_depth, cargo_depth):
        return Cycle("c", "N" * 4000, 4000, 2, ["is+", "cargo+"],
                     [is_depth, cargo_depth], [(1, 1500), (1501, 4000)])

    def test_ratio_is_is_depth_over_cargo_depth(self):
        cycle = self._cycle(120.0, 30.0)
        assert cycle.junction_depth_ratio([(200, 1400)], [(1600, 3900)]) == pytest.approx(4.0)

    def test_equal_depths_give_a_ratio_near_one(self):
        """An IS at the same depth as its cargo is not a multi-copy element."""
        cycle = self._cycle(31.0, 30.0)
        ratio = cycle.junction_depth_ratio([(200, 1400)], [(1600, 3900)])
        assert ratio == pytest.approx(1.033, rel=0.01)

    def test_none_when_depth_unknown(self):
        cycle = self._cycle(None, 30.0)
        assert cycle.junction_depth_ratio([(200, 1400)], [(1600, 3900)]) is None

    def test_none_when_cargo_has_no_depth(self):
        cycle = self._cycle(120.0, None)
        assert cycle.junction_depth_ratio([(200, 1400)], [(1600, 3900)]) is None

    def test_several_is_regions_take_the_deepest(self):
        """Both flanking copies map to one node; the deepest is the shared one."""
        cycle = Cycle("c", "N" * 4000, 4000, 3, ["is+", "cargo+", "is2+"],
                      [120.0, 30.0, 35.0], [(1, 1200), (1201, 3000), (3001, 4000)])
        assert cycle.junction_depth_ratio([(100, 1100), (3100, 3900)],
                                          [(1400, 2900)]) == pytest.approx(4.0)

    def test_empty_regions(self):
        assert self._cycle(120.0, 30.0).junction_depth_ratio([], [(1600, 3900)]) is None
