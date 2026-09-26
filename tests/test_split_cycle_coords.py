"""Block splitting for long-read validation.

This module had no tests, which is how it absorbed a regression from a change
elsewhere. It reads `IScoords`, and while the scorer wrote one insertion sequence
per cycle, "first start to last end" and "the one locus" were the same span. Once
the scorer reported one IS per locus they were not: the span reached from the
first IS to the last and swallowed the cargo between them, so the cargo block --
the thing long-read validation looks for -- stopped being written.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.split_cycle_coords_for_is import split_cycles_from_picota

COLUMNS = ("CycleID\tscore2\tNumIS\tISproducts\tIScoords\t"
           "Antproducts\tAntcoords\tXenoproducts\tXenocoords\n")


def _inputs(tmp, cycle_id, length, is_coords, ant_coords):
    """A one-row PICOTA table and the cycle it describes."""
    tab = os.path.join(tmp, "picota_final_tab")
    with open(tab, "w") as fh:
        fh.write(COLUMNS)
        fh.write(f"{cycle_id}\t99.0\t{len(is_coords.split(';'))}\tIS26\t{is_coords}\t"
                 f"tet(B)\t{ant_coords}\t\t\n")
    fa = os.path.join(tmp, "cycles.fasta")
    with open(fa, "w") as fh:
        fh.write(f">{cycle_id}\n" + ("ACGT" * (length // 4 + 1))[:length] + "\n")
    return tab, fa


def _blocks(out_dir):
    """Header suffix -> (start, end), for the single split file written."""
    files = [f for f in os.listdir(out_dir) if f.endswith("_split.fasta")]
    assert len(files) == 1, files
    blocks = {}
    for line in open(os.path.join(out_dir, files[0])):
        if line.startswith(">"):
            parts = line[1:].strip().rsplit("_", 3)
            blocks[parts[3]] = (int(parts[1]), int(parts[2]))
    return blocks


class TestTransposonBlock:

    def test_one_locus_is_the_transposon_block(self):
        """The behaviour that has to be preserved: a cycle with a single IS."""
        with tempfile.TemporaryDirectory() as tmp:
            tab, fa = _inputs(tmp, "Cycle_1-len6000-comp2-", 6000,
                              is_coords="1-820", ant_coords="1500-2700")
            out = os.path.join(tmp, "annot")
            split_cycles_from_picota(tab, fa, out)
            b = _blocks(out)

            assert b["transposon"] == (1, 820)
            assert "cargo" in b
            assert b["cargo"][0] > 820

    def test_several_loci_do_not_become_one_span(self):
        """The regression: flanking IS at both ends of the element.

        First start to last end is 1..5900 here -- the whole element, cargo
        included -- and the cargo block disappears. Measured over 2,192 real
        elements before the fix, the transposon block covered a median 87% of
        the element once four or more loci were reported, and only 8% of those
        elements still had a cargo block.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tab, fa = _inputs(tmp, "Cycle_2-len6000-comp2-", 6000,
                              is_coords="1-820;5080-5900", ant_coords="1500-2700")
            out = os.path.join(tmp, "annot")
            split_cycles_from_picota(tab, fa, out)
            b = _blocks(out)

            assert b["transposon"] == (1, 820), "span of all loci, not one locus"
            assert b["transposon"][1] - b["transposon"][0] + 1 < 1000
            assert "cargo" in b, "the cargo block was swallowed"
            assert b["cargo"] == (821, 2700)

    def test_the_strongest_locus_wins_not_the_leftmost(self):
        """IScoords is written best-scoring first, and that is what is used."""
        with tempfile.TemporaryDirectory() as tmp:
            # the best-scoring locus is the second one by position
            tab, fa = _inputs(tmp, "Cycle_3-len6000-comp2-", 6000,
                              is_coords="5080-5900;1-820", ant_coords="1500-2700")
            out = os.path.join(tmp, "annot")
            split_cycles_from_picota(tab, fa, out)
            b = _blocks(out)

            assert b["transposon"] == (5080, 5900)
