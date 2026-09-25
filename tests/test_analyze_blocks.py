"""Gap acceptance in check_adjacent.

This module had no tests, which is how an UnboundLocalError survived in it: the
boost floor was assigned only inside the `is_near_mode` branch, immediately
before that branch's `continue`, and read only on the path taken when
`is_near_mode` is False. It therefore raised whenever no earlier pair in the same
read had been near the mode -- and when an earlier pair had been, the value
leaked across loop iterations, so a gap's verdict depended on an unrelated
earlier gap.

Both properties are pinned below: the path is reachable, and a pair's verdict
does not depend on what preceded it.
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.analyze_blocksv3 import check_adjacent
from src.config_loader import ToleranceConfig


class TestGapBoostFloor:
    """The branch that reads the boost floor, reached with nothing before it."""

    def test_a_lone_not_near_mode_gap_is_decided_and_does_not_raise(self):
        """The first pair of a read, not near the mode, reaching the boost check.

        Before the fix this raised UnboundLocalError. The gap is 500 against a
        mode of 1000: too far from the mode to be accepted there, too small to be
        "large", and no reverse statistics exist -- so the final boost check is
        the only thing left to decide it.
        """
        cfg = ToleranceConfig()
        coords = [(0, 1000), (1500, 2500)]          # gap = 500
        types = ["transposon", "cargo"]
        gap_summary = {
            ("transposon", "cargo"): {"mode_gap": 1000, "count": 10, "large_count": 1},
        }

        # floor 1000 vs mode*1.5 = 1500 -> 1500, and 500 <= 1500
        assert check_adjacent(coords, types, gap_summary, cfg) is True

    def test_a_gap_above_the_floor_is_rejected(self):
        """The same path, on the other side of the threshold."""
        cfg = ToleranceConfig()
        coords = [(0, 1000), (3000, 4000)]          # gap = 2000
        types = ["transposon", "cargo"]
        gap_summary = {
            ("transposon", "cargo"): {"mode_gap": 1000, "count": 0, "large_count": 0},
        }

        # floor 1000 vs mode*1.5 = 1500 -> 1500, and 2000 > 1500
        assert check_adjacent(coords, types, gap_summary, cfg) is False

    def test_the_verdict_does_not_depend_on_an_earlier_pair(self):
        """A pair decides the same way alone as it does after a near-mode pair.

        This is the leak, and it is the half of the defect that produced wrong
        answers rather than a traceback. The second pair here is identical in
        both reads; the first read has nothing before it, the second has a pair
        sitting exactly on a mode of 4000. Under the old code the first raised
        and the second returned True, because the earlier pair had left a floor
        of 3000 behind -- three times the configured one.
        """
        cfg = ToleranceConfig()
        gap_summary = {
            ("cargo", "transposon"): {"mode_gap": 4000, "count": 10, "large_count": 5},
            # count 0 keeps the large-gap branch from deciding this pair, so it
            # reaches the boost check where the floor is read.
            ("transposon", "unknown"): {"mode_gap": 600, "count": 0, "large_count": 0},
        }

        alone = check_adjacent(
            [(0, 1000), (2500, 3000)],                      # gap = 1500
            ["transposon", "unknown"], gap_summary, cfg)

        preceded = check_adjacent(
            [(0, 1000), (5000, 6000), (7500, 8000)],        # gaps = 4000, 1500
            ["cargo", "transposon", "unknown"], gap_summary, cfg)

        # floor 1000 vs mode*1.5 = 900 -> 1000, and 1500 > 1000
        assert alone is False
        assert preceded == alone

    def test_a_gap_on_the_mode_is_accepted_without_reaching_the_floor(self):
        """The near-mode branch still short-circuits, as it did before."""
        cfg = ToleranceConfig()
        coords = [(0, 1000), (4000, 5000)]          # gap = 3000, mode 3000
        types = ["transposon", "cargo"]
        gap_summary = {
            ("transposon", "cargo"): {"mode_gap": 3000, "count": 10, "large_count": 1},
        }

        assert check_adjacent(coords, types, gap_summary, cfg) is True
