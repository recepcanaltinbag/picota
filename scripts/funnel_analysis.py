#!/usr/bin/env python3
"""
funnel_analysis.py
------------------
Where elements are lost, stage by stage, and the standard metrics.

The scenario table reports TP/FN/FP but not *where* a missed element was
missed. An element can be lost in three distinct places, and they mean very
different things about the method:

  assembly + traversal   the element never became a candidate cycle. Either the
                         assembler did not preserve the structure or the cycle
                         search did not close a path through it. Nothing
                         downstream can recover it.
  size filter            a cycle covered the element but fell below
                         min_size_of_cycle. A parameter choice, not a failure.
  scoring                a cycle covered the element and was long enough, but
                         scored below the threshold. This is the only stage the
                         scoring model is responsible for.

Reporting them together as "FN" hides which component needs work. Separating
them showed that scoring is responsible for a small minority of losses, and
that is the claim this script exists to support or refute.

Also reports sensitivity and PPV, which the scenario table leaves for the
reader to compute:

  sensitivity = TP / (TP + FN)   of the implanted elements, how many are reported
  PPV         = TP / (TP + FP)   of the reported candidates, how many are real

Usage:
  python scripts/funnel_analysis.py --scenarios bench/scenarios_v2
"""

import argparse
import collections
import glob
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from score_scenarios import SCENARIO_ORDER, evaluate  # noqa: E402

TOOLS = {"path_of_prodigal": "prodigal", "path_of_blastn": "blastn",
         "path_of_makeblastdb": "makeblastdb", "path_of_blastx": "blastx",
         "path_of_blastp": "blastp"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--scenarios", required=True)
    parser.add_argument("--threshold", type=int, default=50)
    parser.add_argument("--score-type", type=int, default=3)
    args = parser.parse_args(argv)

    rows = []
    grand = collections.Counter()
    for name in SCENARIO_ORDER:
        cases = sorted(d for d in glob.glob(os.path.join(args.scenarios, name + "_s*"))
                       if os.path.exists(os.path.join(d, "cycles.fasta")))
        if not cases:
            continue
        pooled = collections.Counter()
        for case in cases:
            r = evaluate(case, args.threshold, TOOLS, args.score_type)
            for key in ("cts", "detected", "detected_tp", "scored", "tp", "fn", "fp"):
                pooled[key] += r[key]
        rows.append((name, len(cases), pooled))
        for key, value in pooled.items():
            grand[key] += value

    print("Where implanted elements are lost")
    print()
    header = ("%-15s %5s %9s %11s %9s %11s %9s"
              % ("scenario", "impl.", "in graph", "lost there", "reported",
                 "lost scoring", "recovery"))
    print(header)
    print("-" * len(header))
    for name, _, p in rows:
        lost_detect = p["cts"] - p["detected_tp"]
        lost_score = p["detected_tp"] - p["tp"]
        print("%-15s %5d %9d %11d %9d %11d %8.0f%%"
              % (name, p["cts"], p["detected_tp"], lost_detect, p["tp"],
                 lost_score, 100 * p["tp"] / p["cts"] if p["cts"] else 0))
    lost_detect = grand["cts"] - grand["detected_tp"]
    lost_score = grand["detected_tp"] - grand["tp"]
    print("-" * len(header))
    print("%-15s %5d %9d %11d %9d %11d %8.0f%%"
          % ("all", grand["cts"], grand["detected_tp"], lost_detect, grand["tp"],
             lost_score, 100 * grand["tp"] / grand["cts"]))
    print()
    print("in graph     = elements a candidate cycle covers at >=95%, before scoring")
    print("lost there   = never became a candidate: assembly or traversal")
    print("lost scoring = covered by a candidate, but scored below the threshold")
    print()

    print("Sensitivity and PPV")
    print()
    header = ("%-15s %5s %5s %5s %13s %13s"
              % ("scenario", "TP", "FN", "FP", "sensitivity", "PPV"))
    print(header)
    print("-" * len(header))
    for name, _, p in rows:
        sens = p["tp"] / (p["tp"] + p["fn"]) if p["tp"] + p["fn"] else 0
        ppv = p["tp"] / (p["tp"] + p["fp"]) if p["tp"] + p["fp"] else 0
        print("%-15s %5d %5d %5d %12.1f%% %12.1f%%"
              % (name, p["tp"], p["fn"], p["fp"], 100 * sens, 100 * ppv))
    print("-" * len(header))
    sens = grand["tp"] / (grand["tp"] + grand["fn"])
    ppv = grand["tp"] / (grand["tp"] + grand["fp"])
    print("%-15s %5d %5d %5d %12.1f%% %12.1f%%"
          % ("all", grand["tp"], grand["fn"], grand["fp"], 100 * sens, 100 * ppv))

    # cargo_is_same is a known structural failure and dominates the pooled
    # figures; reporting the rest separately is not hiding it, it is saying
    # what the method does where it is applicable.
    rest = collections.Counter()
    for name, _, p in rows:
        if name == "cargo_is_same":
            continue
        for key, value in p.items():
            rest[key] += value
    if rest and rest != grand:
        sens = rest["tp"] / (rest["tp"] + rest["fn"])
        ppv = rest["tp"] / (rest["tp"] + rest["fp"])
        print("%-15s %5d %5d %5d %12.1f%% %12.1f%%"
              % ("all but cis", rest["tp"], rest["fn"], rest["fp"],
                 100 * sens, 100 * ppv))
        print()
        print("all but cis  = excluding cargo_is_same, the architecture where the")
        print("               flanking IS also occurs inside the cargo and the graph")
        print("               genuinely cannot distinguish the element from decoys")
    return 0


if __name__ == "__main__":
    sys.exit(main())
