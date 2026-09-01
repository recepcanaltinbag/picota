#!/usr/bin/env python3
"""
summarize_benchmark.py
----------------------
Turn a benchmark summary.tsv into the performance table
(docs/VALIDATION.md, phase 0.5).

Reports, per stage and in aggregate:

  TP   ground-truth composite transposons recovered
  FN   ground-truth elements missed
  FP   reported candidates that match no ground-truth element
  Sensitivity  TP / (TP + FN)   -- of the elements that exist, how many are found
  PPV          TP / (TP + FP)   -- of the candidates reported, how many are real

Specificity is deliberately absent. It needs a count of true negatives, and
there is no enumerable set of "non-composite-transposons" in a genome that a
detector could have reported and correctly did not. Quoting one would mean
inventing a denominator; sensitivity and PPV together already say what a reader
needs, and the FP count says the rest.

Usage:
  python scripts/summarize_benchmark.py bench/summary.tsv
  python scripts/summarize_benchmark.py bench/summary.tsv --stage scored
"""

import argparse
import collections
import csv
import sys


def as_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def fmt(value):
    return "-" if value is None else f"{value * 100:.1f}%"


def load(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def case_metrics(row):
    """TP/FN/FP for one summary row. FP counts reported candidates with no match."""
    tp = as_int(row["CTRecall"])
    total = as_int(row["CTTotal"])
    reported = as_int(row["ReportedCycles"])
    explained = as_int(row["Precision"])
    return {
        "tp": tp,
        "fn": total - tp,
        "fp": max(reported - explained, 0),
        "total": total,
        "reported": reported,
        "sensitivity": ratio(tp, total),
        "ppv": ratio(explained, reported),
        "novel": (as_int(row.get("NovelRecall")), as_int(row.get("NovelTotal"))),
        "amr": (as_int(row.get("AMRRecall")), as_int(row.get("AMRTotal"))),
        "distinct": (as_int(row["CopyDistinct"]), as_int(row["CopyDistinctTotal"])),
    }


def print_per_case(rows, stage):
    header = ("%-9s %5s %6s %6s %8s %8s %4s %4s %4s %12s %8s %10s"
              % ("genome", "CTs", "novel", "seed", "detected", "reported",
                 "TP", "FN", "FP", "sensitivity", "PPV", "shared-IS"))
    print(header)
    print("-" * len(header))
    for row in rows:
        if row.get("Stage", "detection") != stage:
            continue
        m = case_metrics(row)
        print("%-9s %5s %6s %6s %8s %8s %4d %4d %4d %12s %8s %10s"
              % (row["Backbone"], row["CTTotal"], row.get("NovelTotal", "-"),
                 row["Seed"], row["ReportedCycles"] if stage == "detection" else "-",
                 m["reported"], m["tp"], m["fn"], m["fp"],
                 fmt(m["sensitivity"]), fmt(m["ppv"]),
                 "%d/%d" % m["distinct"]))


def print_aggregate(rows, stage):
    totals = collections.Counter()
    for row in rows:
        if row.get("Stage", "detection") != stage:
            continue
        m = case_metrics(row)
        totals["tp"] += m["tp"]
        totals["fn"] += m["fn"]
        totals["fp"] += m["fp"]
        totals["reported"] += m["reported"]
        totals["explained"] += as_int(row["Precision"])
        totals["novel_tp"] += m["novel"][0]
        totals["novel_n"] += m["novel"][1]
        totals["amr_tp"] += m["amr"][0]
        totals["amr_n"] += m["amr"][1]
        totals["dist_tp"] += m["distinct"][0]
        totals["dist_n"] += m["distinct"][1]
        totals["cases"] += 1

    if not totals["cases"]:
        print("  (no rows for this stage)")
        return

    sens = ratio(totals["tp"], totals["tp"] + totals["fn"])
    ppv = ratio(totals["explained"], totals["reported"])
    print("  cases                %d" % totals["cases"])
    print("  ground-truth CTs     %d" % (totals["tp"] + totals["fn"]))
    print("  candidates reported  %d" % totals["reported"])
    print("  TP / FN / FP         %d / %d / %d" % (totals["tp"], totals["fn"], totals["fp"]))
    print("  sensitivity          %s" % fmt(sens))
    print("  PPV                  %s" % fmt(ppv))
    print("  AMR-cargo recall     %d/%d (%s)"
          % (totals["amr_tp"], totals["amr_n"], fmt(ratio(totals["amr_tp"], totals["amr_n"]))))
    print("  novel-cargo recall   %d/%d (%s)"
          % (totals["novel_tp"], totals["novel_n"], fmt(ratio(totals["novel_tp"], totals["novel_n"]))))
    print("  shared-IS distinct   %d/%d (%s)"
          % (totals["dist_tp"], totals["dist_n"], fmt(ratio(totals["dist_tp"], totals["dist_n"]))))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Summarise a benchmark run.")
    parser.add_argument("summary")
    parser.add_argument("--stage", nargs="+", default=None,
                        help="Stages to report. Default: whatever the file holds.")
    args = parser.parse_args(argv)

    rows = load(args.summary)
    if not rows:
        print("empty summary", file=sys.stderr)
        return 1

    stages = args.stage or sorted({r.get("Stage", "detection") for r in rows},
                                  key=lambda s: (s != "detection", s))
    for stage in stages:
        print("\n=== stage: %s ===\n" % stage)
        print_per_case(rows, stage)
        print("\n  aggregate")
        print_aggregate(rows, stage)
    print("\nSpecificity is not reported: it requires a true-negative count, and\n"
          "there is no enumerable set of non-composite-transposons a detector\n"
          "could have reported and correctly did not.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
