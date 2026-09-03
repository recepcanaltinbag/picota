#!/usr/bin/env python3
"""
score_copy_sweep.py
-------------------
Recovery as a function of how repetitive the flanking IS is.

Reports four routes on the same elements, so the comparison is not confounded
by anything except the route itself:

  contig  spades / megahit   an element counted only when ONE contig carries
                             >= --min-coverage of it in a single alignment
  graph   spades / megahit   an element counted when a detected cycle covers
                             >= --min-coverage of it

The contig test is deliberately the strict one described in
compare_contig_recovery.py: summing several alignments from one contig would let
a contig holding only IS+cargo appear to cover an IS-cargo-IS element while
missing an entire IS.

Detection only. Scoring is a separate stage and is measured by
score_scenarios.py; mixing them here would hide which stage lost an element.

Usage:
  python scripts/score_copy_sweep.py --sweep sweep/
"""

import argparse
import collections
import glob
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from compare_contig_recovery import CONTIG_FILES, longest_contiguous_hit  # noqa: E402
from score_picota_benchmark import (  # noqa: E402
    covered_positions,
    read_ground_truth,
    run_blast,
)

ASSEMBLERS = ("spades", "megahit")


def graph_recovered(case, assembler, ground_truth, min_coverage, blastn,
                    makeblastdb):
    """Elements a detected cycle covers to at least min_coverage."""
    cycles = os.path.join(case, "cycles_%s.fasta" % assembler)
    if not os.path.exists(cycles):
        return None
    ground_truth_fasta = os.path.join(case, "ground_truth_cts.fasta")
    coverage = covered_positions(
        run_blast(cycles, ground_truth_fasta, blastn, makeblastdb), 95.0)
    found = set()
    for ct_id, per_cycle in coverage.items():
        length = int(ground_truth[ct_id]["CT_Length"])
        if any(len(p) / length >= min_coverage for p in per_cycle.values()):
            found.add(ct_id)
    return found


def contig_recovered(case, assembler, ground_truth, min_coverage, min_identity,
                     blastn, makeblastdb):
    contigs = os.path.join(case, CONTIG_FILES[assembler])
    if not os.path.exists(contigs):
        return None
    best = longest_contiguous_hit(
        contigs, os.path.join(case, "ground_truth_cts.fasta"),
        min_identity, blastn, makeblastdb)
    return {ct for ct, record in ground_truth.items()
            if best.get(ct, 0) / int(record["CT_Length"]) >= min_coverage}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--sweep", required=True)
    parser.add_argument("--min-coverage", type=float, default=0.95)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--per-case", action="store_true")
    args = parser.parse_args(argv)

    # copies -> route -> [recovered, total]
    totals = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))

    for case in sorted(glob.glob(os.path.join(args.sweep, "*"))):
        tsv = os.path.join(case, "ground_truth.tsv")
        if not os.path.exists(tsv):
            continue
        ground_truth = read_ground_truth(tsv)
        if not ground_truth:
            continue
        copies = int(next(iter(ground_truth.values()))["IS_Genome_Copies"])
        n = len(ground_truth)

        row = {}
        for assembler in ASSEMBLERS:
            for route, fn in (("graph", graph_recovered),
                              ("contig", contig_recovered)):
                if route == "graph":
                    found = fn(case, assembler, ground_truth,
                               args.min_coverage, args.blastn, args.makeblastdb)
                else:
                    found = fn(case, assembler, ground_truth,
                               args.min_coverage, args.min_identity,
                               args.blastn, args.makeblastdb)
                if found is None:
                    continue
                key = "%s_%s" % (route, assembler)
                row[key] = len(found)
                totals[copies][key][0] += len(found)
                totals[copies][key][1] += n

        if args.per_case:
            print("%-16s copies %3d  n %3d  %s"
                  % (os.path.basename(case), copies, n,
                     "  ".join("%s %d" % (k, v) for k, v in sorted(row.items()))))

    routes = ["graph_spades", "graph_megahit", "contig_spades", "contig_megahit"]
    header = "%8s %6s  %s" % ("copies", "n", "  ".join("%16s" % r for r in routes))
    print()
    print(header)
    print("-" * len(header))
    for copies in sorted(totals):
        counts = totals[copies]
        n = max(v[1] for v in counts.values()) if counts else 0
        cells = []
        for route in routes:
            got, total = counts.get(route, [0, 0])
            cells.append("%16s" % ("%d/%d (%3.0f%%)" % (got, total,
                                                        100 * got / total)
                                   if total else "-"))
        print("%8d %6d  %s" % (copies, n, "  ".join(cells)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
