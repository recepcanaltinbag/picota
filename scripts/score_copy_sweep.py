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
    """
    (elements a cycle covers, cycles covering no element, cycles detected).

    A cycle counts for an element only when it covers min_coverage of it. The
    weaker test -- any BLAST hit between the two -- once inflated a precision
    figure here: two candidates were 100% composite-transposon sequence but
    spanned 78% of the element, and counting them as recoveries reported more
    true positives than there were implanted elements.

    The unmatched count matters on this axis in particular. Free-standing IS
    copies give the graph a high-depth node with many edges, and cycles closing
    through it are exactly what a detector reading bubbles would be expected to
    over-report as copy number rises.
    """
    cycles = os.path.join(case, "cycles_%s.fasta" % assembler)
    if not os.path.exists(cycles):
        return None
    detected = {line[1:].strip() for line in open(cycles) if line.startswith(">")}
    ground_truth_fasta = os.path.join(case, "ground_truth_cts.fasta")
    coverage = covered_positions(
        run_blast(cycles, ground_truth_fasta, blastn, makeblastdb), 95.0)

    found, matched_cycles = set(), set()
    for ct_id, per_cycle in coverage.items():
        length = int(ground_truth[ct_id]["CT_Length"])
        for cycle, positions in per_cycle.items():
            if len(positions) / length >= min_coverage:
                found.add(ct_id)
                matched_cycles.add(cycle)
    return found, detected - matched_cycles, detected


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
    parser.add_argument("--include-incomplete", action="store_true",
                        help="Aggregate cases that are missing a route. Off by "
                             "default: a half-finished case otherwise enters "
                             "the totals with a smaller denominator for some "
                             "routes than others, and the percentages stop "
                             "being comparable across the row.")
    args = parser.parse_args(argv)

    # copies -> route -> [recovered, total]
    totals = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    # copies -> assembler -> [cycles matching no element, cycles detected]
    unmatched_totals = collections.defaultdict(
        lambda: collections.defaultdict(lambda: [0, 0]))

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
        found_by_route = {}
        unmatched_by_assembler = {}
        for assembler in ASSEMBLERS:
            for route, fn in (("graph", graph_recovered),
                              ("contig", contig_recovered)):
                if route == "graph":
                    result = fn(case, assembler, ground_truth,
                                args.min_coverage, args.blastn, args.makeblastdb)
                    if result is not None:
                        found, unmatched, detected = result
                        unmatched_by_assembler[assembler] = (len(unmatched),
                                                             len(detected))
                    else:
                        found = None
                else:
                    found = fn(case, assembler, ground_truth,
                               args.min_coverage, args.min_identity,
                               args.blastn, args.makeblastdb)
                key = "%s_%s" % (route, assembler)
                found_by_route[key] = found
                if found is not None:
                    row[key] = len(found)

        missing = [k for k, v in found_by_route.items() if v is None]
        if missing and not args.include_incomplete:
            print("[skip] %s: still missing %s"
                  % (os.path.basename(case), ", ".join(sorted(missing))),
                  file=sys.stderr)
            continue

        for key, found in found_by_route.items():
            if found is None:
                continue
            totals[copies][key][0] += len(found)
            totals[copies][key][1] += n
        for assembler, (unmatched, detected) in unmatched_by_assembler.items():
            unmatched_totals[copies][assembler][0] += unmatched
            unmatched_totals[copies][assembler][1] += detected

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
        denominators = {v[1] for v in counts.values()}
        n = max(denominators) if denominators else 0
        if len(denominators) > 1:
            print("[warn] copies=%d: routes disagree on the denominator %s"
                  % (copies, sorted(denominators)), file=sys.stderr)
        cells = []
        for route in routes:
            got, total = counts.get(route, [0, 0])
            cells.append("%16s" % ("%d/%d (%3.0f%%)" % (got, total,
                                                        100 * got / total)
                                   if total else "-"))
        print("%8d %6d  %s" % (copies, n, "  ".join(cells)))

    print()
    print("Cycles matching no implanted element (graph route)")
    header = "%8s  %-22s %-22s" % ("copies", "spades", "megahit")
    print(header)
    print("-" * len(header))
    for copies in sorted(unmatched_totals):
        cells = []
        for assembler in ASSEMBLERS:
            unmatched, detected = unmatched_totals[copies][assembler]
            cells.append("%-22s" % ("%d/%d (%3.0f%%)"
                                    % (unmatched, detected,
                                       100 * unmatched / detected)
                                    if detected else "-"))
        print("%8d  %s" % (copies, " ".join(cells)))
    print("\nDetection only: these are candidate cycles, before scoring "
          "rejects them.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
