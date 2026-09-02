#!/usr/bin/env python3
"""
compare_contig_recovery.py
--------------------------
How many composite transposons survive assembly into a single contig?

The control that decides whether the assembly graph is worth reading at all: if
the elements can be lifted straight off the contigs, PICOTA adds nothing. Run
against the same simulated ground truth, on the same assemblies PICOTA was given.

An element counts as recovered only when ONE contig carries at least
`--min-coverage` of it in a single contiguous alignment. Two weaker tests are
deliberately rejected:

  * summing several alignments from one contig -- a contig holding only IS+cargo
    aligns twice to an IS-cargo-IS element, once at each flanking copy, and
    appears to cover the whole thing while missing an entire IS;
  * a contig merely longer than the element -- length says nothing about what is
    inside it.

Measured on MG1655 with eight implanted elements: SPAdes carries three intact,
MEGAHIT none, while PICOTA recovers all eight from the graph. The failures are
not spread evenly -- every element built on a 14-copy IS is broken in both
assemblers, and the intact ones are all on elements present twice.

Usage:
  python scripts/compare_contig_recovery.py bench/ [--assembler spades megahit]
"""

import argparse
import collections
import glob
import os
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from score_picota_benchmark import read_ground_truth  # noqa: E402

CONTIG_FILES = {
    "spades": os.path.join("spades", "contigs.fasta"),
    "megahit": os.path.join("megahit", "final.contigs.fa"),
}


def longest_contiguous_hit(contigs, ground_truth_fasta, min_identity,
                           blastn, makeblastdb):
    """
    Longest single alignment per element, in bases.

    One alignment, not the union of several: an element is recovered when a
    contig spans it, not when pieces of it can be found scattered about.
    """
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "contigs")
        subprocess.run([makeblastdb, "-in", contigs, "-dbtype", "nucl", "-out", db],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        result = subprocess.run(
            [blastn, "-query", ground_truth_fasta, "-db", db,
             "-outfmt", "6 qseqid pident length", "-evalue", "1e-20"],
            check=True, capture_output=True, text=True)

    best = collections.defaultdict(int)
    for line in result.stdout.splitlines():
        query, pident, length = line.split("\t")
        if float(pident) < min_identity:
            continue
        ct_id = query.split("_")[0]
        best[ct_id] = max(best[ct_id], int(length))
    return best


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Count composite transposons intact in a single contig.")
    parser.add_argument("runs", help="Directory of benchmark cases.")
    parser.add_argument("--assembler", nargs="+", default=["spades"],
                        choices=sorted(CONTIG_FILES))
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--min-coverage", type=float, default=0.95)
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--per-element", action="store_true",
                        help="List every element rather than one line per case.")
    args = parser.parse_args(argv)

    totals = collections.Counter()
    for case in sorted(glob.glob(os.path.join(args.runs, "*"))):
        ground_truth_tsv = os.path.join(case, "ground_truth.tsv")
        ground_truth_fasta = os.path.join(case, "ground_truth_cts.fasta")
        if not os.path.exists(ground_truth_tsv):
            continue
        ground_truth = read_ground_truth(ground_truth_tsv)

        recovered = {}
        for assembler in args.assembler:
            contigs = os.path.join(case, CONTIG_FILES[assembler])
            if not os.path.exists(contigs):
                continue
            best = longest_contiguous_hit(contigs, ground_truth_fasta,
                                          args.min_identity, args.blastn,
                                          args.makeblastdb)
            fractions = {ct: best.get(ct, 0) / int(rec["CT_Length"])
                         for ct, rec in ground_truth.items()}
            recovered[assembler] = fractions
            intact = sum(1 for f in fractions.values() if f >= args.min_coverage)
            totals[assembler] += intact
            totals[assembler + "_total"] += len(ground_truth)

        if not recovered:
            continue

        summary = "  ".join(
            "%s %d/%d" % (a, sum(1 for f in recovered[a].values() if f >= args.min_coverage),
                          len(ground_truth))
            for a in args.assembler if a in recovered)
        print("%-18s %2d elements | intact in one contig: %s"
              % (os.path.basename(case), len(ground_truth), summary))

        if args.per_element:
            for ct in sorted(ground_truth):
                copies = ground_truth[ct].get("IS_Genome_Copies", "?")
                per = "  ".join("%s %3.0f%%" % (a, recovered[a].get(ct, 0) * 100)
                                for a in args.assembler if a in recovered)
                print("    %-8s IS %-12s copies %-4s %s"
                      % (ct, ground_truth[ct].get("IS_Name", "?")[:12], copies, per))

    print()
    for assembler in args.assembler:
        total = totals.get(assembler + "_total", 0)
        if total:
            intact = totals[assembler]
            print("%-9s intact in one contig: %d/%d (%.1f%%)"
                  % (assembler, intact, total, 100 * intact / total))
    return 0


if __name__ == "__main__":
    sys.exit(main())
