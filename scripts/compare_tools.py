#!/usr/bin/env python3
"""
compare_tools.py
----------------
PICOTA against MobileElementFinder, on the same genomes and the same criterion.

MobileElementFinder (Johansson et al. 2021) predicts composite transposons from
assembled sequence and is the closest published comparator: it reports the same
structure PICOTA reports, so the two can be scored against one ground truth
without either being handicapped by a translation between output formats.

Three inputs per case, which is what makes the comparison fair rather than
flattering:

  genome    the true chromosome -- a perfect assembly. This is what
            MobileElementFinder can do when nothing is lost upstream, and it is
            the number that must be quoted when saying the tool "missed"
            anything. Without it a reader cannot tell a tool limitation from an
            assembly limitation.
  spades    SPAdes contigs from simulated reads
  megahit   MEGAHIT contigs from the same reads

PICOTA reads the assembly graph produced from those same reads, so every route
sees identical input data.

Matching is deliberately identical to the PICOTA benchmark: a predicted element
counts only when one predicted sequence covers >= --min-coverage of a ground
truth element at >= 95% identity. Counting any BLAST hit instead once inflated a
precision figure in this benchmark by treating a 78% fragment as a recovery.

Composite transposon predictions are taken from the `cn_*` records in
MobileElementFinder's mge_sequences.fna, which are the ones its own log calls
putative composite transposons.

Usage:
  python scripts/compare_tools.py --runs copysweep/ --mefinder path/to/mefinder \\
      --threads 2 --jobs 6
"""

import argparse
import collections
import concurrent.futures
import glob
import os
import shutil
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from score_picota_benchmark import (  # noqa: E402
    covered_positions,
    read_ground_truth,
    run_blast,
)

# Input name -> path within the case directory.
INPUTS = {
    "genome": "genome.fasta",
    "spades": os.path.join("spades", "contigs.fasta"),
    "megahit": os.path.join("megahit", "final.contigs.fa"),
}


def mefinder_predictions(case, which, mefinder, threads, force=False):
    """
    Path to the composite transposon predictions for one input, running
    MobileElementFinder if they are not cached.

    Its temp directory is set per run: the default is a fixed /tmp/mge_finder,
    which two concurrent runs would share and overwrite.
    """
    source = os.path.join(case, INPUTS[which])
    if not os.path.exists(source):
        return None

    out_dir = os.path.join(case, "mefinder_%s" % which)
    prefix = os.path.join(out_dir, "mef")
    sequences = prefix + "_mge_sequences.fna"
    if os.path.exists(sequences) and not force:
        return sequences

    os.makedirs(out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="mef_") as tmp:
        result = subprocess.run(
            [mefinder, "find", "-c", source, "-t", str(threads),
             "--temp-dir", tmp, prefix],
            capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(sequences):
        with open(os.path.join(out_dir, "mefinder.log"), "w") as handle:
            handle.write(result.stdout or "")
            handle.write(result.stderr or "")
        return None
    return sequences


def composite_only(sequences, out_path):
    """
    Keep the putative composite transposon records.

    MobileElementFinder emits every mobile element it finds; the composite
    transposon predictions are the `cn_` records. Scoring the whole file would
    credit it for finding the flanking insertion sequences, which is a different
    claim from finding the composite transposon.
    """
    kept, keep = 0, False
    with open(sequences) as handle, open(out_path, "w") as out:
        for line in handle:
            if line.startswith(">"):
                keep = line[1:].startswith("cn_")
                kept += keep
            if keep:
                out.write(line)
    return kept


def recovered(predictions, ground_truth_fasta, ground_truth, min_coverage,
              blastn, makeblastdb):
    if not os.path.exists(predictions) or os.path.getsize(predictions) == 0:
        return set()
    coverage = covered_positions(
        run_blast(predictions, ground_truth_fasta, blastn, makeblastdb), 95.0)
    found = set()
    for ct_id, per_prediction in coverage.items():
        length = int(ground_truth[ct_id]["CT_Length"])
        if any(len(p) / length >= min_coverage for p in per_prediction.values()):
            found.add(ct_id)
    return found


def picota_recovered(case, ground_truth, min_coverage, blastn, makeblastdb):
    cycles = os.path.join(case, "cycles_spades.fasta")
    if not os.path.exists(cycles):
        return None
    return recovered(cycles, os.path.join(case, "ground_truth_cts.fasta"),
                     ground_truth, min_coverage, blastn, makeblastdb)


def run_case(case, args):
    tsv = os.path.join(case, "ground_truth.tsv")
    if not os.path.exists(tsv):
        return None
    ground_truth = read_ground_truth(tsv)
    if not ground_truth:
        return None
    ground_truth_fasta = os.path.join(case, "ground_truth_cts.fasta")
    copies = int(next(iter(ground_truth.values()))["IS_Genome_Copies"])

    row = {"case": os.path.basename(case), "copies": copies,
           "n": len(ground_truth), "counts": {}, "predicted": {}}

    for which in args.inputs:
        sequences = mefinder_predictions(case, which, args.mefinder,
                                         args.threads, args.force)
        if sequences is None:
            continue
        composite = os.path.join(case, "mefinder_%s" % which, "composite.fna")
        predicted = composite_only(sequences, composite)
        found = recovered(composite, ground_truth_fasta, ground_truth,
                          args.min_coverage, args.blastn, args.makeblastdb)
        row["counts"]["mef_" + which] = len(found)
        row["predicted"]["mef_" + which] = predicted

    found = picota_recovered(case, ground_truth, args.min_coverage,
                             args.blastn, args.makeblastdb)
    if found is not None:
        row["counts"]["picota_graph"] = len(found)
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--mefinder", required=True)
    parser.add_argument("--inputs", nargs="+", default=["genome", "spades", "megahit"],
                        choices=sorted(INPUTS))
    parser.add_argument("--min-coverage", type=float, default=0.95)
    parser.add_argument("--threads", type=int, default=2,
                        help="Threads per MobileElementFinder run.")
    parser.add_argument("--jobs", type=int, default=4,
                        help="Concurrent cases.")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--force", action="store_true",
                        help="Re-run MobileElementFinder even when cached.")
    parser.add_argument("--per-case", action="store_true")
    args = parser.parse_args(argv)

    cases = []
    for root in args.runs:
        cases += [c for c in sorted(glob.glob(os.path.join(root, "*")))
                  if os.path.isdir(c)
                  and os.path.exists(os.path.join(c, "ground_truth.tsv"))]
    if not cases:
        print("no cases under %s" % ", ".join(args.runs), file=sys.stderr)
        return 1

    routes = ["picota_graph"] + ["mef_" + i for i in args.inputs]
    totals = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    predicted_totals = collections.defaultdict(lambda: collections.Counter())

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for row in pool.map(lambda c: run_case(c, args), cases):
            if row is None:
                continue
            missing = [r for r in routes if r not in row["counts"]]
            if missing:
                print("[skip] %s: missing %s" % (row["case"], ", ".join(missing)),
                      file=sys.stderr)
                continue
            for route in routes:
                totals[row["copies"]][route][0] += row["counts"][route]
                totals[row["copies"]][route][1] += row["n"]
            for route, count in row["predicted"].items():
                predicted_totals[row["copies"]][route] += count
            if args.per_case:
                print("%-16s copies %3d  n %2d  %s"
                      % (row["case"], row["copies"], row["n"],
                         "  ".join("%s %d" % (r, row["counts"][r]) for r in routes)))

    header = "%8s %5s  %s" % ("copies", "n",
                              "  ".join("%17s" % r for r in routes))
    print()
    print(header)
    print("-" * len(header))
    for copies in sorted(totals):
        cells = []
        for route in routes:
            got, total = totals[copies][route]
            cells.append("%17s" % ("%d/%d (%3.0f%%)" % (got, total, 100 * got / total)
                                   if total else "-"))
        n = max(v[1] for v in totals[copies].values())
        print("%8d %5d  %s" % (copies, n, "  ".join(cells)))

    print()
    print("Composite transposons MobileElementFinder predicted (any, matched or not)")
    for copies in sorted(predicted_totals):
        print("  copies %3d  %s" % (copies, "  ".join(
            "%s %d" % (r, n) for r, n in sorted(predicted_totals[copies].items()))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
