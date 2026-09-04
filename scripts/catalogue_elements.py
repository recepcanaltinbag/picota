#!/usr/bin/env python3
"""
catalogue_elements.py
---------------------
One table of every element implanted across every benchmark case.

Each case records its own elements in ground_truth.tsv, which is what scoring
reads. A reader assessing the benchmark needs the other view: which insertion
sequences were drawn, how long they were, which resistance genes they carried,
and how those choices are distributed across the whole design. Without it the
IS length range is a claim in the methods rather than something anyone can
check, and "sampled from ISfinder, 700-2500 bp" cannot be distinguished from a
run that happened to draw eight short ones.

Writes elements.tsv -- one row per implanted element, prefixed with the case it
belongs to and the axis that case varies -- and prints a summary of the
distributions that the methods section states.

Usage:
  python scripts/catalogue_elements.py --runs scenarios40/ copysweep/ \\
      --out supplementary/elements.tsv
"""

import argparse
import collections
import csv
import glob
import hashlib
import json
import os
import statistics
import sys

# Prefixed to every row so the table stands on its own once cases are pooled.
CONTEXT_COLUMNS = ["Case", "Source", "Axis"]

LENGTH_BINS = [
    ("short (<1000 bp)", 0, 1000),
    ("medium (1000-1499 bp)", 1000, 1500),
    ("long (1500-1999 bp)", 1500, 2000),
    ("very long (>=2000 bp)", 2000, 10 ** 9),
]


def axis_of(case_dir, params):
    """
    What this case varies, in the form the tables index on.

    Read from the recorded parameters rather than parsed out of the directory
    name: a name is a label someone typed, and the two have already disagreed
    once in this benchmark.
    """
    if params.get("is_copies_per_element"):
        return "copies=%d" % (params["is_copies_per_element"] + 2)
    if params.get("shared_is"):
        return "shared_is=%d" % params["shared_is"]
    if params.get("cargo_is_mode", "none") != "none":
        return "cargo_is=%s" % params["cargo_is_mode"]
    if params.get("novel_cts"):
        return "novel_cargo=%d" % params["novel_cts"]
    return "baseline"


def read_cases(roots):
    """
    Cases across every run directory, skipping genomes catalogued twice.

    Two run directories can hold the identical genome: the scenario sweep's
    baseline and the copy sweep's two-copy level are generated from the same
    seed with the same parameters, so their FASTA files are byte-identical and
    their elements are the same elements. Counting both inflated the catalogue
    from 520 elements to 560 and double-weighted their IS families in every
    distribution below.

    Identity is decided on the genome's checksum rather than on the directory
    name, because the names carry no indication that the two coincide.
    """
    seen = {}
    for root in roots:
        source = os.path.basename(os.path.normpath(root))
        for case in sorted(glob.glob(os.path.join(root, "*"))):
            tsv = os.path.join(case, "ground_truth.tsv")
            if not os.path.isfile(tsv):
                continue
            genome = os.path.join(case, "genome.fasta")
            if os.path.exists(genome):
                with open(genome, "rb") as handle:
                    digest = hashlib.md5(handle.read()).hexdigest()
                if digest in seen:
                    print("[skip] %s/%s: same genome as %s"
                          % (source, os.path.basename(case), seen[digest]),
                          file=sys.stderr)
                    continue
                seen[digest] = "%s/%s" % (source, os.path.basename(case))
            params = {}
            payload_path = os.path.join(case, "ground_truth.json")
            if os.path.exists(payload_path):
                params = json.load(open(payload_path)).get("parameters", {})
            with open(tsv) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            yield source, os.path.basename(case), axis_of(case, params), rows


def summarise(rows):
    lengths = [int(r["IS_Length"]) for r in rows if r.get("IS_Length")]
    cargo_lengths = [int(r["Cargo_Length"]) for r in rows if r.get("Cargo_Length")]
    ct_lengths = [int(r["CT_Length"]) for r in rows if r.get("CT_Length")]

    genes = collections.Counter()
    for row in rows:
        for gene in (row.get("Cargo_Genes") or "").split(";"):
            if gene:
                genes[gene] += 1

    print("Elements catalogued: %d across %d case(s)"
          % (len(rows), len({(r["Source"], r["Case"]) for r in rows})))
    print()

    print("Insertion sequences")
    print("  distinct elements   %d" % len({r["IS_Name"] for r in rows}))
    print("  distinct families   %d" % len({r["IS_Family"] for r in rows
                                            if r.get("IS_Family")}))
    if lengths:
        print("  length              %d-%d bp (median %d)"
              % (min(lengths), max(lengths), statistics.median(lengths)))
        for label, low, high in LENGTH_BINS:
            count = sum(1 for value in lengths if low <= value < high)
            if count:
                print("    %-24s %4d  (%4.1f%%)"
                      % (label, count, 100 * count / len(lengths)))
    families = collections.Counter(r["IS_Family"] for r in rows if r.get("IS_Family"))
    print("  commonest families  %s"
          % ", ".join("%s (%d)" % (f, n) for f, n in families.most_common(6)))
    print()

    print("Cargo")
    print("  distinct genes      %d" % len(genes))
    types = collections.Counter(r.get("Cargo_Type", "?") for r in rows)
    print("  types               %s"
          % ", ".join("%s %d" % (t, n) for t, n in sorted(types.items())))
    if cargo_lengths:
        print("  cargo length        %d-%d bp (median %d)"
              % (min(cargo_lengths), max(cargo_lengths),
                 statistics.median(cargo_lengths)))
    print("  commonest genes     %s"
          % ", ".join("%s (%d)" % (g, n) for g, n in genes.most_common(6)))
    print()

    print("Composite transposons")
    if ct_lengths:
        print("  element length      %d-%d bp (median %d)"
              % (min(ct_lengths), max(ct_lengths), statistics.median(ct_lengths)))
    copies = collections.Counter(r.get("IS_Genome_Copies", "?") for r in rows)
    print("  IS genome copies    %s"
          % ", ".join("%s:%d" % (c, n)
                      for c, n in sorted(copies.items(),
                                         key=lambda kv: int(kv[0]))))
    orientations = collections.Counter(r.get("IS_Orientation", "?") for r in rows)
    print("  orientation         %s"
          % ", ".join("%s %d" % (o, n) for o, n in sorted(orientations.items())))
    divergence = {r.get("IS_Divergence_Pct") for r in rows}
    print("  flank divergence    %s%%" % ", ".join(sorted(d for d in divergence if d)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--out", help="Write elements.tsv here.")
    args = parser.parse_args(argv)

    rows = []
    columns = None
    for source, case, axis, case_rows in read_cases(args.runs):
        for row in case_rows:
            if columns is None:
                columns = CONTEXT_COLUMNS + list(row.keys())
            enriched = {"Case": case, "Source": source, "Axis": axis}
            enriched.update(row)
            rows.append(enriched)

    if not rows:
        print("no ground_truth.tsv found under %s" % ", ".join(args.runs),
              file=sys.stderr)
        return 1

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        print("wrote %s (%d rows)" % (args.out, len(rows)))
        print()

    summarise(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
