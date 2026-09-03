#!/usr/bin/env python3
"""
package_supplementary.py
------------------------
Collect the simulated genomes and their ground truth into one directory for
publication as supplementary data.

Everything a reader needs to reproduce or re-score a benchmark case, and nothing
that can be regenerated cheaply. Per case:

  <case>_genome.fasta        the simulated chromosome, implants included
  <case>_ground_truth.tsv    one row per composite transposon
  <case>_ground_truth.fasta  the element sequences, for identity comparison
  <case>_parameters.json     the exact arguments the case was built from

plus a manifest.tsv indexing every case, and README.md stating what the files
are and how to regenerate them. Reads and assemblies are deliberately excluded:
they run to gigabytes and follow from the genome plus the recorded seed.

Usage:
  python scripts/package_supplementary.py --runs bench/ scenarios/ --out supplementary/
"""

import argparse
import glob
import json
import os
import shutil
import sys

MANIFEST_COLUMNS = [
    "Case", "Source", "GenomeLength", "NumCTs", "SharedIS", "CopiesPerElement",
    "ISGenomeCopies", "NovelCargo", "CargoIS", "MinCycleSize", "Backbone", "Seed",
]


def case_dirs(roots):
    for root in roots:
        for path in sorted(glob.glob(os.path.join(root, "*"))):
            if os.path.isfile(os.path.join(path, "ground_truth.json")):
                yield root, path


def summarise(payload, run_payload=None):
    params = payload.get("parameters", {})
    run_params = (run_payload or {}).get("parameters", {})
    cts = payload.get("composite_transposons", [])
    # The realised copy number, not the requested one: cargo_is_mode "same"
    # puts a third copy inside the cargo, so 2 + N understates it by one and a
    # reader comparing a depth ratio against this column would be misled.
    copies = sorted({c.get("IS_Genome_Copies") for c in cts
                     if c.get("IS_Genome_Copies") is not None})
    return {
        "GenomeLength": payload.get("genome_length", ""),
        "NumCTs": len(cts),
        "SharedIS": params.get("shared_is", ""),
        "CopiesPerElement": params.get("is_copies_per_element", 0),
        "ISGenomeCopies": (str(copies[0]) if len(copies) == 1
                           else "%d-%d" % (copies[0], copies[-1]) if copies
                           else ""),
        "NovelCargo": sum(1 for c in cts if c.get("Cargo_Type") == "novel"),
        "CargoIS": params.get("cargo_is_mode", "none"),
        "MinCycleSize": run_params.get("min_cycle_size", ""),
        "Backbone": payload.get("backbone", ""),
        "Seed": params.get("seed", ""),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Package benchmark cases for supplement.")
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--compress", action="store_true",
                        help="Also write a .tar.gz of the packaged directory.")
    args = parser.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    rows = []

    for root, case in case_dirs(args.runs):
        name = os.path.basename(case)
        source = os.path.basename(os.path.normpath(root))
        label = "%s_%s" % (source, name)

        payload = json.load(open(os.path.join(case, "ground_truth.json")))
        for src, suffix in (("genome.fasta", "genome.fasta"),
                            ("ground_truth.tsv", "ground_truth.tsv"),
                            ("ground_truth_cts.fasta", "ground_truth.fasta")):
            source_path = os.path.join(case, src)
            if os.path.exists(source_path):
                shutil.copy2(source_path, os.path.join(args.out, "%s_%s" % (label, suffix)))
        # How the genome was built and how it was then sequenced, assembled and
        # traversed. The second half decides the result as much as the first --
        # the compact scenario reported the wrong number for weeks because the
        # harness used a cycle size the case directory did not record -- so both
        # travel with the supplement.
        run_path = os.path.join(case, "run_parameters.json")
        run_payload = json.load(open(run_path)) if os.path.exists(run_path) else {}
        with open(os.path.join(args.out, "%s_parameters.json" % label), "w") as handle:
            # Genome-build parameters stay at the top level, where readers of
            # earlier supplements already expect them; the run settings and
            # tool versions are added alongside rather than nesting everything
            # and breaking that layout.
            combined = dict(payload.get("parameters", {}))
            combined["run"] = run_payload.get("parameters", {})
            combined["tools"] = run_payload.get("tools", {})
            json.dump(combined, handle, indent=2, sort_keys=True)

        row = {"Case": label, "Source": source}
        row.update(summarise(payload, run_payload))
        rows.append(row)

    with open(os.path.join(args.out, "manifest.tsv"), "w") as handle:
        handle.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(c, "")) for c in MANIFEST_COLUMNS) + "\n")

    with open(os.path.join(args.out, "README.md"), "w") as handle:
        handle.write(
            "# Simulated benchmark data\n\n"
            "Genomes containing a known set of composite transposons, used to measure\n"
            "PICOTA's recall and precision against exact ground truth. IS elements are\n"
            "real sequences from ISfinder and cargo genes are real sequences from CARD;\n"
            "what is simulated is the arrangement, not the biology.\n\n"
            "## Files, per case\n\n"
            "| suffix | contents |\n|---|---|\n"
            "| `_genome.fasta` | the simulated chromosome, implants included |\n"
            "| `_ground_truth.tsv` | one row per composite transposon: IS element and family, cargo genes and type, coordinates, length, genome-wide IS copy number |\n"
            "| `_ground_truth.fasta` | the composite transposon sequences |\n"
            "| `_parameters.json` | the exact arguments the case was generated from |\n\n"
            "`manifest.tsv` indexes every case.\n\n"
            "## Regenerating\n\n"
            "Reads and assemblies are not included; they run to gigabytes and follow\n"
            "deterministically from the genome and the recorded seed:\n\n"
            "```\n"
            "python scripts/simulate_ct_genome.py --out-dir CASE @parameters\n"
            "art_illumina -ss HSXt -i CASE/genome.fasta -p -l 150 -f 50 -m 350 -s 50 -rs SEED -o CASE/art_ -na\n"
            "spades.py -1 CASE/art_1.fq -2 CASE/art_2.fq -o CASE/spades -k 55,77,99 --only-assembler\n"
            "```\n\n"
            "Coordinates in `_ground_truth.tsv` are 1-based and inclusive: for every row,\n"
            "`genome[CT_Start-1:CT_End]` is exactly the sequence in `_ground_truth.fasta`.\n")

    print("packaged %d case(s) into %s" % (len(rows), args.out), file=sys.stderr)

    if args.compress:
        archive = shutil.make_archive(args.out.rstrip("/"), "gztar", args.out)
        print("archive: %s" % archive, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
