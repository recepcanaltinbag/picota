#!/usr/bin/env python3
"""
score_scenarios.py
------------------
Score a graded scenario run and print the success/shortcoming table
(docs/VALIDATION.md, phase 0.5).

Takes the directory produced by run_scenarios.py, runs PICOTA's scoring stage on
each scenario's cycles, and reports for every one:

  detected     candidate cycles the graph traversal produced
  scored       candidates surviving Prodigal, BLAST and the score threshold
  TP / FN / FP against the implanted ground truth
  score range  of the elements that were recovered
  comp         graph components per recovered element, which is what the score
               penalises through sqrt(|comp - 2|)

The point of the table is the comparison down the column, not any single row:
each scenario changes one structural property and holds everything else fixed,
so a drop between rows is attributable to the structure.

Usage:
  python scripts/score_scenarios.py --scenarios scenarios/ [--threshold 50]
"""

import argparse
import csv
import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PICOTA_DIR = os.path.join(SCRIPT_DIR, "..", "picota")
sys.path.insert(0, SCRIPT_DIR)

from score_picota_benchmark import (  # noqa: E402
    covered_positions,
    read_ground_truth,
    run_blast,
)

# Order matters: the table is a gradient, read top to bottom.
SCENARIO_ORDER = ["baseline", "novel_cargo", "cargo_is_diff", "shared_is",
                  "cargo_is_same", "compact"]

SCENARIO_LABEL = {
    "baseline": "unique IS, CARD cargo",
    "novel_cargo": "cargo not in any database",
    "cargo_is_diff": "different IS inside cargo",
    "shared_is": "several CTs on one multi-copy IS",
    "cargo_is_same": "flanking IS also inside cargo",
    "compact": "cycle below min_size_of_cycle",
}


def score_cycles(case, threshold, tools):
    """Run PICOTA scoring; return {CycleID: score0} for candidates above threshold."""
    sys.path.insert(0, os.path.abspath(PICOTA_DIR))
    from src.scoringv4ProtBlast import scoring_main  # noqa: E402

    out_dir = os.path.join(case, "scoring_t%d" % threshold)
    tab = os.path.join(out_dir, "picota_final_tab")
    if not os.path.exists(tab):
        dbs = os.path.join(PICOTA_DIR, "DBs")
        scoring_main(os.path.join(case, "cycles.fasta"), out_dir,
                     os.path.join(dbs, "Antibiotics/protein_fasta_protein_homolog_model.fasta"),
                     os.path.join(dbs, "Xenobiotics/Xenobiotics_classified.fasta"),
                     os.path.join(dbs, "ISes/_tncentral_nointegrall_isfinder-TNs.fasta"),
                     os.path.join(dbs, "CompTns/Known_Tns.fasta"),
                     os.path.join(case, "blastdb_t%d" % threshold),
                     threshold_final_score=threshold, **tools)
    if not os.path.exists(tab):
        return {}
    with open(tab) as handle:
        return {row["CycleID"]: float(row["score0"])
                for row in csv.DictReader(handle, delimiter="\t")}


def components(cycle_id):
    match = re.search(r"-comp(\d+)-", cycle_id)
    return int(match.group(1)) if match else None


def evaluate(case, threshold, tools):
    ground_truth = read_ground_truth(os.path.join(case, "ground_truth.tsv"))
    cycles = os.path.join(case, "cycles.fasta")
    detected = [l[1:].strip() for l in open(cycles) if l.startswith(">")]

    coverage = covered_positions(
        run_blast(cycles, os.path.join(case, "ground_truth_cts.fasta"),
                  tools.get("path_of_blastn", "blastn"),
                  tools.get("path_of_makeblastdb", "makeblastdb")), 95.0)
    cycle_of_ct = {}
    for ct_id, per_cycle in coverage.items():
        length = int(ground_truth[ct_id]["CT_Length"])
        for cycle, positions in per_cycle.items():
            if len(positions) / length >= 0.95:
                cycle_of_ct.setdefault(ct_id, []).append(cycle)

    scores = score_cycles(case, threshold, tools)
    recovered = {ct: [c for c in cs if c in scores] for ct, cs in cycle_of_ct.items()}
    recovered = {ct: cs for ct, cs in recovered.items() if cs}

    hit_scores = [scores[c] for cs in recovered.values() for c in cs]
    comps = [components(c) for cs in recovered.values() for c in cs]
    comps = [c for c in comps if c is not None]
    real_cycles = {c for cs in recovered.values() for c in cs}

    return {
        "cts": len(ground_truth),
        "detected": len(detected),
        "detected_tp": len(cycle_of_ct),
        "scored": len(scores),
        "tp": len(recovered),
        "fn": len(ground_truth) - len(recovered),
        "fp": len(scores) - len(real_cycles),
        "score_min": min(hit_scores) if hit_scores else None,
        "score_max": max(hit_scores) if hit_scores else None,
        "comp_min": min(comps) if comps else None,
        "comp_max": max(comps) if comps else None,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Score a graded scenario run.")
    parser.add_argument("--scenarios", required=True)
    parser.add_argument("--threshold", type=float, default=50.0)
    parser.add_argument("--prodigal", default="prodigal")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--blastx", default="blastx")
    parser.add_argument("--blastp", default="blastp")
    args = parser.parse_args(argv)

    tools = {"path_of_prodigal": args.prodigal, "path_of_blastn": args.blastn,
             "path_of_makeblastdb": args.makeblastdb,
             "path_of_blastx": args.blastx, "path_of_blastp": args.blastp}

    present = [n for n in SCENARIO_ORDER
               if os.path.exists(os.path.join(args.scenarios, n, "cycles.fasta"))]
    if not present:
        print("no scored scenarios in %s" % args.scenarios, file=sys.stderr)
        return 1

    header = ("%-15s %-33s %4s %5s %5s %4s %4s %4s %13s %8s"
              % ("scenario", "structure", "CTs", "det.", "scor.", "TP", "FN", "FP",
                 "score range", "comp"))
    print(header)
    print("-" * len(header))
    for name in present:
        r = evaluate(os.path.join(args.scenarios, name), int(args.threshold), tools)
        score_range = ("%.0f-%.0f" % (r["score_min"], r["score_max"])
                       if r["score_min"] is not None else "-")
        comp_range = ("%d-%d" % (r["comp_min"], r["comp_max"])
                      if r["comp_min"] is not None else "-")
        print("%-15s %-33s %4d %5d %5d %4d %4d %4d %13s %8s"
              % (name, SCENARIO_LABEL.get(name, ""), r["cts"], r["detected"],
                 r["scored"], r["tp"], r["fn"], r["fp"], score_range, comp_range))

    print("\ndet.  = candidate cycles from graph traversal")
    print("scor. = candidates above the score threshold of %g" % args.threshold)
    print("TP/FN = implanted elements recovered / missed in the SCORED output")
    print("FP    = scored candidates matching no implanted element")
    return 0


if __name__ == "__main__":
    sys.exit(main())
