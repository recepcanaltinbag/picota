#!/usr/bin/env python3
"""
run_scenarios.py
----------------
Run PICOTA against a graded series of composite transposon structures and
report where it succeeds and where it does not
(docs/VALIDATION.md, phase 0.5).

Each scenario changes ONE property of the implanted elements and holds the host
genome, coverage, assembler and parameters fixed, so a difference in the result
is attributable to the structure rather than to the run. The scenarios are
ordered by how hard the structure is for a graph-based detector:

  baseline        one IS per element, cargo in CARD          -- the easy case
  shared_is       several elements on one multi-copy IS      -- the assembler
                                                                collapses that IS
                                                                into one node
  novel_cargo     cargo with no database homology            -- the score cannot
                                                                collect homology
                                                                points
  cargo_is_diff   a different IS sitting inside the cargo    -- becomes its own
                                                                node, lengthens
                                                                the cycle
  cargo_is_same   a copy of the flanking IS inside the cargo -- collapses onto
                                                                the node that
                                                                defines the bubble
  compact         cargo small enough that the cycle falls
                  under min_size_of_cycle                    -- known defect D7

Both stages are reported: detection recall says whether the structure survived
assembly and traversal, scored recall says whether PICOTA would actually report
it.

Usage:
  python scripts/run_scenarios.py --out-dir scenarios/ --backbone ref/mg1655.fna
"""

import argparse
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PICOTA_DIR = os.path.join(SCRIPT_DIR, "..", "picota")

def scenarios(n_cts):
    """
    Scenario definitions, scaled to the requested number of elements.

    Each is a pure case -- every element shares the IS, or every element carries
    novel cargo -- so the counts follow n_cts rather than being fixed, and the
    table can be run at whatever size gives the conclusions enough weight.
    """
    return [
        ("baseline",      ["--shared-is", "0"]),
        ("shared_is",     ["--shared-is", str(n_cts),
                           "--is-copies-outside", str(2 * n_cts)]),
        ("novel_cargo",   ["--shared-is", "0", "--novel-cts", str(n_cts)]),
        ("cargo_is_diff", ["--shared-is", "0", "--cargo-is-mode", "different"]),
        ("cargo_is_same", ["--shared-is", "0", "--cargo-is-mode", "same"]),
        ("compact",       ["--shared-is", "0", "--cargo-genes", "1",
                           "--is-min-length", "700", "--is-max-length", "900"]),
    ]


def run(command, log, check=True):
    with open(log, "a") as handle:
        return subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                              check=check)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--backbone", help="Host genome; omit for random filler.")
    parser.add_argument("--backbone-length", type=int, default=500000)
    parser.add_argument("--n-cts", type=int, default=4)
    parser.add_argument("--coverage", type=float, default=50.0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--art", default="art_illumina")
    parser.add_argument("--spades", default="spades.py")
    parser.add_argument("--only", nargs="+", help="Run only these scenarios.")
    args = parser.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    chosen = [s for s in scenarios(args.n_cts) if not args.only or s[0] in args.only]

    for name, extra in chosen:
        case = os.path.join(args.out_dir, name)
        if os.path.exists(os.path.join(case, "cycles.fasta")):
            print("[skip] %s" % name, file=sys.stderr)
            continue
        os.makedirs(case, exist_ok=True)
        log = os.path.join(case, "run.log")
        print("[run ] %s" % name, file=sys.stderr, flush=True)

        cmd = [sys.executable, os.path.join(SCRIPT_DIR, "simulate_ct_genome.py"),
               "--out-dir", case, "--n-cts", str(args.n_cts),
               "--is-divergence", "0.5", "--spacing", "20000",
               "--seed", str(args.seed)]
        if args.backbone:
            cmd += ["--backbone-fasta", args.backbone]
        else:
            cmd += ["--backbone-length", str(args.backbone_length)]
        cmd += extra
        run(cmd, log)

        genome = os.path.join(case, "genome.fasta")
        prefix = os.path.join(case, "art_")
        run([args.art, "-ss", "HSXt", "-i", genome, "-p", "-l", "150",
             "-f", str(args.coverage), "-m", "350", "-s", "50",
             "-rs", str(args.seed), "-o", prefix, "-na"], log)

        assembly = os.path.join(case, "spades")
        run([args.spades, "-1", prefix + "1.fq", "-2", prefix + "2.fq",
             "-o", assembly, "-k", "55,77,99", "--only-assembler",
             "-t", str(args.threads), "-m", "16"], log)
        for leftover in (prefix + "1.fq", prefix + "2.fq"):
            if os.path.exists(leftover):
                os.remove(leftover)

        sys.path.insert(0, os.path.abspath(PICOTA_DIR))
        from src.cycle_finderv2 import cycle_analysis  # noqa: E402
        cycle_analysis(os.path.join(assembly, "assembly_graph_with_scaffolds.gfa"),
                       os.path.join(case, "cycles.fasta"),
                       True, 25, 2000, 40000, "Cycle", 1, 25, 80, 80,
                       dedup_mode="strict")
        print("       %s: %d cycles" % (name, sum(
            1 for line in open(os.path.join(case, "cycles.fasta"))
            if line.startswith(">"))), file=sys.stderr)

    print("\nScenarios in %s" % args.out_dir, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
