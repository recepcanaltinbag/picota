#!/usr/bin/env python3
"""
run_copy_sweep.py
-----------------
How does recovery decay as the flanking IS becomes more repetitive?

The scenario table varies structure; this varies exactly one number. Every
element carries its OWN insertion sequence, and that sequence is planted in the
chromosome a controlled number of additional times, so copy number moves while
how many composite transposons share an IS stays fixed at one. `--shared-is`
cannot do this: raising its copy count also raises how many elements collapse
onto the same node, and the two effects are then inseparable.

Both assemblers run on every genome, and both routes are measured on each --
the assembly graph that PICOTA reads, and the contigs a reader would otherwise
inspect by hand. That is the comparison the sweep exists to make: contigs and
graphs do not decay at the same rate.

Copy levels are 2 (the element itself, nothing else), 4, 8, 16 and 32. Two seeds
per level, twenty elements each, so every level rests on forty elements rather
than a handful.

Usage:
  python scripts/run_copy_sweep.py --out-dir sweep/ --backbone ref/mg1655.fna \\
      --art /path/to/art_illumina --threads 8
"""

import argparse
import os
import shutil
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PICOTA_DIR = os.path.join(SCRIPT_DIR, "..", "picota")
sys.path.insert(0, SCRIPT_DIR)

from run_scenarios import assemble, run  # noqa: E402

# 5 kb rather than the 20 kb the scenarios use: at 32 copies twenty elements
# need 620 insertion points, and 20 kb spacing would want 12 Mb of chromosome.
# 5 kb is still an order of magnitude above the 350 bp fragment size, so the
# copies stay independently placed.
SPACING = 5000

COPY_LEVELS = (2, 4, 8, 16, 32)


def cases(levels, seeds, n_cts):
    for copies in levels:
        for seed in seeds:
            yield ("copies%02d_s%d" % (copies, seed), copies, seed, n_cts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--backbone", required=True)
    parser.add_argument("--n-cts", type=int, default=20)
    parser.add_argument("--levels", type=int, nargs="+", default=list(COPY_LEVELS))
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--coverage", type=float, default=50.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-cycle-size", type=int, default=1000)
    parser.add_argument("--assembler", nargs="+", default=["spades", "megahit"],
                        choices=("spades", "megahit"))
    parser.add_argument("--art", default="art_illumina")
    parser.add_argument("--spades", default="spades.py")
    parser.add_argument("--megahit", default="megahit")
    parser.add_argument("--megahit-toolkit", default="megahit_core")
    parser.add_argument("--fastg2gfa",
                        default=os.path.join(PICOTA_DIR, "tools", "gfaview",
                                             "misc", "fastg2gfa"))
    parser.add_argument("--only", nargs="+")
    args = parser.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    sys.path.insert(0, os.path.abspath(PICOTA_DIR))
    from src.cycle_finderv2 import cycle_analysis  # noqa: E402

    for name, copies, seed, n_cts in cases(args.levels, args.seeds, args.n_cts):
        if args.only and name not in args.only:
            continue
        case = os.path.join(args.out_dir, name)
        if all(os.path.exists(os.path.join(case, "cycles_%s.fasta" % a))
               for a in args.assembler):
            print("[skip] %s" % name, file=sys.stderr)
            continue
        os.makedirs(case, exist_ok=True)
        log = os.path.join(case, "run.log")
        print("[run ] %s  (%d copies, %d elements)" % (name, copies, n_cts),
              file=sys.stderr, flush=True)

        run([sys.executable, os.path.join(SCRIPT_DIR, "simulate_ct_genome.py"),
             "--out-dir", case, "--n-cts", str(n_cts),
             "--shared-is", "0",
             "--is-copies-per-element", str(copies - 2),
             "--is-divergence", "0.5", "--spacing", str(SPACING),
             "--backbone-fasta", args.backbone, "--seed", str(seed)], log)

        prefix = os.path.join(case, "art_")
        run([args.art, "-ss", "HSXt", "-i", os.path.join(case, "genome.fasta"),
             "-p", "-l", "150", "-f", str(args.coverage), "-m", "350", "-s", "50",
             "-rs", str(seed), "-o", prefix, "-na"], log)

        for assembler in args.assembler:
            gfa = assemble(args, case, prefix, assembler, log)
            cycles = os.path.join(case, "cycles_%s.fasta" % assembler)
            cycle_analysis(gfa, cycles, True, 25, args.min_cycle_size, 40000,
                           "Cycle", 1, 25, 80, 80, dedup_mode="strict")
            if assembler == "spades":
                shutil.copyfile(cycles, os.path.join(case, "cycles.fasta"))
            print("       %s / %-7s %d cycles" % (name, assembler, sum(
                1 for line in open(cycles) if line.startswith(">"))),
                file=sys.stderr, flush=True)

        for leftover in (prefix + "1.fq", prefix + "2.fq"):
            if os.path.exists(leftover):
                os.remove(leftover)

    print("\nSweep in %s" % args.out_dir, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
