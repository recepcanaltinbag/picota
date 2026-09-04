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
import shutil
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PICOTA_DIR = os.path.join(SCRIPT_DIR, "..", "picota")
sys.path.insert(0, SCRIPT_DIR)

from run_manifest import write_run_parameters  # noqa: E402

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


def assemble(args, case, prefix, assembler, log):
    """Assemble and return the GFA path. MEGAHIT emits no GFA of its own, so
    its largest-k contigs go through contig2fastg and fastg2gfa."""
    out = os.path.join(case, assembler)
    if assembler == "spades":
        run([args.spades, "-1", prefix + "1.fq", "-2", prefix + "2.fq",
             "-o", out, "-k", "55,77,99", "--only-assembler",
             "-t", str(args.threads), "-m", "16"], log)
        return os.path.join(out, "assembly_graph_with_scaffolds.gfa")

    run([args.megahit, "-1", prefix + "1.fq", "-2", prefix + "2.fq",
         "-o", out, "-t", str(args.threads), "--k-list", "55,77,99"], log)
    contigs = os.path.join(out, "intermediate_contigs", "k99.contigs.fa")
    if not os.path.exists(contigs):
        raise RuntimeError("MEGAHIT produced no k99 contigs in %s" % out)
    gfa = os.path.join(out, "assembly_graph.gfa")
    with open(log, "a") as handle:
        with open(gfa + ".fastg", "w") as fastg:
            subprocess.run([args.megahit_toolkit, "contig2fastg", "99", contigs],
                           stdout=fastg, stderr=handle, check=True)
        with open(gfa, "w") as out_gfa:
            subprocess.run([args.fastg2gfa, gfa + ".fastg"],
                           stdout=out_gfa, stderr=handle, check=True)
    return gfa


def run(command, log, check=True):
    with open(log, "a") as handle:
        return subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT,
                              check=check)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--backbone", help="Host genome; omit for random filler.")
    parser.add_argument("--backbone-length", type=int, default=500000)
    parser.add_argument("--is-copies-per-element", type=int, default=0,
                        help="Extra free-standing copies of each element's own "
                             "IS, so a scenario can be run at a copy number "
                             "other than the two its flanks provide. The "
                             "structural scenarios are otherwise all two-copy, "
                             "which is the easiest case for a contig method.")
    parser.add_argument("--n-cts", type=int, default=4,
                        help="Elements per genome. Kept small deliberately: "
                             "forty elements add eighty IS copies on top of the "
                             "host's own, which makes the chromosome more "
                             "repetitive than any real isolate and changes "
                             "assembly globally rather than at the element.")
    parser.add_argument("--coverage", type=float, default=50.0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(1, 11)),
                        help="One genome per seed. Elements within a genome "
                             "share its assembly and read set, so they are not "
                             "independent measurements; genomes are.")
    parser.add_argument("--art", default="art_illumina")
    parser.add_argument("--spades", default="spades.py")
    parser.add_argument("--megahit", default="megahit")
    parser.add_argument("--megahit-toolkit", default="megahit_core",
                        help="contig2fastg provider; PICOTA itself calls "
                             "megahit_core (src/assembly.py).")
    parser.add_argument("--fastg2gfa",
                        default=os.path.join(PICOTA_DIR, "tools", "gfaview",
                                             "misc", "fastg2gfa"),
                        help="Defaults to the copy PICOTA ships and config.yaml "
                             "already points at, so the benchmark converts "
                             "MEGAHIT graphs exactly as the pipeline does.")
    parser.add_argument("--assembler", nargs="+", default=["spades"],
                        choices=("spades", "megahit"),
                        help="Assemblers to run. Each writes its own "
                             "cycles fasta; spades also writes cycles.fasta.")
    parser.add_argument("--min-cycle-size", type=int, default=1000,
                        help="Passed to cycle_analysis. Must track "
                             "min_size_of_cycle in config.yaml -- this was "
                             "hardcoded to 2000 and silently held the compact "
                             "scenario at the pre-2026-09 threshold.")
    parser.add_argument("--only", nargs="+", help="Run only these scenarios.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    chosen = [s for s in scenarios(args.n_cts) if not args.only or s[0] in args.only]

    plan = [(name, extra, seed) for name, extra in chosen for seed in args.seeds]
    for name, extra, seed in plan:
        case = os.path.join(args.out_dir, "%s_s%d" % (name, seed))
        done = all(os.path.exists(os.path.join(case, "cycles_%s.fasta" % a))
                   for a in args.assembler)
        if done:
            print("[skip] %s" % name, file=sys.stderr)
            continue
        os.makedirs(case, exist_ok=True)
        log = os.path.join(case, "run.log")
        print("[run ] %s seed %d" % (name, seed), file=sys.stderr, flush=True)

        cmd = [sys.executable, os.path.join(SCRIPT_DIR, "simulate_ct_genome.py"),
               "--out-dir", case, "--n-cts", str(args.n_cts),
               "--is-copies-per-element", str(args.is_copies_per_element),
               "--is-divergence", "0.5", "--spacing", "20000",
               "--seed", str(seed)]
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
             "-rs", str(seed), "-o", prefix, "-na"], log)

        sys.path.insert(0, os.path.abspath(PICOTA_DIR))
        from src.cycle_finderv2 import cycle_analysis  # noqa: E402

        for assembler in args.assembler:
            gfa = assemble(args, case, prefix, assembler, log)
            cycles = os.path.join(case, "cycles_%s.fasta" % assembler)
            cycle_analysis(gfa, cycles, True, 25, args.min_cycle_size, 40000,
                           "Cycle", 1, 25, 80, 80, dedup_mode="strict")
            if assembler == "spades":
                shutil.copyfile(cycles, os.path.join(case, "cycles.fasta"))
                # The depth sidecar has to travel with the FASTA it describes.
                # Scoring looks for <cycles>.depths.tsv beside the file it was
                # handed, so copying only the FASTA left every case without
                # depth: depth_ratio came back None and score3's multi-copy term
                # sat at its unknown-depth default of 0.5 for every candidate,
                # silently contributing nothing to any ranking.
                depths = os.path.splitext(cycles)[0] + ".depths.tsv"
                if os.path.exists(depths):
                    shutil.copyfile(depths, os.path.join(case, "cycles.depths.tsv"))
            print("       %s s%d / %s: %d cycles" % (name, seed, assembler, sum(
                1 for line in open(cycles) if line.startswith(">"))),
                file=sys.stderr, flush=True)

        write_run_parameters(case, {
            "is_copies_per_element": args.is_copies_per_element,
            "min_cycle_size": args.min_cycle_size,
            "coverage": args.coverage,
            "assemblers": list(args.assembler),
            "kmers": "55,77,99",
            "dedup_mode": "strict",
            "read_length": 150,
            "read_simulator": "art_illumina",
            "art_profile": "HSXt",
            "scenario": name, "seed": seed, "n_cts": args.n_cts,
        }, tool_overrides={"art_illumina": args.art})

        for leftover in (prefix + "1.fq", prefix + "2.fq"):
            if os.path.exists(leftover):
                os.remove(leftover)

    print("\nScenarios in %s" % args.out_dir, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
