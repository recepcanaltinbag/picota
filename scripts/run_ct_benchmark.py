#!/usr/bin/env python3
"""
run_ct_benchmark.py
-------------------
Drive the whole PICOTA benchmark: simulate genomes with known composite
transposons, sequence them, assemble, run detection, and score
(docs/VALIDATION.md, phase 0.5).

One case is one (backbone, CT count, seed) combination. For each it runs

    simulate_ct_genome.py -> wgsim -> SPAdes -> cycle_analysis -> score

and appends a row to summary.tsv. Sweeping backbones and CT counts is the point:
a single genome at a single density says little, and a recall number that holds
across two host chromosomes and 3 to 12 implanted CTs says a great deal more.

Cases are independent and resumable -- a case whose summary row already exists
is skipped, so an interrupted sweep can simply be re-run.

Usage:
  python scripts/run_ct_benchmark.py \\
      --backbone ref/mg1655.fna ref/sakai.fna \\
      --n-cts 3 6 12 --seeds 1 2 \\
      --coverage 50 --out-dir bench/

Requires wgsim, spades.py, blastn and makeblastdb on PATH; override with the
--wgsim / --spades / --blastn / --makeblastdb flags (e.g. to point at a conda
environment).
"""

import argparse
import csv
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PICOTA_DIR = os.path.join(SCRIPT_DIR, "..", "picota")

SUMMARY_COLUMNS = [
    "Case", "Backbone", "NCTs", "SharedIS", "Seed", "Mode", "GenomeLength",
    "Segments", "Links", "ReportedCycles", "TruncatedSearches",
    "CTRecall", "CTTotal", "Precision", "PrecisionTotal",
    "CopyDistinct", "CopyDistinctTotal",
]


def run(command, log_path=None, check=True):
    """Run a command, tee-ing output to a log file when one is given."""
    with open(log_path, "a") if log_path else open(os.devnull, "w") as log:
        return subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                              check=check)


def genome_length(fasta):
    total = 0
    with open(fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def count_gfa(path):
    segments = links = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith("S"):
                segments += 1
            elif line.startswith("L"):
                links += 1
    return segments, links


def simulate_case(args, case_dir, backbone, n_cts, shared_is, seed):
    log = os.path.join(case_dir, "simulate.log")
    run([sys.executable, os.path.join(SCRIPT_DIR, "simulate_ct_genome.py"),
         "--out-dir", case_dir, "--backbone-fasta", backbone,
         "--n-cts", str(n_cts), "--shared-is", str(shared_is),
         "--is-copies-outside", str(args.is_copies_outside),
         "--is-divergence", str(args.is_divergence),
         "--cargo-genes", str(args.cargo_genes),
         "--spacing", str(args.spacing), "--seed", str(seed),
         "--is-fasta", args.is_fasta, "--cargo-fasta", args.cargo_fasta], log)


def sequence_and_assemble(args, case_dir):
    genome = os.path.join(case_dir, "genome.fasta")
    length = genome_length(genome)
    # wgsim counts read PAIRS, and each pair contributes 2 x read length.
    pairs = int(args.coverage * length / (2 * args.read_length))
    reads = [os.path.join(case_dir, "r1.fq"), os.path.join(case_dir, "r2.fq")]
    log = os.path.join(case_dir, "reads.log")
    run([args.wgsim, "-N", str(pairs),
         "-1", str(args.read_length), "-2", str(args.read_length),
         "-e", str(args.error_rate), "-r", "0", "-R", "0", "-X", "0",
         "-S", "1", genome] + reads, log)

    assembly_dir = os.path.join(case_dir, "spades")
    run([args.spades, "-1", reads[0], "-2", reads[1], "-o", assembly_dir,
         "-k", args.kmers, "--only-assembler", "-t", str(args.threads),
         "-m", str(args.memory_gb)], os.path.join(case_dir, "spades.log"))

    gfa = os.path.join(assembly_dir, "assembly_graph_with_scaffolds.gfa")
    if not os.path.exists(gfa):
        raise RuntimeError(f"SPAdes produced no GFA in {assembly_dir}")
    return gfa, length


def detect_cycles(args, gfa, out_fasta, mode):
    """Run cycle_analysis in-process and report how many searches were cut."""
    sys.path.insert(0, os.path.abspath(PICOTA_DIR))
    from src.cycle_finderv2 import GraphWork, cycle_analysis  # noqa: E402

    cycle_analysis(gfa, out_fasta, True, args.path_limit,
                   args.min_cycle, args.max_cycle, "Cycle", 1,
                   args.max_components, args.k_mer_sim, args.threshold_sim,
                   dedup_mode=mode)

    graph_work = GraphWork()
    graph_work.find_all_path = True
    graph_work.path_limit = args.path_limit
    node_dict, edge_dict = graph_work.parse_gfa(gfa)
    graph_work.dfs_iterative(graph_work.generate_genome_graph(node_dict, edge_dict))

    reported = sum(1 for line in open(out_fasta) if line.startswith(">"))
    return reported, graph_work.truncated_searches


def score_case(args, case_dir, cycles_fasta):
    sys.path.insert(0, SCRIPT_DIR)
    from score_picota_benchmark import read_ground_truth, run_blast, score  # noqa: E402

    ground_truth = read_ground_truth(os.path.join(case_dir, "ground_truth.tsv"))
    cycle_ids = [line[1:].strip() for line in open(cycles_fasta)
                 if line.startswith(">")]
    rows = run_blast(cycles_fasta,
                     os.path.join(case_dir, "ground_truth_cts.fasta"),
                     args.blastn, args.makeblastdb)
    return score(ground_truth, rows, cycle_ids, args.min_identity, args.min_coverage)


def existing_cases(summary_path):
    if not os.path.exists(summary_path):
        return set()
    with open(summary_path) as handle:
        return {(row["Case"], row["Mode"]) for row in csv.DictReader(handle, delimiter="\t")}


def append_row(summary_path, row):
    is_new = not os.path.exists(summary_path)
    with open(summary_path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t")
        if is_new:
            writer.writeheader()
        writer.writerow(row)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Sweep the PICOTA benchmark over backbones and CT counts.")
    parser.add_argument("--backbone", nargs="+", required=True,
                        help="Host chromosome FASTA(s) to implant into.")
    parser.add_argument("--n-cts", nargs="+", type=int, default=[3, 6, 12])
    parser.add_argument("--seeds", nargs="+", type=int, default=[1])
    parser.add_argument("--shared-fraction", type=float, default=0.5,
                        help="Fraction of each case's CTs built on ONE shared IS "
                             "element. Default: %(default)s")
    parser.add_argument("--modes", nargs="+", default=["legacy", "strict"])
    parser.add_argument("--out-dir", required=True)

    parser.add_argument("--coverage", type=float, default=50.0)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--error-rate", type=float, default=0.001)
    parser.add_argument("--kmers", default="55,77,99")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-gb", type=int, default=24)

    parser.add_argument("--is-copies-outside", type=int, default=6)
    parser.add_argument("--is-divergence", type=float, default=0.5)
    parser.add_argument("--cargo-genes", type=int, default=2)
    parser.add_argument("--spacing", type=int, default=20000)
    parser.add_argument("--is-fasta", default=os.path.join(PICOTA_DIR, "DBs", "ISes", "IS.fna"))
    parser.add_argument("--cargo-fasta",
                        default=os.path.join(PICOTA_DIR, "DBs", "Antibiotics",
                                             "nucleotide_fasta_protein_homolog_model.fasta"))

    parser.add_argument("--path-limit", type=int, default=25)
    parser.add_argument("--min-cycle", type=int, default=2000)
    parser.add_argument("--max-cycle", type=int, default=40000)
    parser.add_argument("--max-components", type=int, default=25)
    parser.add_argument("--k-mer-sim", type=int, default=80)
    parser.add_argument("--threshold-sim", type=int, default=80)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--min-coverage", type=float, default=0.95)

    parser.add_argument("--wgsim", default="wgsim")
    parser.add_argument("--spades", default="spades.py")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    args = parser.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    summary_path = os.path.join(args.out_dir, "summary.tsv")
    done = existing_cases(summary_path)

    for backbone in args.backbone:
        label = os.path.splitext(os.path.basename(backbone))[0]
        for n_cts in args.n_cts:
            shared = max(2, int(round(n_cts * args.shared_fraction)))
            shared = min(shared, n_cts)
            for seed in args.seeds:
                case = f"{label}_n{n_cts}_s{seed}"
                if all((case, mode) in done for mode in args.modes):
                    print(f"[skip] {case}", file=sys.stderr)
                    continue

                case_dir = os.path.join(args.out_dir, case)
                os.makedirs(case_dir, exist_ok=True)
                print(f"[run ] {case}: {n_cts} CTs, {shared} sharing one IS",
                      file=sys.stderr)

                try:
                    if not os.path.exists(os.path.join(case_dir, "ground_truth.tsv")):
                        simulate_case(args, case_dir, backbone, n_cts, shared, seed)
                    gfa = os.path.join(case_dir, "spades",
                                       "assembly_graph_with_scaffolds.gfa")
                    if os.path.exists(gfa):
                        length = genome_length(os.path.join(case_dir, "genome.fasta"))
                    else:
                        gfa, length = sequence_and_assemble(args, case_dir)
                    segments, links = count_gfa(gfa)
                except Exception as error:  # noqa: BLE001 - one bad case must not stop the sweep
                    print(f"[fail] {case}: {error}", file=sys.stderr)
                    continue

                for mode in args.modes:
                    if (case, mode) in done:
                        continue
                    cycles = os.path.join(case_dir, f"cycles_{mode}.fasta")
                    try:
                        reported, truncated = detect_cycles(args, gfa, cycles, mode)
                        result = score_case(args, case_dir, cycles)
                    except Exception as error:  # noqa: BLE001
                        print(f"[fail] {case}/{mode}: {error}", file=sys.stderr)
                        continue

                    recall, recall_total = result["ct_recall"]
                    precision, precision_total = result["precision"]
                    distinct, distinct_total = result["copy_distinctness"]
                    append_row(summary_path, {
                        "Case": case, "Backbone": label, "NCTs": n_cts,
                        "SharedIS": shared, "Seed": seed, "Mode": mode,
                        "GenomeLength": length, "Segments": segments, "Links": links,
                        "ReportedCycles": reported, "TruncatedSearches": truncated,
                        "CTRecall": recall, "CTTotal": recall_total,
                        "Precision": precision, "PrecisionTotal": precision_total,
                        "CopyDistinct": distinct, "CopyDistinctTotal": distinct_total,
                    })
                    print(f"       {mode}: recall {recall}/{recall_total}, "
                          f"precision {precision}/{precision_total}, "
                          f"copy-distinct {distinct}/{distinct_total}", file=sys.stderr)

    print(f"\nSummary: {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
