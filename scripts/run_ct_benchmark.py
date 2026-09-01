#!/usr/bin/env python3
"""
run_ct_benchmark.py
-------------------
Drive the whole PICOTA benchmark: simulate genomes with known composite
transposons, sequence them, assemble, run detection, and score
(docs/VALIDATION.md, phase 0.5).

One case is one (backbone, CT count, seed) combination, scored under every
requested assembler, deduplication mode and min_size_of_cycle. For each it runs

    simulate_ct_genome.py -> ART -> SPAdes/MEGAHIT -> cycle_analysis -> score

and appends a row to summary.tsv. The sweeping is the point: a single genome, at
one CT density, from one assembler, with one seed, measures a coincidence. A
recall figure that holds across two host chromosomes, 3 to 12 implanted CTs,
several seeds and more than one assembler measures the tool -- and the assembly
graph is PICOTA's input, so a single-assembler number describes the pair rather
than PICOTA.

Cases are independent and resumable -- a case whose summary row already exists
is skipped, so an interrupted sweep can simply be re-run.

Usage:
  python scripts/run_ct_benchmark.py \\
      --backbone ref/mg1655.fna ref/sakai.fna \\
      --n-cts 3 6 12 --seeds 1 2 \\
      --coverage 50 --out-dir bench/

Reads are simulated with ART by default (empirical, platform-specific Illumina
quality profiles -- the choice a methods section can defend). --simulator wgsim
is faster but applies one uniform error rate with invented quality strings; use
it for development loops, not for published numbers.

Requires art_illumina (or wgsim), spades.py, blastn and makeblastdb on PATH;
MEGAHIT additionally needs megahit, megahit_toolkit and fastg2gfa. Override any
of them with --art / --wgsim / --spades / --megahit / --blastn / --makeblastdb.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PICOTA_DIR = os.path.join(SCRIPT_DIR, "..", "picota")

SUMMARY_COLUMNS = [
    "Case", "Backbone", "NCTs", "SharedIS", "Seed", "Assembler", "Mode",
    "MinCycle", "GenomeLength",
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
         "--is-min-length", str(args.is_min_length),
         "--is-max-length", str(args.is_max_length),
         "--is-fasta", args.is_fasta, "--cargo-fasta", args.cargo_fasta], log)


def simulate_reads(args, case_dir, genome, seed):
    """
    Generate paired-end Illumina reads.

    ART is the default because it is what a methods section can defend: it
    carries empirical, platform-specific quality profiles, so error rate rises
    along the read the way it really does. wgsim applies one uniform error rate
    with invented quality strings, which is fine for a fast development loop and
    not fine for a published benchmark.
    """
    reads = [os.path.join(case_dir, "r1.fq"), os.path.join(case_dir, "r2.fq")]
    log = os.path.join(case_dir, "reads.log")

    if args.simulator == "art":
        prefix = os.path.join(case_dir, "art_")
        run([args.art, "-ss", args.art_profile, "-i", genome, "-p",
             "-l", str(args.read_length), "-f", str(args.coverage),
             "-m", str(args.fragment_mean), "-s", str(args.fragment_sd),
             "-rs", str(seed), "-o", prefix, "-na"], log)
        for produced, wanted in zip([prefix + "1.fq", prefix + "2.fq"], reads):
            if not os.path.exists(produced):
                raise RuntimeError(f"ART produced no {produced}")
            os.replace(produced, wanted)
        return reads

    length = genome_length(genome)
    # wgsim counts read PAIRS, and each pair contributes 2 x read length.
    pairs = int(args.coverage * length / (2 * args.read_length))
    run([args.wgsim, "-N", str(pairs),
         "-1", str(args.read_length), "-2", str(args.read_length),
         "-e", str(args.error_rate), "-r", "0", "-R", "0", "-X", "0",
         "-S", str(seed), genome] + reads, log)
    return reads


def assemble(args, case_dir, reads, assembler):
    """
    Assemble and return the GFA path.

    Two assemblers rather than one because the assembly graph is PICOTA's input:
    a recall figure from a single assembler measures the pair, not the tool.
    SPAdes and MEGAHIT differ in how aggressively they simplify repeats, which
    is exactly the structure a composite transposon lives in.
    """
    assembly_dir = os.path.join(case_dir, assembler)
    log = os.path.join(case_dir, f"{assembler}.log")

    if assembler == "spades":
        run([args.spades, "-1", reads[0], "-2", reads[1], "-o", assembly_dir,
             "-k", args.kmers, "--only-assembler", "-t", str(args.threads),
             "-m", str(args.memory_gb)], log)
        gfa = os.path.join(assembly_dir, "assembly_graph_with_scaffolds.gfa")
    else:
        run([args.megahit, "-1", reads[0], "-2", reads[1], "-o", assembly_dir,
             "-t", str(args.threads), "--k-list", args.kmers.replace(",", ",")], log)
        contigs_dir = os.path.join(assembly_dir, "intermediate_contigs")
        largest_k = max(int(k) for k in args.kmers.split(","))
        fastg = os.path.join(contigs_dir, f"k{largest_k}.contigs.fa")
        if not os.path.exists(fastg):
            raise RuntimeError(f"MEGAHIT produced no k{largest_k} contigs")
        gfa = os.path.join(assembly_dir, "assembly_graph.gfa")
        with open(os.path.join(case_dir, f"{assembler}_gfa.log"), "a") as handle:
            fastg_path = gfa + ".fastg"
            subprocess.run([args.megahit_toolkit, "contig2fastg",
                            str(largest_k), fastg],
                           stdout=open(fastg_path, "w"), stderr=handle, check=True)
            subprocess.run([args.fastg2gfa, fastg_path],
                           stdout=open(gfa, "w"), stderr=handle, check=True)

    if not os.path.exists(gfa):
        raise RuntimeError(f"{assembler} produced no GFA in {assembly_dir}")
    return gfa


def sequence_and_assemble(args, case_dir, seed, assembler):
    genome = os.path.join(case_dir, "genome.fasta")
    length = genome_length(genome)
    reads = simulate_reads(args, case_dir, genome, seed)
    return assemble(args, case_dir, reads, assembler), length


def detect_cycles(args, gfa, out_fasta, mode, min_cycle):
    """Run cycle_analysis in-process and report how many searches were cut."""
    sys.path.insert(0, os.path.abspath(PICOTA_DIR))
    from src.cycle_finderv2 import GraphWork, cycle_analysis  # noqa: E402

    cycle_analysis(gfa, out_fasta, True, args.path_limit,
                   min_cycle, args.max_cycle, "Cycle", 1,
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


def free_intermediates(case_dir):
    """
    Delete the regenerable bulk once a case is scored.

    A single 5 Mb genome at 50x is roughly 250 MB of reads plus a SPAdes working
    directory, so a sweep of a few dozen cases fills a disk. Everything removed
    here is reproducible from genome.fasta and the recorded seed; the ground
    truth, the reported cycles and the summary row all stay.
    """
    for name in ("r1.fq", "r2.fq"):
        path = os.path.join(case_dir, name)
        if os.path.exists(path):
            os.remove(path)

    keep = {"assembly_graph_with_scaffolds.gfa", "assembly_graph.fastg",
            "assembly_graph.gfa", "contigs.fasta", "final.contigs.fa",
            "spades.log", "log"}
    for assembler in ("spades", "megahit"):
        assembly = os.path.join(case_dir, assembler)
        if not os.path.isdir(assembly):
            continue
        for entry in os.listdir(assembly):
            path = os.path.join(assembly, entry)
            if entry in keep:
                continue
            if os.path.isdir(path):
                shutil.rmtree(path, ignore_errors=True)
            else:
                os.remove(path)


def existing_cases(summary_path):
    if not os.path.exists(summary_path):
        return set()
    with open(summary_path) as handle:
        # Summaries written before the assembler dimension existed have no
        # Assembler column; those runs were all SPAdes.
        return {(row["Case"], row.get("Assembler") or "spades", row["Mode"],
                 int(row["MinCycle"]))
                for row in csv.DictReader(handle, delimiter="\t")}


def append_row(summary_path, row):
    """
    Append one result row, rewriting the header if the schema has grown.

    A summary written by an older version lacks columns this one produces, and
    silently appending wider rows to it would misalign every field. Rather than
    fail, the old file is migrated: existing rows keep their values and gain
    empty cells for the new columns.
    """
    if os.path.exists(summary_path):
        with open(summary_path) as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))
        current = existing[0].keys() if existing else SUMMARY_COLUMNS
        if set(current) != set(SUMMARY_COLUMNS):
            with open(summary_path, "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS,
                                        delimiter="\t", restval="")
                writer.writeheader()
                for old in existing:
                    writer.writerow({k: old.get(k, "") for k in SUMMARY_COLUMNS})
    else:
        with open(summary_path, "w", newline="") as handle:
            csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS,
                           delimiter="\t").writeheader()

    with open(summary_path, "a", newline="") as handle:
        csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t",
                       restval="").writerow(row)


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
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="Keep simulated reads and the full SPAdes working "
                             "directory. Off by default: they are the bulk of "
                             "the disk cost and are regenerable from "
                             "genome.fasta and the recorded seed.")

    parser.add_argument("--simulator", choices=("art", "wgsim"), default="art",
                        help="Read simulator. ART carries empirical "
                             "platform-specific quality profiles and is the "
                             "publication-grade choice; wgsim is faster but "
                             "applies one uniform error rate with invented "
                             "quality strings. Default: %(default)s")
    parser.add_argument("--art", default="art_illumina",
                        help="Path to art_illumina.")
    parser.add_argument("--art-profile", default="HSXt",
                        help="ART built-in sequencing system: HS25 (HiSeq 2500), "
                             "HSXt (HiSeqX TruSeq, 150bp), HSXn (HiSeqX PCR "
                             "free), MSv3 (MiSeq v3). Default: %(default)s")
    parser.add_argument("--fragment-mean", type=int, default=350,
                        help="Mean DNA fragment length for ART. Default: %(default)s")
    parser.add_argument("--fragment-sd", type=int, default=50,
                        help="Fragment length SD for ART. Default: %(default)s")
    parser.add_argument("--coverage", type=float, default=50.0)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--error-rate", type=float, default=0.001)
    parser.add_argument("--assemblers", nargs="+", default=["spades"],
                        choices=("spades", "megahit"),
                        help="Assemblers to test. The assembly graph is PICOTA's "
                             "input, so a number from one assembler measures the "
                             "pair rather than the tool. Default: %(default)s")
    parser.add_argument("--megahit", default="megahit")
    parser.add_argument("--megahit-toolkit", default="megahit_toolkit")
    parser.add_argument("--fastg2gfa",
                        default=os.path.join(PICOTA_DIR, "..", "picota", "tools",
                                             "gfaview", "misc", "fastg2gfa"))
    parser.add_argument("--kmers", default="55,77,99")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-gb", type=int, default=24)

    parser.add_argument("--is-copies-outside", type=int, default=6)
    parser.add_argument("--is-divergence", type=float, default=0.5)
    parser.add_argument("--cargo-genes", type=int, default=2)
    parser.add_argument("--spacing", type=int, default=20000)
    parser.add_argument("--is-min-length", type=int, default=700)
    parser.add_argument("--is-max-length", type=int, default=2500)
    parser.add_argument("--is-fasta", default=os.path.join(PICOTA_DIR, "DBs", "ISes", "IS.fna"))
    parser.add_argument("--cargo-fasta",
                        default=os.path.join(PICOTA_DIR, "DBs", "Antibiotics",
                                             "nucleotide_fasta_protein_homolog_model.fasta"))

    parser.add_argument("--path-limit", type=int, default=25)
    parser.add_argument("--min-cycles", nargs="+", type=int, default=[2000],
                        help="min_size_of_cycle value(s) to test. Sweeping this "
                             "measures defect D7 directly: the graph cycle of a "
                             "composite transposon is IS + cargo, so a compact "
                             "CT falls under the shipped default of 2000 and is "
                             "dropped before scoring. Default: %(default)s")
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
                if all((case, mode, mc) in done
                       for mode in args.modes for mc in args.min_cycles):
                    print(f"[skip] {case}", file=sys.stderr)
                    continue

                case_dir = os.path.join(args.out_dir, case)
                os.makedirs(case_dir, exist_ok=True)
                print(f"[run ] {case}: {n_cts} CTs, {shared} sharing one IS",
                      file=sys.stderr)

                try:
                    if not os.path.exists(os.path.join(case_dir, "ground_truth.tsv")):
                        simulate_case(args, case_dir, backbone, n_cts, shared, seed)
                    length = genome_length(os.path.join(case_dir, "genome.fasta"))
                    reads = [os.path.join(case_dir, "r1.fq"),
                             os.path.join(case_dir, "r2.fq")]
                    if not all(os.path.exists(path) for path in reads):
                        reads = simulate_reads(
                            args, case_dir,
                            os.path.join(case_dir, "genome.fasta"), seed)
                except Exception as error:  # noqa: BLE001 - one bad case must not stop the sweep
                    print(f"[fail] {case}: {error}", file=sys.stderr)
                    if not args.keep_intermediates:
                        free_intermediates(case_dir)
                    continue

                for assembler in args.assemblers:
                    try:
                        gfa = assemble(args, case_dir, reads, assembler)
                        segments, links = count_gfa(gfa)
                    except Exception as error:  # noqa: BLE001
                        print(f"[fail] {case}/{assembler}: {error}", file=sys.stderr)
                        continue

                    for mode in args.modes:
                      for min_cycle in args.min_cycles:
                        if (case, assembler, mode, min_cycle) in done:
                            continue
                        cycles = os.path.join(
                            case_dir, f"cycles_{assembler}_{mode}_{min_cycle}.fasta")
                        try:
                            reported, truncated = detect_cycles(
                                args, gfa, cycles, mode, min_cycle)
                            result = score_case(args, case_dir, cycles)
                        except Exception as error:  # noqa: BLE001
                            print(f"[fail] {case}/{assembler}/{mode}/{min_cycle}: "
                                  f"{error}", file=sys.stderr)
                            continue

                        recall, recall_total = result["ct_recall"]
                        precision, precision_total = result["precision"]
                        distinct, distinct_total = result["copy_distinctness"]
                        append_row(summary_path, {
                            "Case": case, "Backbone": label, "NCTs": n_cts,
                            "SharedIS": shared, "Seed": seed,
                            "Assembler": assembler, "Mode": mode,
                            "MinCycle": min_cycle, "GenomeLength": length,
                            "Segments": segments, "Links": links,
                            "ReportedCycles": reported,
                            "TruncatedSearches": truncated,
                            "CTRecall": recall, "CTTotal": recall_total,
                            "Precision": precision,
                            "PrecisionTotal": precision_total,
                            "CopyDistinct": distinct,
                            "CopyDistinctTotal": distinct_total,
                        })
                        print(f"       {assembler}/{mode} min_cycle={min_cycle}: "
                              f"recall {recall}/{recall_total}, "
                              f"precision {precision}/{precision_total}, "
                              f"copy-distinct {distinct}/{distinct_total}",
                              file=sys.stderr)

                if not args.keep_intermediates:
                    free_intermediates(case_dir)

    
    print(f"\nSummary: {summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
