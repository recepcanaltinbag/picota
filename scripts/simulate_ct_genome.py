#!/usr/bin/env python3
"""
simulate_ct_genome.py
---------------------
Build a synthetic bacterial genome containing a known set of composite
transposons, for measuring PICOTA's recall against exact ground truth
(docs/VALIDATION.md, phase 0.5).

Downloading closed genomes gives realism but never certainty: the true CT
catalogue of a real isolate is itself an annotation, with its own errors. A
simulated genome inverts that trade -- every composite transposon is known by
construction, so recall and precision are exact rather than estimated. The
sequences themselves are real: IS elements come from the bundled ISfinder set
and cargo genes from CARD, so what is synthetic is the arrangement, not the
biology.

The arrangement is what the benchmark is about. The case PICOTA has never been
tested on is several composite transposons built on the SAME IS element with
DIFFERENT cargo: the assembler collapses every copy of that IS into one node,
and each cargo becomes a separate branch through it. --shared-is sets how many
CTs share one element, which is the knob that matters.

Usage:
  python scripts/simulate_ct_genome.py --out-dir sim/ \
      --backbone-length 4500000 --n-cts 8 --shared-is 4 \
      --is-copies-outside 6 --is-divergence 0.5 --seed 1

Outputs, all in --out-dir:
  genome.fasta            the simulated chromosome
  ground_truth.tsv        one row per composite transposon
  ground_truth_cts.fasta  the CT sequences, for identity comparison
  ground_truth.json       the same, machine-readable, plus run parameters
"""

import argparse
import json
import os
import random
import sys

DNA = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

# Elements outside this band are not insertion sequences in the useful sense:
# the short end is fragments, the long end is whole transposons.
IS_MIN_LENGTH = 700
IS_MAX_LENGTH = 2500

GROUND_TRUTH_COLUMNS = [
    "CT_ID", "IS_Name", "IS_Family", "IS_Length", "IS_Divergence_Pct",
    "IS_Orientation", "Cargo_Genes", "Cargo_Length", "CT_Start", "CT_End",
    "CT_Length", "IS_Genome_Copies",
]


def reverse_complement(seq):
    return "".join(COMPLEMENT.get(base, "N") for base in reversed(seq))


def read_fasta(path, min_length=0, max_length=None):
    """Yield (header, sequence) pairs, filtered by length."""
    header, chunks = None, []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(chunks)
                    if len(seq) >= min_length and (max_length is None or len(seq) <= max_length):
                        yield header, seq
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip().upper())
    if header is not None:
        seq = "".join(chunks)
        if len(seq) >= min_length and (max_length is None or len(seq) <= max_length):
            yield header, seq


def is_clean_dna(seq):
    """Reject sequences with ambiguity codes -- they confuse assembly and BLAST."""
    return bool(seq) and set(seq) <= set(DNA)


def random_backbone(length, gc_content, rng):
    """
    Filler sequence with the requested GC content.

    Random filler is deliberate: it contains no repeats, so every repeat in the
    finished genome is one this script put there. That is what makes the ground
    truth exact. It also makes the genome easier to assemble than a real one, so
    recall measured here is an upper bound -- state it that way in any write-up.
    """
    gc = "GC"
    at = "AT"
    return "".join(rng.choice(gc) if rng.random() < gc_content else rng.choice(at)
                   for _ in range(length))


def mutate(seq, divergence_pct, rng):
    """Apply independent substitutions at the given percentage."""
    if divergence_pct <= 0:
        return seq
    rate = divergence_pct / 100.0
    out = []
    for base in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in DNA if b != base]))
        else:
            out.append(base)
    return "".join(out)


def parse_is_name(header):
    """ISfinder headers are ISname_family_group."""
    parts = header.split("_")
    if len(parts) >= 3:
        return parts[0], parts[-2]
    return header, "unknown"


def parse_card_name(header):
    """CARD headers are gb|ACC|strand|coords|ARO:id|GeneName [organism]."""
    fields = header.split("|")
    if len(fields) >= 6:
        return fields[5].split("[")[0].strip()
    return header.split()[0]


def load_is_elements(path, rng, count, min_length=IS_MIN_LENGTH,
                     max_length=IS_MAX_LENGTH):
    """Pick distinct IS elements of usable length, deterministically."""
    pool = [(parse_is_name(h), seq)
            for h, seq in read_fasta(path, min_length, max_length)
            if is_clean_dna(seq)]
    if len(pool) < count:
        raise SystemExit(f"Only {len(pool)} usable IS elements in {path}, need {count}")
    return rng.sample(pool, count)


def load_cargo_genes(path, rng, count):
    """Pick distinct AMR genes to serve as cargo."""
    pool = [(parse_card_name(h), seq)
            for h, seq in read_fasta(path, 300, 4000)
            if is_clean_dna(seq)]
    if len(pool) < count:
        raise SystemExit(f"Only {len(pool)} usable cargo genes in {path}, need {count}")
    return rng.sample(pool, count)


def build_composite_transposon(is_seq, cargo_seqs, divergence_pct, orientation, rng):
    """
    IS - cargo - IS, the defining structure.

    Each flanking copy is mutated independently, because real copies of one IS
    within a genome are near-identical rather than identical, and that small
    divergence is exactly what a k-mer identity threshold has to cope with.
    """
    left = mutate(is_seq, divergence_pct, rng)
    right = mutate(is_seq, divergence_pct, rng)
    if orientation == "inverted":
        right = reverse_complement(right)

    spacer = lambda: "".join(rng.choice(DNA) for _ in range(rng.randint(50, 200)))
    cargo = spacer().join(cargo_seqs) if len(cargo_seqs) > 1 else cargo_seqs[0]
    return left + spacer() + cargo + spacer() + right, cargo


def simulate(args):
    rng = random.Random(args.seed)

    if args.shared_is > args.n_cts:
        raise SystemExit("--shared-is cannot exceed --n-cts")

    # One IS element is reused by --shared-is composite transposons; the rest
    # each get their own.
    distinct_is_needed = 1 + (args.n_cts - args.shared_is) if args.shared_is else args.n_cts
    is_elements = load_is_elements(args.is_fasta, rng, distinct_is_needed,
                                   args.is_min_length, args.is_max_length)
    cargo_pool = load_cargo_genes(args.cargo_fasta, rng, args.n_cts * args.cargo_genes)

    assignments = []
    shared = is_elements[0] if args.shared_is else None
    others = iter(is_elements[1:] if args.shared_is else is_elements)
    for index in range(args.n_cts):
        element = shared if index < args.shared_is else next(others)
        cargo = cargo_pool[index * args.cargo_genes:(index + 1) * args.cargo_genes]
        assignments.append((element, cargo))

    # Count how many copies of each IS end up in the genome: two per composite
    # transposon that uses it, plus any free-standing copies requested.
    is_copy_counts = {}
    for (name_family, _), _ in assignments:
        is_copy_counts[name_family[0]] = is_copy_counts.get(name_family[0], 0) + 2
    if shared:
        is_copy_counts[shared[0][0]] = is_copy_counts.get(shared[0][0], 0) + args.is_copies_outside

    pieces = []
    records = []
    position = 0

    def append(seq):
        nonlocal position
        pieces.append(seq)
        start = position
        position += len(seq)
        return start

    for index, ((name_family, is_seq), cargo) in enumerate(assignments, start=1):
        append(random_backbone(args.spacing, args.gc_content, rng))
        ct_seq, cargo_seq = build_composite_transposon(
            is_seq, [seq for _, seq in cargo], args.is_divergence,
            args.is_orientation, rng)
        start = append(ct_seq)
        is_name, is_family = name_family
        records.append({
            "CT_ID": f"CT{index:03d}",
            "IS_Name": is_name,
            "IS_Family": is_family,
            "IS_Length": len(is_seq),
            "IS_Divergence_Pct": args.is_divergence,
            "IS_Orientation": args.is_orientation,
            "Cargo_Genes": ";".join(name for name, _ in cargo),
            "Cargo_Length": len(cargo_seq),
            "CT_Start": start + 1,
            "CT_End": start + len(ct_seq),
            "CT_Length": len(ct_seq),
            "IS_Genome_Copies": is_copy_counts[is_name],
            "Sequence": ct_seq,
        })

    # Free-standing copies of the shared IS, so its assembly-graph node carries
    # the depth of a genuinely multi-copy element rather than of a pair.
    if shared and args.is_copies_outside:
        for _ in range(args.is_copies_outside):
            append(random_backbone(args.spacing, args.gc_content, rng))
            append(mutate(shared[1], args.is_divergence, rng))

    append(random_backbone(max(args.backbone_length - position, args.spacing),
                           args.gc_content, rng))

    return "".join(pieces), records


def write_outputs(out_dir, genome, records, params):
    os.makedirs(out_dir, exist_ok=True)

    genome_path = os.path.join(out_dir, "genome.fasta")
    with open(genome_path, "w") as handle:
        handle.write(">simulated_chromosome\n")
        for i in range(0, len(genome), 70):
            handle.write(genome[i:i + 70] + "\n")

    tsv_path = os.path.join(out_dir, "ground_truth.tsv")
    with open(tsv_path, "w") as handle:
        handle.write("\t".join(GROUND_TRUTH_COLUMNS) + "\n")
        for record in records:
            handle.write("\t".join(str(record[c]) for c in GROUND_TRUTH_COLUMNS) + "\n")

    cts_path = os.path.join(out_dir, "ground_truth_cts.fasta")
    with open(cts_path, "w") as handle:
        for record in records:
            handle.write(f">{record['CT_ID']}_{record['IS_Name']}_len{record['CT_Length']}\n")
            seq = record["Sequence"]
            for i in range(0, len(seq), 70):
                handle.write(seq[i:i + 70] + "\n")

    json_path = os.path.join(out_dir, "ground_truth.json")
    with open(json_path, "w") as handle:
        json.dump({"parameters": params,
                   "genome_length": len(genome),
                   "composite_transposons": records}, handle, indent=2)

    return genome_path, tsv_path, cts_path, json_path


def build_parser():
    parser = argparse.ArgumentParser(
        description="Simulate a bacterial genome with a known set of composite "
                    "transposons (docs/VALIDATION.md phase 0.5).")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--is-fasta", default="picota/DBs/ISes/IS.fna",
                        help="ISfinder nucleotide FASTA. Default: %(default)s")
    parser.add_argument("--cargo-fasta",
                        default="picota/DBs/Antibiotics/nucleotide_fasta_protein_homolog_model.fasta",
                        help="CARD nucleotide FASTA. Default: %(default)s")
    parser.add_argument("--backbone-length", type=int, default=4500000,
                        help="Approximate genome length. Default: %(default)s")
    parser.add_argument("--n-cts", type=int, default=8,
                        help="Composite transposons to implant. Default: %(default)s")
    parser.add_argument("--shared-is", type=int, default=4,
                        help="How many of them are built on the SAME IS element, "
                             "each with different cargo. This is the case PICOTA "
                             "has never been tested on. Default: %(default)s")
    parser.add_argument("--is-copies-outside", type=int, default=6,
                        help="Extra free-standing copies of the shared IS. "
                             "Default: %(default)s")
    parser.add_argument("--is-divergence", type=float, default=0.5,
                        help="Percent divergence applied to each IS copy "
                             "independently. Default: %(default)s")
    parser.add_argument("--is-orientation", choices=("direct", "inverted"),
                        default="direct",
                        help="Orientation of the flanking copies. Default: %(default)s")
    parser.add_argument("--is-min-length", type=int, default=IS_MIN_LENGTH,
                        help="Shortest IS element to draw. Raising this together "
                             "with a small --cargo-genes is what makes the shared "
                             "repeat dominate each cycle, which is the regime "
                             "where deduplication decides recall. Default: %(default)s")
    parser.add_argument("--is-max-length", type=int, default=IS_MAX_LENGTH,
                        help="Longest IS element to draw. Default: %(default)s")
    parser.add_argument("--cargo-genes", type=int, default=2,
                        help="AMR genes per cargo. Default: %(default)s")
    parser.add_argument("--spacing", type=int, default=20000,
                        help="Backbone between implanted elements. Default: %(default)s")
    parser.add_argument("--gc-content", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=1)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    genome, records = simulate(args)
    paths = write_outputs(args.out_dir, genome, records, vars(args))

    print(f"Genome: {len(genome):,} bp with {len(records)} composite transposons",
          file=sys.stderr)
    shared = [r for r in records if r["IS_Name"] == records[0]["IS_Name"]] if records else []
    if args.shared_is:
        print(f"  {args.shared_is} of them share {records[0]['IS_Name']} "
              f"({records[0]['IS_Family']}), present in "
              f"{records[0]['IS_Genome_Copies']} genome copies", file=sys.stderr)
        print(f"  distinct cargo across those: "
              f"{len({r['Cargo_Genes'] for r in shared})}", file=sys.stderr)
    for path in paths:
        print(f"  wrote {path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
