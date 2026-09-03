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
    "IS_Orientation", "Cargo_Genes", "Cargo_Type", "Cargo_IS", "Cargo_Repeats",
    "Cargo_Length", "CT_Start",
    "CT_End", "CT_Length", "IS_Genome_Copies",
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


def load_backbone(path):
    """
    Read a real genome to implant into, taking its longest record.

    A real chromosome carries its own IS elements -- K-12 MG1655 alone has
    dozens -- so the implanted CTs are exact ground truth but the background is
    not empty. Native IS pairs can form cycles of their own, which show up as
    apparent false positives. That is a property of the test, not of PICOTA, and
    any report using a real backbone has to say so.
    """
    records = list(read_fasta(path))
    if not records:
        raise SystemExit(f"No sequences in {path}")
    header, seq = max(records, key=lambda record: len(record[1]))
    return header.split()[0], seq


def choose_positions(backbone_length, count, min_spacing, rng, placement):
    """
    Pick `count` insertion points at least `min_spacing` apart.

    'random' draws uniformly over all valid configurations: sample offsets in
    the length left over once the mandatory gaps are reserved, sort them, then
    add the gaps back. Even spacing is a systematic artefact -- it puts every
    composite transposon the same distance apart, which is not how a genome
    looks and could flatter or punish an assembler in ways nothing else would
    reveal -- so it is available but not the default.
    """
    if not count:
        return []

    reserved = (count + 1) * min_spacing
    available = backbone_length - reserved
    if available < 0:
        raise SystemExit(
            f"Backbone of {backbone_length:,} bp cannot hold {count} inserts at "
            f"--spacing {min_spacing:,}; use a longer backbone or fewer elements")

    if placement == "even":
        span = backbone_length // (count + 1)
        return [span * (i + 1) for i in range(count)]

    offsets = sorted(rng.randint(0, available) for _ in range(count))
    return [min_spacing * (i + 1) + offset for i, offset in enumerate(offsets)]


def place_inserts(backbone, inserts, rng, min_spacing, placement="random"):
    """
    Splice inserts into the backbone.

    Returns (genome, starts) where starts[i] is the 0-based offset of inserts[i]
    in the finished genome. Positions are chosen against the ORIGINAL backbone
    and the offsets accumulated as we splice, because every insert shifts
    everything downstream -- getting that wrong is the one way to produce ground
    truth that looks plausible and is silently wrong.
    """
    count = len(inserts)
    if not count:
        return backbone, []

    positions = choose_positions(len(backbone), count, min_spacing, rng, placement)

    pieces = []
    starts = []
    output_length = 0
    previous = 0
    for position, sequence in zip(positions, inserts):
        chunk = backbone[previous:position]
        pieces.append(chunk)
        output_length += len(chunk)
        starts.append(output_length)
        pieces.append(sequence)
        output_length += len(sequence)
        previous = position
    pieces.append(backbone[previous:])
    return "".join(pieces), starts


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


STOP_CODONS = {"TAA", "TAG", "TGA"}


def make_orf(length, rng):
    """
    A random open reading frame: start codon, no internal stops, stop codon.

    This is cargo for a *novel* composite transposon -- one whose payload is not
    an antibiotic resistance or xenobiotic gene and therefore has no homology to
    any database PICOTA searches. Prodigal will still call it as a gene, so the
    candidate reaches scoring exactly as a real novel element would; what it
    cannot do is collect homology points. Testing with only CARD cargo would
    measure PICOTA on the easy half of its own problem.
    """
    codons = ["ATG"]
    for _ in range(max(length // 3 - 2, 1)):
        while True:
            codon = "".join(rng.choice(DNA) for _ in range(3))
            if codon not in STOP_CODONS:
                codons.append(codon)
                break
    codons.append(rng.choice(sorted(STOP_CODONS)))
    return "".join(codons)


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


def build_composite_transposon(is_seq, cargo_seqs, divergence_pct, orientation,
                               rng, cargo_repeats=(), cargo_is=None):
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

    # Repeats embedded inside the cargo are what fragment a composite transposon
    # in the assembly graph. Each one also occurs elsewhere in the genome, so the
    # assembler cannot resolve it and breaks the cargo into separate nodes --
    # raising the cycle's component count, which the score penalises through
    # sqrt(|comp - 2|). Real cargo carries exactly this kind of thing: integron
    # cassettes, small mobile elements, duplicated segments.
    pieces = []
    for index, gene in enumerate(cargo_seqs):
        if index:
            pieces.append(spacer())
        pieces.append(gene)
        if index < len(cargo_repeats):
            pieces.append(spacer())
            pieces.append(cargo_repeats[index])
    # An IS landing inside the cargo, which happens constantly in real elements.
    # Where it goes in the graph depends entirely on which element it is: a
    # different IS becomes its own node and merely lengthens the cycle, while a
    # copy of the flanking IS collapses onto the very node that defines the
    # bubble, so one graph node now sits at three places in the structure.
    if cargo_is is not None and len(pieces) > 1:
        middle = len(pieces) // 2
        pieces[middle:middle] = [spacer(), cargo_is, spacer()]
    cargo = "".join(pieces)
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

    if args.novel_cts > args.n_cts:
        raise SystemExit("--novel-cts cannot exceed --n-cts")

    assignments = []
    shared = is_elements[0] if args.shared_is else None
    others = iter(is_elements[1:] if args.shared_is else is_elements)
    for index in range(args.n_cts):
        element = shared if index < args.shared_is else next(others)
        # The last --novel-cts elements carry cargo with no database homology.
        if index >= args.n_cts - args.novel_cts:
            cargo = [("novel_orf_%d_%d" % (index + 1, g + 1),
                      make_orf(rng.randint(600, 1400), rng))
                     for g in range(args.cargo_genes)]
        else:
            cargo = cargo_pool[index * args.cargo_genes:(index + 1) * args.cargo_genes]
        assignments.append((element, cargo))

    # Count how many copies of each IS end up in the genome: two per composite
    # transposon that uses it, plus any free-standing copies requested.
    # Two flanking copies per composite transposon, plus one more inside the
    # cargo when the interrupting element is the flanking one, plus any
    # free-standing copies. Omitting the cargo-internal copy understates the
    # count by one per element, which is a ground-truth error rather than a
    # cosmetic one: it is what a reader would compare a depth ratio against.
    per_element = 3 if args.cargo_is_mode == "same" else 2
    is_copy_counts = {}
    for (name_family, _), _ in assignments:
        is_copy_counts[name_family[0]] = is_copy_counts.get(name_family[0], 0) + per_element
    if shared:
        is_copy_counts[shared[0][0]] = is_copy_counts.get(shared[0][0], 0) + args.is_copies_outside
    # Decoy copies of every element's own IS. --is-copies-outside raises the
    # copy number of one shared element and so confounds copy number with how
    # many composite transposons share it; this raises each element's own IS
    # independently, which is what isolating copy number as a variable needs.
    if args.is_copies_per_element:
        for name in {nf[0] for nf, _ in (a[0] for a in assignments)}:
            is_copy_counts[name] += args.is_copies_per_element

    # A single extra element used to interrupt every cargo in "different" mode,
    # so it too is multi-copy across the genome. It must not be one of the
    # flanking elements: an independent draw can collide, and a collision turns
    # this scenario silently into the "same" one for that composite transposon.
    interrupting_is = None
    if args.cargo_is_mode == "different":
        flanking = {name for (name, _), _ in is_elements}
        candidates = load_is_elements(args.is_fasta, rng, distinct_is_needed + 1,
                                      args.is_min_length, args.is_max_length)
        interrupting_is = next((c for c in candidates if c[0][0] not in flanking),
                               None)
        if interrupting_is is None:
            raise SystemExit("could not draw an IS element distinct from the "
                             "flanking ones; widen --is-min-length/--is-max-length")

    # One pool of repeat elements shared by every composite transposon, so each
    # is genuinely multi-copy across the genome rather than unique to one cargo.
    repeat_pool = [make_dna(args.cargo_repeat_length, rng.randint(0, 10 ** 6))
                   for _ in range(args.cargo_repeats)]

    inserts = []
    pending = []
    for index, ((name_family, is_seq), cargo) in enumerate(assignments, start=1):
        cargo_is = None
        if args.cargo_is_mode == "same":
            cargo_is = mutate(is_seq, args.is_divergence, rng)
        elif args.cargo_is_mode == "different":
            cargo_is = mutate(interrupting_is[1], args.is_divergence, rng)

        ct_seq, cargo_seq = build_composite_transposon(
            is_seq, [seq for _, seq in cargo], args.is_divergence,
            args.is_orientation, rng, repeat_pool, cargo_is)
        is_name, is_family = name_family
        inserts.append(ct_seq)
        pending.append({
            "CT_ID": f"CT{index:03d}",
            "IS_Name": is_name,
            "IS_Family": is_family,
            "IS_Length": len(is_seq),
            "IS_Divergence_Pct": args.is_divergence,
            "IS_Orientation": args.is_orientation,
            "Cargo_Genes": ";".join(name for name, _ in cargo),
            "Cargo_Type": "novel" if cargo[0][0].startswith("novel_orf_") else "AMR",
            "Cargo_Repeats": args.cargo_repeats,
            "Cargo_IS": (is_name if args.cargo_is_mode == "same"
                         else interrupting_is[0][0] if interrupting_is else "none"),
            "Cargo_Length": len(cargo_seq),
            "CT_Length": len(ct_seq),
            "IS_Genome_Copies": is_copy_counts[is_name],
            "Sequence": ct_seq,
        })

    # Free-standing copies of the shared IS, so its assembly-graph node carries
    # the depth of a genuinely multi-copy element rather than of a pair.
    if shared and args.is_copies_outside:
        for _ in range(args.is_copies_outside):
            inserts.append(mutate(shared[1], args.is_divergence, rng))

    # Free-standing copies of each element's own IS, mutated independently so
    # they diverge from the flanking pair exactly as those diverge from one
    # another. Placed once per distinct element, not once per composite
    # transposon, so a shared IS is not counted twice.
    if args.is_copies_per_element:
        seen = set()
        for (name_family, is_seq), _ in assignments:
            if name_family[0] in seen:
                continue
            seen.add(name_family[0])
            for _ in range(args.is_copies_per_element):
                inserts.append(mutate(is_seq, args.is_divergence, rng))

    # Extra copies of each cargo repeat, so the assembler really cannot place
    # them and the cargo fragments as intended.
    for repeat in repeat_pool:
        for _ in range(args.cargo_repeat_copies):
            inserts.append(repeat)

    if args.backbone_fasta:
        backbone_name, backbone = load_backbone(args.backbone_fasta)
    else:
        backbone_name = "simulated_chromosome"
        needed = (len(pending) * (1 + args.is_copies_per_element)
                  + args.is_copies_outside + 1) * args.spacing
        backbone = random_backbone(max(args.backbone_length, needed),
                                   args.gc_content, rng)

    genome, starts = place_inserts(backbone, inserts, rng, args.spacing,
                                   args.placement)

    records = []
    for record, start in zip(pending, starts):
        record["CT_Start"] = start + 1
        record["CT_End"] = start + record["CT_Length"]
        records.append(record)

    return genome, records, backbone_name


def write_outputs(out_dir, genome, records, params, backbone_name="simulated_chromosome"):
    os.makedirs(out_dir, exist_ok=True)

    genome_path = os.path.join(out_dir, "genome.fasta")
    with open(genome_path, "w") as handle:
        handle.write(f">{backbone_name}_with_{len(records)}_implanted_CTs\n")
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
                   "backbone": backbone_name,
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
    parser.add_argument("--backbone-fasta",
                        help="Real genome to implant into, instead of random "
                             "filler. Its longest record is used. Note that a "
                             "real chromosome carries its own IS elements, which "
                             "can form cycles of their own and appear as false "
                             "positives -- say so when reporting.")
    parser.add_argument("--backbone-length", type=int, default=4500000,
                        help="Length of the random backbone, when "
                             "--backbone-fasta is not given. Default: %(default)s")
    parser.add_argument("--n-cts", type=int, default=8,
                        help="Composite transposons to implant. Default: %(default)s")
    parser.add_argument("--shared-is", type=int, default=4,
                        help="How many of them are built on the SAME IS element, "
                             "each with different cargo. This is the case PICOTA "
                             "has never been tested on. Default: %(default)s")
    parser.add_argument("--is-copies-outside", type=int, default=6,
                        help="Extra free-standing copies of the shared IS. "
                             "Default: %(default)s")
    parser.add_argument("--is-copies-per-element", type=int, default=0,
                        help="Free-standing copies of EACH element's own IS, so "
                             "every element sits at 2 + N genome copies "
                             "(3 + N under --cargo-is-mode same). Unlike "
                             "--is-copies-outside, which raises one shared "
                             "element, this varies copy number without varying "
                             "how many composite transposons share an IS.")
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
    parser.add_argument("--cargo-is-mode", choices=("none", "different", "same"),
                        default="none",
                        help="Insert an IS element into each cargo. 'different' "
                             "uses another element, which becomes its own graph "
                             "node. 'same' uses a copy of the flanking IS, which "
                             "collapses onto the node that defines the bubble and "
                             "is the harder case by far. Default: %(default)s")
    parser.add_argument("--cargo-repeats", type=int, default=0,
                        help="Repeat elements embedded inside each cargo. Each "
                             "one fragments the composite transposon in the "
                             "assembly graph, raising the cycle's component "
                             "count, which the score penalises. Use this to test "
                             "candidates that do not assemble into a clean "
                             "two-node bubble. Default: %(default)s")
    parser.add_argument("--cargo-repeat-length", type=int, default=400,
                        help="Length of each embedded repeat. Default: %(default)s")
    parser.add_argument("--cargo-repeat-copies", type=int, default=3,
                        help="Free-standing copies of each embedded repeat "
                             "elsewhere in the genome, which is what makes it "
                             "unresolvable. Default: %(default)s")
    parser.add_argument("--novel-cts", type=int, default=0,
                        help="How many of the composite transposons carry cargo "
                             "with NO database homology -- random open reading "
                             "frames rather than CARD genes. These are the novel "
                             "elements PICOTA claims to find, and they cannot "
                             "collect homology score. Default: %(default)s")
    parser.add_argument("--cargo-genes", type=int, default=2,
                        help="AMR genes per cargo. Default: %(default)s")
    parser.add_argument("--placement", choices=("random", "even"), default="random",
                        help="Where to implant. 'random' scatters the elements, "
                             "subject to --spacing; 'even' puts them at equal "
                             "intervals, which is a systematic artefact. "
                             "Default: %(default)s")
    parser.add_argument("--spacing", type=int, default=20000,
                        help="Minimum backbone between implanted elements. "
                             "Default: %(default)s")
    parser.add_argument("--gc-content", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=1)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    genome, records, backbone_name = simulate(args)
    paths = write_outputs(args.out_dir, genome, records, vars(args), backbone_name)

    print(f"Genome: {len(genome):,} bp ({backbone_name}) with {len(records)} "
          f"composite transposons", file=sys.stderr)
    if args.backbone_fasta:
        print("  backbone is a real genome: its native IS elements can form "
              "cycles of their own, so apparent false positives are expected",
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
