#!/usr/bin/env python3
"""
score_picota_benchmark.py
-------------------------
Score PICOTA's reported cycles against a simulated ground truth
(docs/VALIDATION.md, phase 0.5).

Pairs with simulate_ct_genome.py: that script builds a genome whose composite
transposons are known exactly, this one measures how many of them came back.

The headline number is copy-distinctness recall. Overall recall can look fine
while the interesting failure is invisible: when several composite transposons
share one IS element the assembler collapses that element into a single node,
and a deduplication step that treats a shared node as evidence of duplication
will report one candidate where the genome holds four. Only a metric that asks
"were CTs sharing an IS reported SEPARATELY" catches that.

A reported cycle is normally shorter than the ground-truth CT it matches, and
that is correct rather than a miss: both flanking copies of the IS collapse into
one graph node, so the cycle is IS + cargo where the CT is IS + cargo + IS.
Recall is therefore measured as coverage of the ground-truth CT, not as a length
match.

Usage:
  python scripts/score_picota_benchmark.py \
      --ground-truth sim/ground_truth.tsv \
      --ground-truth-fasta sim/ground_truth_cts.fasta \
      --cycles sim/cycles_strict.fasta \
      --blastn blastn --makeblastdb makeblastdb
"""

import argparse
import collections
import os
import subprocess
import sys
import tempfile

BLAST_FORMAT = "6 qseqid sseqid pident length qlen slen qstart qend sstart send"


def read_ground_truth(path):
    """Load ground_truth.tsv into {CT_ID: {column: value}}."""
    records = {}
    header = None
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            records[fields[0]] = dict(zip(header, fields))
    return records


def run_blast(cycles_fasta, ground_truth_fasta, blastn, makeblastdb, evalue="1e-20"):
    """BLAST reported cycles against the ground-truth CTs; return raw rows."""
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "gt")
        subprocess.run([makeblastdb, "-in", ground_truth_fasta, "-dbtype", "nucl",
                        "-out", db], check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        result = subprocess.run(
            [blastn, "-query", cycles_fasta, "-db", db, "-outfmt", BLAST_FORMAT,
             "-evalue", evalue],
            check=True, capture_output=True, text=True)
    return [line.split("\t") for line in result.stdout.splitlines() if line.strip()]


def ct_id_from_subject(subject):
    """Ground-truth FASTA ids are CT001_ISname_lenNNNN."""
    return subject.split("_")[0]


def covered_positions(rows, min_identity):
    """
    {CT_ID: {cycle_id: set of covered subject positions}}.

    Positions rather than alignment lengths, so overlapping HSPs are not
    double-counted -- a cycle containing two copies of the same IS produces
    exactly that, and summing lengths would inflate coverage past 100%.
    """
    coverage = collections.defaultdict(lambda: collections.defaultdict(set))
    for query, subject, pident, _, _, _, _, _, sstart, send in rows:
        if float(pident) < min_identity:
            continue
        low, high = sorted((int(sstart), int(send)))
        coverage[ct_id_from_subject(subject)][query].update(range(low, high + 1))
    return coverage


def query_covered_fraction(rows, min_identity):
    """{cycle_id: fraction of the CYCLE explained by ground-truth CTs}."""
    positions = collections.defaultdict(set)
    lengths = {}
    for query, _, pident, _, qlen, _, qstart, qend, _, _ in rows:
        if float(pident) < min_identity:
            continue
        lengths[query] = int(qlen)
        low, high = sorted((int(qstart), int(qend)))
        positions[query].update(range(low, high + 1))
    return {query: len(pos) / lengths[query] for query, pos in positions.items()
            if lengths.get(query)}


def read_cycle_ids(path):
    return [line[1:].strip() for line in open(path) if line.startswith(">")]


def score(ground_truth, rows, cycle_ids, min_identity, min_coverage):
    coverage = covered_positions(rows, min_identity)

    matches = {}
    for ct_id, record in ground_truth.items():
        length = int(record["CT_Length"])
        best_cycle, best_fraction = None, 0.0
        for cycle, positions in coverage.get(ct_id, {}).items():
            fraction = len(positions) / length
            if fraction > best_fraction:
                best_fraction, best_cycle = fraction, cycle
        matches[ct_id] = (best_cycle if best_fraction >= min_coverage else None,
                          best_fraction)

    recovered = [ct for ct, (cycle, _) in matches.items() if cycle]

    # Precision: a reported cycle counts as a true positive only when it
    # represents a COMPLETE element -- it must cover at least `min_coverage` of
    # some ground-truth composite transposon.
    #
    # An earlier version asked only whether the candidate's own sequence was
    # CT-derived, which is a weaker and misleading test: a fragment of a
    # composite transposon is made entirely of that element's sequence and
    # passed, so partial candidates covering 78% of an element were counted as
    # true positives and precision came out too high. Being built from the right
    # sequence is not the same as being the right thing.
    covers_an_element = set()
    for ct_id, per_cycle in coverage.items():
        length = int(ground_truth[ct_id]["CT_Length"])
        for cycle, positions in per_cycle.items():
            if len(positions) / length >= min_coverage:
                covers_an_element.add(cycle)
    true_positives = [c for c in cycle_ids if c in covers_an_element]

    # Copy-distinctness: among CTs sharing an IS, how many came back as
    # SEPARATE cycles. Two ground-truth CTs mapped onto one cycle count once.
    by_element = collections.defaultdict(list)
    for ct_id, record in ground_truth.items():
        by_element[record["IS_Name"]].append(ct_id)
    shared_groups = {name: cts for name, cts in by_element.items() if len(cts) > 1}

    distinct_expected = sum(len(cts) for cts in shared_groups.values())
    distinct_found = 0
    for cts in shared_groups.values():
        distinct_found += len({matches[ct][0] for ct in cts if matches[ct][0]})

    # Recall split by whether the cargo was findable by homology at all. PICOTA's
    # score is largely a sum of database hits, so a headline recall that mixes
    # CARD cargo with novel cargo can hide a systematic blind spot.
    by_cargo = collections.defaultdict(lambda: [0, 0])
    for ct_id, record in ground_truth.items():
        cargo_type = record.get("Cargo_Type", "AMR")
        by_cargo[cargo_type][1] += 1
        if matches[ct_id][0]:
            by_cargo[cargo_type][0] += 1

    return {
        "matches": matches,
        "by_cargo": dict(by_cargo),
        "ct_recall": (len(recovered), len(ground_truth)),
        "precision": (len(true_positives), len(cycle_ids)),
        "copy_distinctness": (distinct_found, distinct_expected),
        "shared_groups": shared_groups,
    }


def report(result, ground_truth, handle=sys.stdout):
    print(f"{'CT':<8}{'IS':<14}{'len':>7}{'covered':>10}  status", file=handle)
    for ct_id in sorted(ground_truth):
        cycle, fraction = result["matches"][ct_id]
        status = f"RECOVERED  {cycle}" if cycle else "MISSED"
        print(f"{ct_id:<8}{ground_truth[ct_id]['IS_Name']:<14}"
              f"{ground_truth[ct_id]['CT_Length']:>7}{fraction * 100:>9.1f}%  {status}",
              file=handle)

    got, total = result["ct_recall"]
    print(f"\nCT recall:           {got}/{total}", file=handle)
    got, total = result["precision"]
    print(f"Precision:           {got}/{total} reported cycles explained by a CT",
          file=handle)
    for cargo_type, (got, total) in sorted(result.get("by_cargo", {}).items()):
        print(f"  {cargo_type + ' cargo:':<20} {got}/{total} recovered", file=handle)

    got, total = result["copy_distinctness"]
    if total:
        print(f"Copy-distinctness:   {got}/{total} CTs sharing an IS reported "
              f"as separate cycles", file=handle)
        for name, cts in sorted(result["shared_groups"].items()):
            print(f"  {name}: {', '.join(sorted(cts))}", file=handle)
    else:
        print("Copy-distinctness:   n/a (no IS shared between CTs)", file=handle)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Score PICOTA cycles against a simulated ground truth.")
    parser.add_argument("--ground-truth", required=True)
    parser.add_argument("--ground-truth-fasta", required=True)
    parser.add_argument("--cycles", required=True)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--min-coverage", type=float, default=0.95)
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    args = parser.parse_args(argv)

    ground_truth = read_ground_truth(args.ground_truth)
    cycle_ids = read_cycle_ids(args.cycles)
    rows = run_blast(args.cycles, args.ground_truth_fasta,
                     args.blastn, args.makeblastdb)
    report(score(ground_truth, rows, cycle_ids,
                 args.min_identity, args.min_coverage), ground_truth)
    return 0


if __name__ == "__main__":
    sys.exit(main())
