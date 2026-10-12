#!/usr/bin/env python3
"""
build_is_db.py
--------------
Builds the insertion-sequence reference set PICOTA searches, from the merged
TnCentral + ISfinder FASTA and the ISfinder nucleotide release (IS.fna).

Two things are wrong with the merged FASTA as distributed
(_tncentral_nointegrall_isfinder-TNs.fasta), and this script fixes both:

1. Glued records. In 21 places the newline before a header, and the one after
   it, were lost before the file was wrapped at 60 columns, so a header sits
   in the middle of the previous record's sequence:

       ...gtgcagccgtcttctgaaaacgaca>In1223-KX784502tgtcattt
       ...gtaaactctgatttatttgtggaactgaca>T
       n6023-GU562437.2ggcactgttgcaaagttagcgatgaggcagcc

   BLAST reads each of these as one longer sequence, so the second record is
   never searched under its own name and the first carries foreign sequence.
   The header is recovered as everything before the trailing run of lowercase
   a/c/g/t/n: accessions end in a digit or an uppercase letter, so the boundary
   is unambiguous. A record that does not fit that shape stops the build.

2. Missing ISfinder elements. Elements present in IS.fna but absent from the
   merged set, by name and by sequence (either strand), are appended as
   ">ISFinder_<name>", the convention the merged set already uses, so
   output_formatter.infer_is_family() resolves their family from
   isfinder_taxonomy.tsv.

   MITEs are not added (--include-mites adds them). They are short,
   non-autonomous elements that encode no transposase, so they cannot be the
   flanking insertion sequences of a composite transposon. Counted as IS hits,
   they let a cycle through score3's IS gate on the strength of a MITE: on the
   simulated benchmark one non-transposon cycle rose from 17 to 29 that way.
   The merged set carries 2 of ISfinder's 60 and keeps them.

Usage (from the repository root):
  python scripts/build_is_db.py
  python scripts/build_is_db.py --merged A.fasta --isfinder IS.fna --out B.fasta

The output is deterministic: the merged records in their original order,
then the added ISfinder records in IS.fna order. A summary goes to stderr.
"""

import argparse
import os
import re
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_IS_DIR     = os.path.join(_SCRIPT_DIR, '..', 'picota', 'DBs', 'ISes')
MERGED_PATH = os.path.join(_IS_DIR, '_tncentral_nointegrall_isfinder-TNs.fasta')
IS_FNA_PATH = os.path.join(_IS_DIR, 'IS.fna')
OUT_PATH    = os.path.join(_IS_DIR, 'picota_IS_db.fasta')

ISFINDER_PREFIX = 'ISFinder_'
WRAP = 60

_TRAILING_SEQ = re.compile(r'[acgtn]+$')
# Name-ACCESSION: the accession, after the last '-', is uppercase, digits,
# '_' and '.', and ends in a digit or capital.
_ACCESSION    = re.compile(r'-[A-Z0-9_.]*[A-Z0-9]$')
_COMPLEMENT   = str.maketrans('ACGTN', 'TGCAN')


class BuildError(Exception):
    pass


def read_fasta_lines(path):
    with open(path, encoding='utf-8', errors='replace') as fh:
        return [line.rstrip('\r\n') for line in fh]


def parse_merged(lines):
    """Records of the merged set as (header, sequence), glued records repaired.

    Returns (records, repaired_headers).
    """
    records, repaired = [], []
    header, seq = None, []
    # Text that began mid-line with '>' and whose header has not been
    # separated from its sequence yet.
    glued = None

    def close():
        nonlocal header, seq, glued
        if glued is not None:
            header_text, seq_text = split_glued(''.join(glued))
            repaired.append(header_text)
            records.append((header_text, seq_text))
            glued = None
        elif header is not None:
            records.append((header, ''.join(seq)))
        header, seq = None, []

    for line in lines:
        if line.startswith('>'):
            close()
            header = line[1:].strip()
            continue
        if '>' in line:
            before, after = line.split('>', 1)
            if '>' in after:
                raise BuildError(f"two headers on one line: {line[:80]}")
            if glued is not None:
                glued.append(before)
            else:
                seq.append(before)
            close()
            glued = [after]
            continue
        if glued is not None:
            glued.append(line.strip())
        else:
            seq.append(line.strip())
    close()
    return records, repaired


def split_glued(text):
    """'Name-ACC' + lowercase sequence, run together, back into its two parts."""
    m = _TRAILING_SEQ.search(text)
    if not m or m.start() == 0:
        raise BuildError(f"cannot find the sequence in glued record: {text[:80]}")
    header, seq = text[:m.start()], text[m.start():]
    if not _ACCESSION.search(header):
        raise BuildError(f"recovered header is not Name-ACCESSION: {header!r}")
    return header, seq


def parse_is_fna(lines):
    """IS.fna as (name, sequence); header is >Name[_variant]_Family_Group."""
    records, name, seq = [], None, []
    for line in lines:
        if line.startswith('>'):
            if name is not None:
                records.append((name, ''.join(seq)))
            parts = line[1:].strip().split('_')
            if len(parts) < 3:
                raise BuildError(f"IS.fna header without family and group: {line}")
            name, seq = '_'.join(parts[:-2]), []
        else:
            seq.append(line.strip())
    if name is not None:
        records.append((name, ''.join(seq)))
    return records


def canonical(seq):
    """One key for a sequence and its reverse complement."""
    s = seq.upper()
    rc = s.translate(_COMPLEMENT)[::-1]
    return min(s, rc)


def merged_name(header):
    """The element name a merged-set header refers to."""
    if header.startswith(ISFINDER_PREFIX):
        return header[len(ISFINDER_PREFIX):]
    return header.rsplit('-', 1)[0] if '-' in header else header


def is_mite(name):
    return name.upper().startswith('MITE')


def build(merged_lines, is_fna_lines, include_mites=False):
    """Return (records, report) for the combined set."""
    merged, repaired = parse_merged(merged_lines)
    names = {merged_name(h).lower() for h, _ in merged}
    seqs = {canonical(s) for _, s in merged if s}

    added, skipped_mites = [], []
    for name, seq in parse_is_fna(is_fna_lines):
        if not seq or name.lower() in names or canonical(seq) in seqs:
            continue
        if is_mite(name) and not include_mites:
            skipped_mites.append(name)
            continue
        added.append((ISFINDER_PREFIX + name, seq))
        names.add(name.lower())
        seqs.add(canonical(seq))

    empty = [h for h, s in merged if not s]
    if empty:
        raise BuildError(f"{len(empty)} records without sequence, e.g. {empty[0]}")
    report = {
        'merged_records': len(merged),
        'repaired': repaired,
        'added': [h for h, _ in added],
        'skipped_mites': skipped_mites,
        'total': len(merged) + len(added),
    }
    return merged + added, report


def write_fasta(records, path):
    with open(path, 'w') as fh:
        for header, seq in records:
            fh.write(f'>{header}\n')
            for i in range(0, len(seq), WRAP):
                fh.write(seq[i:i + WRAP] + '\n')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--merged', default=MERGED_PATH, help='Merged TnCentral + ISfinder FASTA')
    ap.add_argument('--isfinder', default=IS_FNA_PATH, help='ISfinder IS.fna')
    ap.add_argument('--out', default=OUT_PATH, help='Output FASTA. Default: %(default)s')
    ap.add_argument('--include-mites', action='store_true',
                    help='Also add the MITEs missing from the merged set')
    args = ap.parse_args(argv)

    try:
        records, report = build(read_fasta_lines(args.merged), read_fasta_lines(args.isfinder),
                                include_mites=args.include_mites)
    except BuildError as e:
        sys.exit(f"build_is_db: {e}")
    write_fasta(records, args.out)

    print(f"merged set:          {report['merged_records']} records "
          f"({len(report['repaired'])} glued records repaired)", file=sys.stderr)
    print(f"ISfinder added:      {len(report['added'])} "
          f"({len(report['skipped_mites'])} MITEs left out)", file=sys.stderr)
    print(f"written:             {report['total']} records -> {args.out}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
