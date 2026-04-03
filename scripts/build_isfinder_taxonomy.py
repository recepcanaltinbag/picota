#!/usr/bin/env python3
"""
build_isfinder_taxonomy.py
--------------------------
Parses the bundled IS.fna (ISFinder nucleotide sequences) to extract
IS element -> superfamily/group mapping and writes a TSV lookup table.

Header format in IS.fna:
  >ISname_ISfamily_ISgroup          (3 parts, the common case)
  >ISname_variant_ISfamily_ISgroup  (4 parts, e.g. ISBj2_B_IS5_IS5)

Output: DBs/ISes/isfinder_taxonomy.tsv
  Columns: IS_Name  IS_Family  IS_Group

Usage (run from the picota/ directory):
  python scripts/build_isfinder_taxonomy.py

This script needs to be re-run only if IS.fna is updated.
The generated TSV is used by output_formatter.infer_is_family() as the
primary lookup source (~5900 IS elements), falling back to the curated
_IS_EXACT dict for elements not in IS.fna.
"""

import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PICOTA_DIR = os.path.join(_SCRIPT_DIR, '..', 'picota')
FASTA_PATH  = os.path.join(_PICOTA_DIR, 'DBs', 'ISes', 'IS.fna')
OUT_PATH    = os.path.join(_PICOTA_DIR, 'DBs', 'ISes', 'isfinder_taxonomy.tsv')


def parse_is_fna(fasta_path: str) -> list:
    """
    Parse IS.fna headers and return list of {name, family, group} dicts.

    Header format: >IS_Name[_variant]_IS_Family_IS_Group
    The last two underscore-fields are always family and group.
    """
    records = []
    seen = set()

    with open(fasta_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            header = line[1:].strip()
            parts = header.split('_')
            if len(parts) < 3:
                continue

            # Last two parts are always family / group
            group  = parts[-1]
            family = parts[-2]
            # Name = everything before the last two fields
            name   = '_'.join(parts[:-2])

            # Normalise "unknown" -> use family value for both
            if group.lower()  == 'unknown':
                group = family
            if family.lower() == 'unknown':
                family = group

            if not name or not family:
                continue

            key = name.lower()
            if key not in seen:
                seen.add(key)
                records.append({'IS_Name': name, 'IS_Family': family, 'IS_Group': group})

    return records


def write_tsv(records: list, out_path: str):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write('IS_Name\tIS_Family\tIS_Group\n')
        for r in records:
            fh.write(f"{r['IS_Name']}\t{r['IS_Family']}\t{r['IS_Group']}\n")


def main():
    if not os.path.exists(FASTA_PATH):
        print(f"[ERROR] IS.fna not found at: {FASTA_PATH}", file=sys.stderr)
        sys.exit(1)

    records = parse_is_fna(FASTA_PATH)
    write_tsv(records, OUT_PATH)
    print(f"[OK] {len(records)} IS elements written to {OUT_PATH}")

    # Sanity check on known problematic cases (the classic prefix-collision traps)
    lookup = {r['IS_Name'].lower(): r for r in records}
    tests = [
        ('IS26',   'IS6',          'IS6'),
        ('IS10',   'IS4',          'IS4'),
        ('IS21',   'IS21',         'IS21'),
        ('IS2',    'IS3',          'IS3'),
        ('IS1380', 'IS1380',       'IS1380'),
        ('ISEcp1', 'IS1380',       'IS1380'),
        ('ISCR1',  'IS91',         'IS91'),
    ]
    print("\nSanity checks (IS_Name → IS_Family):")
    all_ok = True
    for name, exp_fam, exp_grp in tests:
        r = lookup.get(name.lower())
        if r:
            ok = r['IS_Family'] == exp_fam
            mark = 'OK' if ok else 'FAIL'
            if not ok:
                all_ok = False
            print(f"  [{mark}] {name:<12} family={r['IS_Family']:<14} group={r['IS_Group']}"
                  + ('' if ok else f"  <-- expected {exp_fam}"))
        else:
            print(f"  [?]  {name:<12} not found in IS.fna")

    if all_ok:
        print("\nAll sanity checks passed.")
    else:
        print("\nSome checks failed — review IS.fna header format.")
        sys.exit(1)


if __name__ == '__main__':
    main()
