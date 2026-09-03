#!/usr/bin/env python3
"""
run_manifest.py
---------------
Record what actually produced a benchmark case.

ground_truth.json records how the genome was built. It says nothing about how
the genome was then sequenced, assembled and traversed, and those choices decide
the result as much as the genome does: the compact scenario reported 8 of 12
elements for weeks because the harness passed min_size_of_cycle=2000 while the
shipped config said 1000, and nothing in the case directory recorded which had
been used.

So each case also gets run_parameters.json:

  parameters   the harness settings that reach the pipeline -- cycle size,
               dedup mode, coverage, k-mers, which assemblers ran
  tools        the version string each external program reported, captured at
               run time rather than assumed. megahit_core is omitted: it prints
               usage for every version flag and its build travels with megahit,
               whose version is recorded.
  command      the exact argv

Usable as a CLI to backfill a directory of cases whose run predates this file:

  python scripts/run_manifest.py --cases sweep/ --min-cycle-size 1000 \\
      --coverage 50 --assembler spades megahit
"""

import argparse
import json
import os
import subprocess
import sys

VERSION_COMMANDS = {
    "spades": (["spades.py", "--version"], None),
    "megahit": (["megahit", "--version"], None),
    # ART's banner names the program on one line and the version on the next,
    # so match the version line: "ART_Illumina" alone records a date range.
    "art_illumina": (["art_illumina"], "Version"),
    "prodigal": (["prodigal", "-v"], None),
    "blastn": (["blastn", "-version"], None),
    "makeblastdb": (["makeblastdb", "-version"], None),
}


def tool_version(command, marker=None, executable=None):
    """
    First informative line the program prints, or None when it is absent.

    ART writes its banner to stderr and exits non-zero with no arguments, so
    output is captured from both streams and the return code ignored; a version
    string that only appears on failure is still the version string.
    """
    argv = list(command)
    if executable:
        argv[0] = executable
    try:
        result = subprocess.run(argv, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (result.stdout or "") + "\n" + (result.stderr or "")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if marker:
        for line in lines:
            if marker in line:
                return line
    return lines[0] if lines else None


def collect_tool_versions(overrides=None):
    overrides = overrides or {}
    versions = {}
    for name, (command, marker) in VERSION_COMMANDS.items():
        versions[name] = tool_version(command, marker, overrides.get(name))
    return versions


def write_run_parameters(case_dir, parameters, tool_overrides=None,
                         command=None, versions=None):
    payload = {
        "parameters": parameters,
        "tools": versions if versions is not None
                 else collect_tool_versions(tool_overrides),
        "command": command if command is not None else sys.argv,
    }
    path = os.path.join(case_dir, "run_parameters.json")
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("--cases", required=True,
                        help="Directory whose immediate subdirectories are cases.")
    parser.add_argument("--min-cycle-size", type=int, required=True)
    parser.add_argument("--coverage", type=float, required=True)
    parser.add_argument("--assembler", nargs="+", required=True)
    parser.add_argument("--kmers", default="55,77,99")
    parser.add_argument("--dedup-mode", default="strict")
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--art", default="art_illumina")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)

    parameters = {
        "min_cycle_size": args.min_cycle_size,
        "coverage": args.coverage,
        "assemblers": args.assembler,
        "kmers": args.kmers,
        "dedup_mode": args.dedup_mode,
        "read_length": args.read_length,
        "read_simulator": "art_illumina",
        "art_profile": "HSXt",
    }
    # Probe once rather than per case: the tools do not change mid-directory,
    # and each probe spawns a process.
    versions = collect_tool_versions({"art_illumina": args.art})

    written = 0
    for name in sorted(os.listdir(args.cases)):
        case = os.path.join(args.cases, name)
        if not os.path.isdir(case):
            continue
        if not os.path.exists(os.path.join(case, "ground_truth.json")):
            continue
        target = os.path.join(case, "run_parameters.json")
        if os.path.exists(target) and not args.overwrite:
            continue
        write_run_parameters(case, parameters, command=["backfilled"],
                             versions=versions)
        written += 1
        print("wrote %s" % target)

    print("\n%d case(s) recorded" % written)
    for name, version in sorted(versions.items()):
        print("  %-14s %s" % (name, version or "NOT FOUND"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
