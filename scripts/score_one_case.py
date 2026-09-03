#!/usr/bin/env python3
"""
score_one_case.py
-----------------
Score one benchmark case, so the cache can be filled in parallel.

score_scenarios.py walks its cases in sequence and each one runs Prodigal and
four BLAST passes. On a 24-core machine that left 18 cores idle and cost 2.7
minutes per case -- two hours for a 60-case run. The work is independent per
case and already cached on picota_final_tab, so warming that cache from separate
processes and letting score_scenarios.py read it back produces the same numbers
from the same code path, about three times faster.

This is a scheduling wrapper, not a second implementation: it calls
score_scenarios.score_cycles, which is what the aggregate would call anyway.

Usage:
  python scripts/score_one_case.py <case-dir> [--threshold 50] [--score-type 3]
  ls -d bench/scenarios_v2/*_s*/ | xargs -P 4 -I{} python scripts/score_one_case.py {}
"""

import argparse
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from score_scenarios import score_cycles  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    parser.add_argument("case")
    parser.add_argument("--threshold", type=int, default=50)
    parser.add_argument("--score-type", type=int, default=3, choices=(0, 1, 2, 3))
    parser.add_argument("--prodigal", default="prodigal")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--blastx", default="blastx")
    parser.add_argument("--blastp", default="blastp")
    args = parser.parse_args(argv)

    tools = {"path_of_prodigal": args.prodigal, "path_of_blastn": args.blastn,
             "path_of_makeblastdb": args.makeblastdb,
             "path_of_blastx": args.blastx, "path_of_blastp": args.blastp}

    scores = score_cycles(args.case.rstrip("/"), args.threshold, tools,
                          args.score_type)
    print("%s: %d above threshold" % (os.path.basename(args.case.rstrip("/")),
                                      len(scores)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
