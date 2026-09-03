#!/usr/bin/env bash
#
# build_report.sh
# ---------------
# Every published number, from the benchmark directories, in one command.
#
# The point is not convenience. Each table in the paper was produced by a
# different script with its own flags, and a flag that drifts between a run and
# a re-run is exactly how the compact scenario came to report the wrong number
# for weeks -- the harness passed min_size_of_cycle=2000 while the shipped
# config said 1000, and nothing downstream could tell. Here the flags are
# written down once and every table comes from the same invocation.
#
# Usage:
#   scripts/build_report.sh <benchmark-dir> <output-dir>
#
# Expects <benchmark-dir> to contain scenarios_v2/ and copysweep/ as produced by
# run_scenarios.py and run_copy_sweep.py, and ref/ with the host genomes.

set -euo pipefail

BENCH="${1:?usage: build_report.sh <benchmark-dir> <output-dir>}"
OUT="${2:?usage: build_report.sh <benchmark-dir> <output-dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

THRESHOLD="${THRESHOLD:-50}"
SCORE_TYPE="${SCORE_TYPE:-3}"
MIN_COVERAGE="${MIN_COVERAGE:-0.95}"

mkdir -p "$OUT"

echo "==> scenario table (score${SCORE_TYPE}, threshold ${THRESHOLD})"
python3 "$HERE/score_scenarios.py" \
    --scenarios "$BENCH/scenarios_v2" \
    --threshold "$THRESHOLD" --score-type "$SCORE_TYPE" \
    | tee "$OUT/table_scenarios.txt"

echo
echo "==> copy-number sweep (detection only, ${MIN_COVERAGE} coverage)"
python3 "$HERE/score_copy_sweep.py" \
    --sweep "$BENCH/copysweep" --min-coverage "$MIN_COVERAGE" \
    | tee "$OUT/table_copy_sweep.txt"

echo
echo "==> contig baseline on the same assemblies"
python3 "$HERE/compare_contig_recovery.py" \
    "$BENCH/scenarios_v2" --assembler spades megahit \
    --min-coverage "$MIN_COVERAGE" \
    | tee "$OUT/table_contig_baseline.txt"

echo
echo "==> element catalogue"
python3 "$HERE/catalogue_elements.py" \
    --runs "$BENCH/scenarios_v2" "$BENCH/copysweep" \
    --out "$OUT/elements.tsv" \
    | tee "$OUT/catalogue_summary.txt"

echo
echo "==> supplementary package"
python3 "$HERE/package_supplementary.py" \
    --runs "$BENCH/scenarios_v2" "$BENCH/copysweep" \
    --out "$OUT/supplementary"
cp "$OUT/elements.tsv" "$OUT/supplementary/elements.tsv"

echo
echo "Report written to $OUT"
ls -1 "$OUT"
