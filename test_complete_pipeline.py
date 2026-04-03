#!/usr/bin/env python3
"""
test_complete_pipeline.py – Real end-to-end PICOTA pipeline runner.

Reads a CSV file with SRA accession pairs (short-read + optional long-read),
downloads data, assembles, detects cycles, scores, and exports enriched results.

Usage
-----
  # Using bundled test GFA (no SRA download, no assembly):
  conda run -n evobiomig python3 test_complete_pipeline.py --gfa_mode

  # Full run from an SRA list:
  conda run -n evobiomig python3 test_complete_pipeline.py \\
      --sra_list my_samples.csv --config picota/config.yaml --output results/

SRA list format (CSV, header required):
  sra_short_id,sra_long_id
  SRR12345678,SRR98765432
  SRR11111111,-

Use '-' in sra_long_id when no long-read data is available.
"""

import argparse
import csv
import os
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

# Add picota src to path
SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'picota', 'src')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'picota'))
sys.path.insert(0, SRC_DIR)

from src.logger_professional import PICOTALogger, AnalysisProgress, ResultsFormatter

# ─── Constants ───────────────────────────────────────────────────────────────
BANNER_WIDTH = 70
TEST_GFA     = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'picota', 'test_data', 'testNitro.gfa')
DB_BASE      = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'picota', 'DBs')


def db_path(*parts):
    return os.path.join(DB_BASE, *parts)


def tool_available(name: str) -> bool:
    return shutil.which(name) is not None


# ─── Real pipeline steps ─────────────────────────────────────────────────────

def step_cycle_detection(gfa_file: str, out_fasta: str, logger) -> int:
    """Run cycle detection on a GFA file. Returns number of cycles found."""
    from src.cycle_finderv2 import cycle_analysis

    if os.path.exists(out_fasta):
        n = sum(1 for l in open(out_fasta) if l.startswith('>'))
        logger.info(f"  Cycle FASTA already exists — {n} cycles (skipping)")
        return n

    cycle_analysis(
        gfa_file, out_fasta,
        find_all_path=False, path_limit=15,
        min_size_of_cycle=3000, max_size_of_cycle=100_000,
        name_prefix_cycle='Cycle',
        min_component_number=1, max_component_number=25,
        k_mer_sim=200, threshold_sim=99
    )
    n = sum(1 for l in open(out_fasta) if l.startswith('>')) if os.path.exists(out_fasta) else 0
    logger.info(f"  Detected {n} candidate cycles")
    return n


def step_scoring(cycle_fasta: str, out_dir: str, logger) -> str:
    """Run BLAST scoring. Returns path to picota_final_tab."""
    from src.scoringv4ProtBlast import scoring_main

    final_tab = os.path.join(out_dir, 'picota_final_tab')
    if os.path.exists(final_tab):
        n = sum(1 for _ in open(final_tab)) - 1
        logger.info(f"  Scoring already done — {n} hits (skipping)")
        return final_tab

    scoring_main(
        cycle_fasta, out_dir,
        db_path('Antibiotics', 'protein_fasta_protein_homolog_model.fasta'),
        db_path('Xenobiotics', 'Xenobiotics_classified.fasta'),
        db_path('ISes', '_tncentral_nointegrall_isfinder-TNs.fasta'),
        db_path('CompTns', 'Known_Tns.fasta'),
        out_dir,
        mean_of_CompTns=5850, std_of_CompTns=2586,
        total_score_type=0, threshold_final_score=0,
        max_z=20, dist_type=1,
        path_of_prodigal='prodigal',
        path_of_blastn='blastn',
        path_of_makeblastdb='makeblastdb',
        path_of_blastx='blastx',
        path_of_blastp='blastp',
        logger_name='picota_complete'
    )
    n = (sum(1 for _ in open(final_tab)) - 1) if os.path.exists(final_tab) else 0
    logger.info(f"  Scored {n} cycles above threshold")
    return final_tab


def step_load_enriched(out_dir: str, logger) -> list:
    """Load picota_enriched.csv rows (generated automatically by scoring_main)."""
    from src.output_formatter import ENRICHED_COLS

    enriched_csv = os.path.join(out_dir, 'picota_enriched.csv')
    if not os.path.exists(enriched_csv):
        logger.warning("  picota_enriched.csv not found — enriched output missing")
        return []

    rows = []
    with open(enriched_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(dict(row))
    logger.info(f"  Loaded {len(rows)} enriched rows from {enriched_csv}")
    return rows


def step_long_read_mapping(fasta_record: str, long_fastq: str,
                            map_dir: str, sra_id: str, threads: int, logger) -> str:
    """Run minimap2 mapping of long reads to a CT FASTA. Returns sorted BAM path."""
    import subprocess

    os.makedirs(map_dir, exist_ok=True)
    bam_name   = f"{sra_id}_{os.path.basename(fasta_record)}_mapping.bam"
    sorted_bam = os.path.join(map_dir, bam_name.replace('.bam', '_sorted.bam'))
    bai        = sorted_bam + '.bai'

    if os.path.exists(sorted_bam) and os.path.exists(bai):
        logger.info(f"  BAM already exists: {sorted_bam}")
        return sorted_bam

    sam_path = sorted_bam.replace('_sorted.bam', '.sam')
    with open(sam_path, 'w') as fh:
        subprocess.run(
            ['minimap2', '-ax', 'map-ont', '-t', str(threads), fasta_record, long_fastq],
            stdout=fh, check=True
        )
    bam_path = sorted_bam.replace('_sorted.bam', '.bam')
    subprocess.run(['samtools', 'view', '-bS', sam_path, '-o', bam_path], check=True)
    subprocess.run(['samtools', 'sort', bam_path, '-o', sorted_bam], check=True)
    subprocess.run(['samtools', 'index', sorted_bam], check=True)
    if os.path.exists(sam_path):
        os.remove(sam_path)
    logger.info(f"  Mapped: {sorted_bam}")
    return sorted_bam


# ─── Main pipeline ───────────────────────────────────────────────────────────

def run_pipeline(sra_list_file: str, output_dir: str, gfa_mode: bool = False,
                 long_read_threads: int = 4):
    """Execute the complete PICOTA pipeline."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    log_file = output_path / 'picota_complete.log'
    logger   = PICOTALogger.setup(str(log_file), level='INFO')
    progress = AnalysisProgress(logger)

    # Define steps
    progress.add_step('setup',      'Environment check')
    progress.add_step('cycles',     'Cycle detection')
    progress.add_step('scoring',    'BLAST scoring + enriched CSV')
    progress.add_step('long_reads', 'Long-read validation')
    progress.add_step('export',     'Results export')

    logger.info("=" * BANNER_WIDTH)
    logger.info("  PICOTA — Complete Analysis Pipeline")
    logger.info(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"  Output : {output_path.absolute()}")
    logger.info("=" * BANNER_WIDTH)

    # ── Step 1: Environment ───────────────────────────────────────────────────
    progress.start_step('setup')
    missing = [t for t in ('prodigal', 'blastn', 'blastp', 'makeblastdb')
               if not tool_available(t)]
    if missing:
        logger.warning(f"  Missing tools: {', '.join(missing)}")
        logger.warning("  Install via:  conda install -c bioconda " + " ".join(missing))
    else:
        logger.info("  All required tools available")
    progress.complete_step('setup')

    # ── Step 2: Load samples ──────────────────────────────────────────────────
    if not os.path.exists(sra_list_file):
        logger.error(f"SRA list not found: {sra_list_file}")
        sys.exit(1)

    with open(sra_list_file, newline='') as fh:
        reader = csv.DictReader(fh)
        samples = []
        for row in reader:
            short_id = row.get('sra_short_id', '').strip()
            long_id  = row.get('sra_long_id', '').strip()
            if short_id:
                long_id = None if long_id in ('', '-', 'null', 'None') else long_id
                samples.append((short_id, long_id))

    logger.info(f"  Loaded {len(samples)} sample(s) from {sra_list_file}")

    all_enriched = []
    t_pipeline   = time.time()

    for short_id, long_id in samples:
        logger.info("\n" + "─" * BANNER_WIDTH)
        logger.info(f"  Sample: {short_id}  |  long-read: {long_id or 'none'}")
        logger.info("─" * BANNER_WIDTH)

        sample_dir = output_path / short_id
        sample_dir.mkdir(exist_ok=True)

        # ── Cycle detection ───────────────────────────────────────────────────
        progress.start_step('cycles')
        if gfa_mode:
            gfa_file = TEST_GFA
            logger.info(f"  [GFA mode] Using bundled test GFA: {gfa_file}")
        else:
            # In full mode the GFA comes from the assembly step
            # (assembly is handled by picota_testv3.py; here we expect it ready)
            gfa_candidates = list((output_path / 'assembly' / short_id).glob('*.gfa'))
            if not gfa_candidates:
                logger.error(f"  No GFA found for {short_id} — run assembly first or use --gfa_mode")
                continue
            gfa_file = str(gfa_candidates[0])

        cycle_fasta = str(sample_dir / f'{short_id}_cycles.fasta')
        n_cycles = step_cycle_detection(gfa_file, cycle_fasta, logger)
        if n_cycles == 0:
            logger.warning(f"  No cycles detected for {short_id} — skipping scoring")
            progress.complete_step('cycles', f"0 cycles — skipped")
            continue
        progress.complete_step('cycles', f"{n_cycles} cycles")

        # ── BLAST scoring ─────────────────────────────────────────────────────
        if missing:
            logger.warning("  Scoring skipped — required tools missing")
            continue

        progress.start_step('scoring')
        score_dir = str(sample_dir / 'scoring')
        os.makedirs(score_dir, exist_ok=True)
        final_tab = step_scoring(cycle_fasta, score_dir, logger)
        enriched  = step_load_enriched(score_dir, logger)
        for row in enriched:
            row['SRA_ID'] = short_id
        all_enriched.extend(enriched)
        progress.complete_step('scoring', f"{len(enriched)} rows in enriched CSV")

        # ── Long-read validation ──────────────────────────────────────────────
        progress.start_step('long_reads')
        if long_id and tool_available('minimap2') and tool_available('samtools'):
            long_fastq = output_path / 'raw_long' / long_id / f'{long_id}_1.fastq'
            if long_fastq.exists():
                map_dir = str(output_path / 'mapping' / long_id)
                from src.split_cycle_coords_for_is import split_cycles_from_picota
                annotated = split_cycles_from_picota(final_tab, cycle_fasta,
                                                     str(sample_dir / 'annot'), 0)
                for fa in annotated:
                    try:
                        step_long_read_mapping(fa, str(long_fastq), map_dir,
                                               long_id, long_read_threads, logger)
                    except Exception as e:
                        logger.error(f"  Mapping error for {fa}: {e}")
                progress.complete_step('long_reads', f"{len(annotated)} FASTA(s) mapped")
            else:
                progress.complete_step('long_reads', f"FASTQ not found: {long_fastq}")
        else:
            progress.complete_step('long_reads', "skipped (no long-read ID or tools missing)")

    # ── Step 5: Export combined results ──────────────────────────────────────
    progress.start_step('export')
    if all_enriched:
        csv_out  = output_path / 'picota_enriched_combined.csv'
        json_out = output_path / 'picota_enriched_combined.json'
        ResultsFormatter.to_csv(all_enriched,  str(csv_out))
        ResultsFormatter.to_json(all_enriched, str(json_out))
        logger.info(f"  CSV  → {csv_out}")
        logger.info(f"  JSON → {json_out}")
        progress.complete_step('export', f"{len(all_enriched)} total rows written")

        # Print summary table
        ResultsFormatter.print_summary(
            [{
                'Category':          r.get('Category', ''),
                'Score':             r.get('Score', '0'),
                'CT_ID':             r.get('CT_Tag', ''),
                'Antibiotic_Classes': r.get('Antibiotic_Class', ''),
            } for r in all_enriched],
            logger
        )
    else:
        progress.complete_step('export', "no results to export")

    # ── Final summary ─────────────────────────────────────────────────────────
    progress.report_summary()
    elapsed = time.time() - t_pipeline
    logger.info(f"\nTotal pipeline time: {elapsed:.1f}s")
    logger.info("=" * BANNER_WIDTH)
    logger.info("  PICOTA pipeline completed successfully")
    logger.info("=" * BANNER_WIDTH)


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="PICOTA — Complete Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick test using bundled testNitro.gfa (no download/assembly needed):
  conda run -n evobiomig python3 test_complete_pipeline.py --gfa_mode

  # Full run from SRA list:
  conda run -n evobiomig python3 test_complete_pipeline.py \\
      --sra_list sra_ids.csv --output /data/picota_results/
        """
    )
    parser.add_argument('--sra_list', '-s',
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'picota', 'test_sra_ids.csv'),
                        help='CSV with sra_short_id,sra_long_id columns')
    parser.add_argument('--output', '-o', default='picota_results',
                        help='Output directory (default: picota_results)')
    parser.add_argument('--gfa_mode', action='store_true',
                        help='Use bundled testNitro.gfa instead of real assembly (for testing)')
    parser.add_argument('--threads', '-t', type=int, default=4,
                        help='Threads for minimap2 (default: 4)')

    args = parser.parse_args()

    # In GFA mode create a minimal SRA list if none exists
    if args.gfa_mode:
        default_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'picota', 'test_sra_ids.csv')
        if not os.path.exists(default_csv):
            with open(default_csv, 'w') as fh:
                fh.write('sra_short_id,sra_long_id\n')
                fh.write('testNitro,-\n')
        args.sra_list = default_csv

    run_pipeline(args.sra_list, args.output,
                 gfa_mode=args.gfa_mode,
                 long_read_threads=args.threads)


if __name__ == '__main__':
    main()
