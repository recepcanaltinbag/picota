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


# ─── SRA + Assembly helpers (full mode) ──────────────────────────────────────

def step_sra_download(short_id: str, raw_dir: str, sra_dir: str,
                      fastq_dump_path: str, logger) -> list:
    """Download SRA FASTQ files if not already present. Returns list of FASTQ paths."""
    from src.sra_download import run_sra_down

    expected = [os.path.join(raw_dir, f"{short_id}_{i}.fastq") for i in (1, 2)]
    missing  = [f for f in expected if not os.path.exists(f)]
    if not missing:
        logger.info(f"  FASTQ already present — skipping download")
        return [f for f in expected if os.path.exists(f)]

    logger.info(f"  Downloading {short_id} from SRA ...")
    run_sra_down(short_id, raw_dir, sra_dir, fastq_dump_path,
                 keep_sra_file=True, the_force=False, logger_name='picota_complete')
    return [f for f in expected if os.path.exists(f)]


def step_assembly(short_id: str, raw_files: list, asm_dir: str,
                  threads: int, logger) -> list:
    """Run MEGAHIT assembly if no GFA exists yet. Returns list of GFA paths."""
    from src.assembly import assembly_main

    gfa_files = [f for f in (Path(asm_dir).glob('*.gfa'))]
    if gfa_files:
        logger.info(f"  Assembly already exists — skipping")
        return [str(f) for f in gfa_files]

    logger.info(f"  Running assembly for {short_id} ...")
    assembly_main(
        short_id, raw_files, asm_dir,
        threads, "99",
        True,   # assembly_quiet
        False,  # assembly_keep_temp_files
        "spades.py", "fastp", False,
        "megahit", "megahit",
        "", "", 'picota_complete'
    )
    return [str(f) for f in Path(asm_dir).glob('*.gfa')]


# ─── Main pipeline ───────────────────────────────────────────────────────────

def _process_sample(short_id, long_id, output_path, gfa_mode,
                    missing_tools, long_read_threads, logger):
    """Run the full pipeline for a single sample. Returns list of enriched rows."""
    STEPS = 5
    logger.info("\n" + "═" * BANNER_WIDTH)
    logger.info(f"  Sample {short_id}  |  long-read: {long_id or 'none'}")
    logger.info("═" * BANNER_WIDTH)

    sample_dir = output_path / short_id
    sample_dir.mkdir(exist_ok=True)

    # ── Step 1/5: SRA download + Assembly ────────────────────────────────────
    logger.info(f"\n  [1/{STEPS}] SRA Download + Assembly")
    t0 = time.time()
    if gfa_mode:
        gfa_file = TEST_GFA
        logger.info(f"  ⤼ GFA mode — using bundled testNitro.gfa")
    else:
        raw_dir = str(output_path / 'raw' / short_id)
        sra_dir = str(output_path / 'sra'  / short_id)
        asm_dir = str(output_path / 'assembly' / short_id)
        os.makedirs(raw_dir, exist_ok=True)
        os.makedirs(sra_dir, exist_ok=True)
        os.makedirs(asm_dir, exist_ok=True)

        raw_files = step_sra_download(short_id, raw_dir, sra_dir,
                                      'parallel-fastq-dump', logger)
        if not raw_files:
            logger.error(f"  ✗ No FASTQ files found after download — aborting {short_id}")
            return []

        gfa_list = step_assembly(short_id, raw_files, asm_dir,
                                 threads=4, logger=logger)
        if not gfa_list:
            logger.error(f"  ✗ Assembly produced no GFA — aborting {short_id}")
            return []
        gfa_file = gfa_list[0]

    logger.info(f"  ✓ GFA ready: {gfa_file}  ({time.time()-t0:.1f}s)")

    # ── Step 2/5: Cycle detection ─────────────────────────────────────────────
    logger.info(f"\n  [2/{STEPS}] Cycle Detection")
    t0 = time.time()
    cycle_fasta = str(sample_dir / f'{short_id}_cycles.fasta')
    n_cycles = step_cycle_detection(gfa_file, cycle_fasta, logger)
    if n_cycles == 0:
        logger.warning(f"  ⚠  No cycles detected — skipping scoring")
        return []
    logger.info(f"  ✓ {n_cycles} candidate cycles  ({time.time()-t0:.1f}s)")

    # ── Step 3/5: BLAST scoring ───────────────────────────────────────────────
    logger.info(f"\n  [3/{STEPS}] BLAST Scoring")
    score_dir = str(sample_dir / 'scoring')
    final_tab = os.path.join(score_dir, 'picota_final_tab')
    enriched_csv = os.path.join(score_dir, 'picota_enriched.csv')

    if os.path.exists(final_tab) and os.path.exists(enriched_csv):
        logger.info(f"  ⤼ Scoring already done — loading existing results")
        enriched = step_load_enriched(score_dir, logger)
    elif missing_tools:
        logger.warning(f"  ⚠  Scoring skipped — missing: {', '.join(missing_tools)}")
        return []
    else:
        t0 = time.time()
        os.makedirs(score_dir, exist_ok=True)
        final_tab = step_scoring(cycle_fasta, score_dir, logger)
        enriched  = step_load_enriched(score_dir, logger)
        logger.info(f"  ✓ {len(enriched)} enriched rows  ({time.time()-t0:.1f}s)")

    for row in enriched:
        row['SRA_ID'] = short_id

    # ── Step 4/5: Annotation split ────────────────────────────────────────────
    logger.info(f"\n  [4/{STEPS}] Annotation Split")
    annot_dir = str(sample_dir / 'annot')
    existing_annot = list(Path(annot_dir).glob('*.fasta')) if os.path.isdir(annot_dir) else []
    if existing_annot:
        logger.info(f"  ⤼ Annotation already exists — {len(existing_annot)} FASTA(s)")
        annotated = [str(f) for f in existing_annot]
    else:
        t0 = time.time()
        from src.split_cycle_coords_for_is import split_cycles_from_picota
        annotated = split_cycles_from_picota(final_tab, cycle_fasta, annot_dir, 0)
        logger.info(f"  ✓ {len(annotated)} annotated FASTA(s)  ({time.time()-t0:.1f}s)")

    # ── Step 5/5: Long-read validation ────────────────────────────────────────
    logger.info(f"\n  [5/{STEPS}] Long-read Validation")
    if long_id and tool_available('minimap2') and tool_available('samtools'):
        long_fastq = output_path / 'raw_long' / long_id / f'{long_id}_1.fastq'
        if long_fastq.exists():
            t0 = time.time()
            map_dir = str(output_path / 'mapping' / long_id)
            n_mapped = 0
            for fa in annotated:
                try:
                    step_long_read_mapping(fa, str(long_fastq), map_dir,
                                           long_id, long_read_threads, logger)
                    n_mapped += 1
                except Exception as e:
                    logger.error(f"  Mapping error: {e}")
            logger.info(f"  ✓ {n_mapped} FASTA(s) mapped  ({time.time()-t0:.1f}s)")
        else:
            logger.info(f"  ⤼ Long-read FASTQ not found: {long_fastq}")
    else:
        logger.info(f"  ⤼ Long-read validation skipped (no long-read ID or tools missing)")

    return enriched


def run_pipeline(sra_list_file: str, output_dir: str, gfa_mode: bool = False,
                 long_read_threads: int = 4):
    """Execute the complete PICOTA pipeline for all samples in sra_list_file."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    log_file = output_path / 'picota_complete.log'
    logger   = PICOTALogger.setup(str(log_file), level='INFO')
    progress = AnalysisProgress(logger)

    progress.add_step('setup',   'Environment check')
    progress.add_step('samples', 'Per-sample analysis')
    progress.add_step('export',  'Results export')

    logger.info("=" * BANNER_WIDTH)
    logger.info("  PICOTA — Complete Analysis Pipeline")
    logger.info(f"  Mode   : {'GFA test (bundled testNitro.gfa)' if gfa_mode else 'Full (SRA download + assembly)'}")
    logger.info(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"  Output : {output_path.absolute()}")
    logger.info("=" * BANNER_WIDTH)

    # ── Environment check ─────────────────────────────────────────────────────
    progress.start_step('setup')
    missing = [t for t in ('prodigal', 'blastn', 'blastp', 'makeblastdb')
               if not tool_available(t)]
    if missing:
        logger.warning(f"  Missing tools: {', '.join(missing)}")
        logger.warning("  Install via: conda install -c bioconda " + " ".join(missing))
    else:
        logger.info("  All required tools available ✓")
    progress.complete_step('setup')

    # ── Load SRA list ─────────────────────────────────────────────────────────
    if not os.path.exists(sra_list_file):
        logger.error(f"SRA list not found: {sra_list_file}")
        sys.exit(1)

    with open(sra_list_file, newline='') as fh:
        samples = []
        for row in csv.DictReader(fh):
            short_id = row.get('sra_short_id', '').strip()
            long_id  = row.get('sra_long_id', '').strip()
            if short_id:
                long_id = None if long_id in ('', '-', 'null', 'None') else long_id
                samples.append((short_id, long_id))

    logger.info(f"  Loaded {len(samples)} sample(s)\n")

    # ── Process samples ───────────────────────────────────────────────────────
    progress.start_step('samples')
    all_enriched = []
    t_pipeline   = time.time()

    for i, (short_id, long_id) in enumerate(samples, 1):
        logger.info(f"\n{'█'*BANNER_WIDTH}")
        logger.info(f"  [{i}/{len(samples)}]  {short_id}")
        logger.info(f"{'█'*BANNER_WIDTH}")
        try:
            rows = _process_sample(short_id, long_id, output_path,
                                   gfa_mode, missing, long_read_threads, logger)
            all_enriched.extend(rows)
            logger.info(f"\n  ✅ {short_id} done — {len(rows)} enriched rows")
        except Exception as e:
            logger.error(f"\n  ✗ {short_id} failed: {e}")
            import traceback
            logger.debug(traceback.format_exc())

    progress.complete_step('samples', f"{len(samples)} sample(s), {len(all_enriched)} total CT rows")

    # ── Export combined results ───────────────────────────────────────────────
    progress.start_step('export')
    if all_enriched:
        csv_out  = output_path / 'picota_enriched_combined.csv'
        json_out = output_path / 'picota_enriched_combined.json'
        ResultsFormatter.to_csv(all_enriched,  str(csv_out))
        ResultsFormatter.to_json(all_enriched, str(json_out))
        logger.info(f"\n  CSV  → {csv_out}")
        logger.info(f"  JSON → {json_out}")
        progress.complete_step('export', f"{len(all_enriched)} rows written")

        ResultsFormatter.print_summary(
            [{'Category': r.get('Category',''), 'Score': r.get('Score','0'),
              'CT_ID': r.get('CT_Tag',''), 'Antibiotic_Classes': r.get('Antibiotic_Class','')}
             for r in all_enriched], logger
        )
    else:
        progress.complete_step('export', "no results to export")

    progress.report_summary()
    logger.info(f"\nTotal time: {time.time()-t_pipeline:.1f}s")
    logger.info("=" * BANNER_WIDTH)
    logger.info("  PICOTA pipeline finished")
    logger.info("=" * BANNER_WIDTH)


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    _here = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description="PICOTA — Complete Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick test with bundled testNitro.gfa (no internet, no assembly):
  conda run -n evobiomig python3 test_complete_pipeline.py --gfa_mode

  # Full run from real SRA accessions:
  conda run -n evobiomig python3 test_complete_pipeline.py \\
      --sra_list picota/test_sra_ids.csv --output /data/picota_results/
        """
    )
    parser.add_argument('--sra_list', '-s',
                        default=os.path.join(_here, 'picota', 'test_sra_ids.csv'),
                        help='CSV with sra_short_id,sra_long_id columns')
    parser.add_argument('--output', '-o', default='picota_results',
                        help='Output directory (default: picota_results)')
    parser.add_argument('--gfa_mode', action='store_true',
                        help='Quick test: use bundled testNitro.gfa, skip SRA/assembly')
    parser.add_argument('--threads', '-t', type=int, default=4,
                        help='Assembly / mapping threads (default: 4)')

    args = parser.parse_args()

    # --gfa_mode: always use a single testNitro entry (ignores real SRA IDs)
    if args.gfa_mode:
        gfa_csv = os.path.join(_here, 'picota', 'testNitro_sra_ids.csv')
        with open(gfa_csv, 'w') as fh:
            fh.write('sra_short_id,sra_long_id\ntestNitro,-\n')
        args.sra_list = gfa_csv

    run_pipeline(args.sra_list, args.output,
                 gfa_mode=args.gfa_mode,
                 long_read_threads=args.threads)


if __name__ == '__main__':
    main()
