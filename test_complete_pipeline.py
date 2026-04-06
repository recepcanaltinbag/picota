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

# Ensure picota/src is importable regardless of working directory.
# setup_picota_paths() may extend these further at runtime via --picota_dir.
_here = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(_here, 'picota', 'src'),
           os.path.join(_here, 'picota')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from src.logger_professional import PICOTALogger, AnalysisProgress, ResultsFormatter
from src.config_loader import load_config
from src.output_formatter import query_sra_organism

# ─── Constants ───────────────────────────────────────────────────────────────
BANNER_WIDTH = 70
TEST_GFA     = None  # Will be set by setup_picota_paths
DB_BASE      = None  # Will be set by setup_picota_paths


def setup_picota_paths(picota_dir: str):
    """Set up global paths based on picota_dir."""
    global TEST_GFA, DB_BASE, SRC_DIR

    SRC_DIR = os.path.join(picota_dir, 'picota', 'src')
    sys.path.insert(0, os.path.join(picota_dir, 'picota'))
    sys.path.insert(0, SRC_DIR)

    TEST_GFA = os.path.join(picota_dir, 'picota', 'test_data', 'testNitro.gfa')
    DB_BASE = os.path.join(picota_dir, 'picota', 'DBs')


def db_path(*parts):
    return os.path.join(DB_BASE, *parts)


def tool_available(name: str) -> bool:
    return shutil.which(name) is not None


# ─── Real pipeline steps ─────────────────────────────────────────────────────

def step_cycle_detection(gfa_file: str, out_fasta: str, logger, cfg=None) -> int:
    """Run cycle detection on a GFA file. Returns number of cycles found."""
    from src.cycle_finderv2 import cycle_analysis

    if os.path.exists(out_fasta):
        with open(out_fasta) as _f:
            n = sum(1 for l in _f if l.startswith('>'))
        logger.info(f"  Cycle FASTA already exists — {n} cycles (skipping)")
        return n

    # Read params from config if available, otherwise use safe defaults
    opt = getattr(cfg, 'options', None) if cfg else None
    pth = getattr(cfg, 'paths',   None) if cfg else None
    cycle_analysis(
        gfa_file, out_fasta,
        find_all_path     = getattr(pth, 'find_all_path',        False),
        path_limit        = getattr(pth, 'path_limit',           15),
        min_size_of_cycle = getattr(opt, 'min_size_of_cycle',    3000),
        max_size_of_cycle = getattr(opt, 'max_size_of_cycle',    100_000),
        name_prefix_cycle = getattr(opt, 'name_prefix_cycle',    'Cycle'),
        min_component_number = getattr(opt, 'min_component_number', 1),
        max_component_number = getattr(opt, 'max_component_number', 25),
        k_mer_sim         = getattr(opt, 'k_mer_sim',            200),
        threshold_sim     = getattr(opt, 'threshold_sim',        99),
    )
    if os.path.exists(out_fasta):
        with open(out_fasta) as _f:
            n = sum(1 for l in _f if l.startswith('>'))
    else:
        n = 0
    logger.info(f"  Detected {n} candidate cycles")
    return n


def step_scoring(cycle_fasta: str, out_dir: str, logger, cfg=None) -> str:
    """Run BLAST scoring. Returns path to picota_final_tab."""
    from src.scoringv4ProtBlast import scoring_main

    final_tab = os.path.join(out_dir, 'picota_final_tab')
    if os.path.exists(final_tab):
        with open(final_tab) as _f:
            n = sum(1 for _ in _f) - 1
        logger.info(f"  Scoring already done — {n} hits (skipping)")
        return final_tab

    # Read scoring params from config if available
    opt = getattr(cfg, 'options', None) if cfg else None
    pth = getattr(cfg, 'paths',   None) if cfg else None
    scoring_main(
        cycle_fasta, out_dir,
        db_path('Antibiotics', 'protein_fasta_protein_homolog_model.fasta'),
        db_path('Xenobiotics', 'Xenobiotics_classified.fasta'),
        db_path('ISes', '_tncentral_nointegrall_isfinder-TNs.fasta'),
        db_path('CompTns', 'Known_Tns.fasta'),
        out_dir,
        mean_of_CompTns       = getattr(opt, 'mean_of_CompTns',       5850),
        std_of_CompTns        = getattr(opt, 'std_of_CompTns',        2586),
        total_score_type      = getattr(opt, 'total_score_type',      0),
        threshold_final_score = getattr(opt, 'threshold_final_score', 0),
        max_z                 = getattr(opt, 'max_z',                 20),
        dist_type             = getattr(opt, 'dist_type',             1),
        path_of_prodigal      = getattr(pth, 'path_of_prodigal',     'prodigal'),
        path_of_blastn        = getattr(pth, 'path_of_blastn',       'blastn'),
        path_of_makeblastdb   = getattr(pth, 'path_of_makeblastdb',  'makeblastdb'),
        path_of_blastx        = getattr(pth, 'path_of_blastx',       'blastx'),
        path_of_blastp        = getattr(pth, 'path_of_blastp',       'blastp'),
        logger_name           = 'picota_complete',
    )
    if os.path.exists(final_tab):
        with open(final_tab) as _f:
            n = sum(1 for _ in _f) - 1
    else:
        n = 0
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
                      fastq_dump_path: str, logger) -> tuple:
    """Download SRA FASTQ files if not already present.

    Returns tuple: (list of valid FASTQ paths, organism name)
    """
    from src.sra_download import run_sra_down

    expected = [os.path.join(raw_dir, f"{short_id}_{i}.fastq") for i in (1, 2)]

    def _valid_fastqs(paths):
        return [f for f in paths if os.path.exists(f) and os.path.getsize(f) > 0]

    current = _valid_fastqs(expected)

    if len(current) == 2:
        logger.info(f"  FASTQ already present and valid — skipping download")
        sra_organism = query_sra_organism(short_id)  # cached lookup — fast
        return current, sra_organism

    sra_organism = query_sra_organism(short_id)
    if sra_organism and sra_organism != 'Unknown':
        logger.info(f"  SRA organism: {sra_organism}")
        os.makedirs(sra_dir, exist_ok=True)
        with open(os.path.join(sra_dir, 'sra_organism.txt'), 'w', encoding='utf-8') as fh:
            fh.write(sra_organism + '\n')

    logger.info(f"  Downloading {short_id} from SRA ...")
    run_sra_down(short_id, raw_dir, sra_dir, fastq_dump_path,
                 keep_sra_file=True, the_force=False, logger_name='picota_complete')

    current = _valid_fastqs(expected)
    if len(current) != 2:
        logger.warning(f"  FASTQ download incomplete for {short_id}; retrying with force")
        # Clear possibly-corrupted fastq files and retry
        for f in expected:
            try:
                if os.path.exists(f):
                    os.remove(f)
            except Exception:
                pass
        run_sra_down(short_id, raw_dir, sra_dir, fastq_dump_path,
                     keep_sra_file=True, the_force=True, logger_name='picota_complete')
        current = _valid_fastqs(expected)

    if len(current) != 2:
        logger.error(f"  Unable to obtain both FASTQ files for {short_id}")

    return current, sra_organism


def step_assembly(short_id: str, raw_files: list, asm_dir: str,
                  threads: int, logger, cfg) -> list:
    """Run MEGAHIT assembly if no GFA exists yet. Returns list of GFA paths."""
    from src.assembly import assembly_main

    gfa_files = [f for f in (Path(asm_dir).glob('*.gfa'))]
    if gfa_files:
        logger.info(f"  Assembly already exists — skipping")

        # Kullanıcıya skorları ve en iyi GFA'yı göster
        try:
            from src.assembly import process_gfa_files
            best = process_gfa_files([str(f) for f in gfa_files], cfg.paths.path_of_bandage)
            if best:
                logger.info(f"  Existing best GFA: {best}")
                print(f"Existing best GFA: {best}")
        except Exception as e:
            logger.warning(f"  Unable to evaluate existing GFAs: {e}")

        return [str(f) for f in gfa_files]

    assembly_kmer_list = cfg.paths.assembly_k_mer_list
    assembly_keep_temp_files = cfg.paths.keep_temp_files
    assembly_path_of_spades = cfg.paths.path_of_spades
    assembly_path_of_fastp = cfg.paths.path_of_fastp
    assembly_skip_filtering = cfg.paths.skip_filtering
    assembler_type = cfg.paths.assembler_type
    assembly_path_of_megahit = cfg.paths.path_of_megahit
    gfa_tools_path = cfg.paths.gfa_tools_path
    path_of_bandage = cfg.paths.path_of_bandage

    logger.info(f"  Running assembly for {short_id} ...")
    assembly_main(
        short_id, raw_files, asm_dir,
        threads, assembly_kmer_list,
        getattr(cfg.paths, 'quiet', True),
        assembly_keep_temp_files,
        assembly_path_of_spades, assembly_path_of_fastp, assembly_skip_filtering,
        assembler_type, assembly_path_of_megahit,
        gfa_tools_path, path_of_bandage, 'picota_complete'
    )
    return [str(f) for f in Path(asm_dir).glob('*.gfa')]


# ─── Main pipeline ───────────────────────────────────────────────────────────

def _process_sample(short_id, long_id, output_path, gfa_mode,
                    missing_tools, long_read_threads, logger, cfg):
    """Run the full pipeline for a single sample. Returns list of enriched rows."""
    STEPS = 5
    logger.info("\n" + "═" * BANNER_WIDTH)
    logger.info(f"  Sample {short_id}  |  long-read: {long_id or 'none'}")
    logger.info("═" * BANNER_WIDTH)

    sample_dir = output_path / short_id
    sample_dir.mkdir(exist_ok=True)

    scoring_dir = sample_dir / 'scoring'
    enriched_csv_path = scoring_dir / 'picota_enriched.csv'

    required_columns = {'SRA_Organism', 'ARO_ID', 'ARO_Name', 'ARO_Drug_Class'}
    if enriched_csv_path.exists():
        try:
            with open(enriched_csv_path, newline='', encoding='utf-8') as fh:
                existing_cols = set(csv.DictReader(fh).fieldnames or [])
        except Exception:
            existing_cols = set()

        if not required_columns.issubset(existing_cols):
            logger.info(f"  ⤼ Sample {short_id} enriched file present but missing new columns; rerunning pipeline")
        else:
            logger.info(f"  ⤼ Sample {short_id} already has enriched output with required columns; skipping full pipeline")
            return step_load_enriched(str(scoring_dir), logger)

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

        # Önce FASTQ indir ve her koşulda assembly işlemini yeniden çalıştır (yeni sonuç üret).
        fastq_dump = getattr(getattr(cfg, 'paths', None), 'fastq_dump', 'parallel-fastq-dump') if cfg else 'parallel-fastq-dump'
        raw_files, sra_organism = step_sra_download(short_id, raw_dir, sra_dir,
                                                    fastq_dump, logger)

        if not raw_files:
            logger.warning(f"  ⤼ No FASTQ files after download, checking assembly folder")
            existing_gfa = list(Path(asm_dir).glob('*.gfa'))
            if existing_gfa:
                gfa_file = str(existing_gfa[0])
                logger.info(f"  ⤼ Found existing assembly despite missing FASTQ, using existing {gfa_file}")
            else:
                logger.error(f"  ✗ No FASTQ files found after download and no existing assembly — aborting {short_id}")
                return []
        else:
            # Assembly yeniden hesaplanıp gfa üzerine yazılacak.
            asm_threads = getattr(getattr(cfg, 'paths', None), 'assembly_threads', 4) if cfg else 4
            gfa_list = step_assembly(short_id, raw_files, asm_dir,
                                     threads=asm_threads, logger=logger, cfg=cfg)
            if not gfa_list:
                logger.error(f"  ✗ Assembly produced no GFA — aborting {short_id}")
                return []
            gfa_file = gfa_list[0]

            # Delete raw FASTQs after successful assembly if config says so
            delete_fastq = getattr(getattr(cfg, 'options', None), 'delete_fastq_files', False) if cfg else False
            if delete_fastq:
                for f in raw_files:
                    try:
                        os.remove(f)
                    except Exception:
                        pass
                logger.info(f"  ✓ Raw FASTQs deleted (delete_fastq_files=true)")

    logger.info(f"  ✓ GFA ready: {gfa_file}  ({time.time()-t0:.1f}s)")

    # ── Step 2/5: Cycle detection ─────────────────────────────────────────────
    logger.info(f"\n  [2/{STEPS}] Cycle Detection")
    t0 = time.time()
    cycle_fasta = str(sample_dir / f'{short_id}_cycles.fasta')
    n_cycles = step_cycle_detection(gfa_file, cycle_fasta, logger, cfg)
    if n_cycles == 0:
        logger.warning(f"  ⚠  No cycles detected — skipping scoring")
        return []
    logger.info(f"  ✓ {n_cycles} candidate cycles  ({time.time()-t0:.1f}s)")

    # ── Step 3/5: BLAST scoring ───────────────────────────────────────────────
    logger.info(f"\n  [3/{STEPS}] BLAST Scoring")
    score_dir    = str(sample_dir / 'scoring')
    final_tab    = os.path.join(score_dir, 'picota_final_tab')
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
        final_tab = step_scoring(cycle_fasta, score_dir, logger, cfg)
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
        long_raw_dir  = output_path / 'raw_long' / long_id
        long_sra_dir  = output_path / 'sra_long'  / long_id
        long_fastq    = long_raw_dir / f'{long_id}_1.fastq'

        # Download long-read FASTQ if not already present
        if not (long_fastq.exists() and long_fastq.stat().st_size > 0):
            logger.info(f"  Downloading long-read {long_id} ...")
            try:
                from src.sra_download import run_sra_down
                long_raw_dir.mkdir(parents=True, exist_ok=True)
                long_sra_dir.mkdir(parents=True, exist_ok=True)
                fastq_dump = getattr(cfg.paths, 'fastq_dump', 'parallel-fastq-dump') if cfg else 'parallel-fastq-dump'
                run_sra_down(long_id, str(long_raw_dir), str(long_sra_dir), fastq_dump)
                logger.info(f"  Downloaded {long_id}")
            except Exception as e:
                logger.warning(f"  Long-read download failed: {e}")

        if long_fastq.exists() and long_fastq.stat().st_size > 0:
            t0 = time.time()
            map_dir = str(output_path / 'mapping' / long_id)
            n_mapped = 0
            mapping_threads = getattr(getattr(cfg, 'options', None), 'mapping_threads', long_read_threads) if cfg else long_read_threads
            for fa in annotated:
                try:
                    sorted_bam = step_long_read_mapping(fa, str(long_fastq), map_dir,
                                                        long_id, mapping_threads, logger)
                    # BAM analysis — same as picota_testv3.py
                    from src.bam_analyse import bam_file_analyze
                    from src.analyze_blocksv3 import analyze_blocks
                    fa_base = os.path.basename(fa)
                    out_full    = os.path.join(map_dir, f"{long_id}_{fa_base}_full")
                    out_partial = os.path.join(map_dir, f"{long_id}_{fa_base}_partial")
                    bam_file_analyze(sorted_bam, out_full, out_partial)
                    analyze_blocks(out_partial, out_full, map_dir, fa_base, cfg)
                    n_mapped += 1
                except Exception as e:
                    logger.error(f"  Mapping/BAM error: {e}")
            logger.info(f"  ✓ {n_mapped} FASTA(s) mapped + analysed  ({time.time()-t0:.1f}s)")
        else:
            logger.info(f"  ⤼ Long-read FASTQ not available after download attempt: {long_fastq}")
    else:
        logger.info(f"  ⤼ Long-read validation skipped (no long-read ID or tools missing)")

    return enriched


def run_pipeline(sra_list_file: str, output_dir: str, gfa_mode: bool = False,
                 long_read_threads: int = 4, cfg=None):
    """Execute the complete PICOTA pipeline for all samples in sra_list_file."""

    _here = os.path.dirname(os.path.abspath(__file__))
    if cfg is None:
        cfg = load_config(os.path.join(_here, 'picota', 'config.yaml'))

    # Resolve relative tool paths against the script directory (_here).
    # Config paths like "./picota/tools/..." are written relative to the
    # repo root (same dir as this script), so we resolve against _here.
    # Bare command names ("megahit", "fastp") are left untouched.
    for attr in ('gfa_tools_path', 'path_of_bandage', 'path_of_spades',
                 'path_of_megahit', 'path_of_fastp'):
        val = getattr(cfg.paths, attr, '')
        if val and not os.path.isabs(val) and (val.startswith('./') or os.sep in val):
            setattr(cfg.paths, attr, os.path.normpath(os.path.join(_here, val)))

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
                                   gfa_mode, missing, long_read_threads, logger,
                                   cfg)
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

        # Önce varsa eski sonuçları sil, sonra yaz (overwrite garanti)
        if csv_out.exists():
            csv_out.unlink()
        if json_out.exists():
            json_out.unlink()

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

  # Run with PICOTA code from a different directory:
  conda run -n evobiomig python3 /path/to/test_complete_pipeline.py \\
      --picota_dir /path/to/picota --output /data/picota_results/
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
    parser.add_argument('--picota_dir', '-p',
                        default=_here,
                        help='Directory containing PICOTA source code (default: script dir)')
    parser.add_argument('--config', '-c',
                        default=None,
                        help='Path to config YAML (default: picota_dir/picota/config.yaml)')

    args = parser.parse_args()

    # Set up paths based on picota_dir
    setup_picota_paths(args.picota_dir)

    # Load config if given
    cfg = None
    if args.config:
        cfg = load_config(args.config)

    # Config can supply defaults for --sra_list and --output;
    # CLI args take priority when explicitly provided by the user.
    # argparse sets defaults, so we check if the user actually passed the flag.
    _defaults = parser.parse_args([])  # parse with no args to get bare defaults
    if cfg is not None:
        cfg_paths = getattr(cfg, 'paths', cfg)
        if args.sra_list == _defaults.sra_list:
            cfg_sra = getattr(cfg_paths, 'sra_id_file', None)
            if cfg_sra:
                args.sra_list = cfg_sra
        if args.output == _defaults.output:
            cfg_out = getattr(cfg_paths, 'outdir', None)
            if cfg_out:
                args.output = cfg_out
        if args.threads == _defaults.threads:
            cfg_threads = getattr(cfg_paths, 'assembly_threads', None)
            if cfg_threads:
                args.threads = int(cfg_threads)

    # --gfa_mode: always use a single testNitro entry (ignores real SRA IDs)
    if args.gfa_mode:
        gfa_csv = os.path.join(args.picota_dir, 'picota', 'testNitro_sra_ids.csv')
        with open(gfa_csv, 'w') as fh:
            fh.write('sra_short_id,sra_long_id\ntestNitro,-\n')
        args.sra_list = gfa_csv

    run_pipeline(args.sra_list, args.output,
                 gfa_mode=args.gfa_mode,
                 long_read_threads=args.threads,
                 cfg=cfg)


if __name__ == '__main__':
    main()
