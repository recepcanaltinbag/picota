#!/usr/bin/env python3
"""
PICOTA - Pipeline for Identification of Composite Transposons from Assembly graphs

A bioinformatics pipeline for detecting composite transposons in bacterial genomes
using assembly graph analysis and homology-based scoring.

Author: Recep Canaltinbag
Version: 1.0.0-rc1
License: MIT
"""

import argparse
import sys
import os
from pathlib import Path
import yaml
import logging
import shutil

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

def setup_logging(log_file: str = None, level: str = "INFO"):
    """Setup logging configuration"""
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file) if log_file else logging.NullHandler()
        ]
    )
    return logging.getLogger("picota")

def load_config(config_path: str) -> dict:
    """Load configuration from YAML file"""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file not found: {config_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing configuration file: {e}")
        sys.exit(1)

def validate_environment():
    """Check if required tools are available"""
    required_tools = ['blastn', 'blastp', 'prodigal', 'spades.py', 'megahit']
    missing_tools = []

    for tool in required_tools:
        if not shutil.which(tool):
            missing_tools.append(tool)

    if missing_tools:
        print("Warning: Missing required tools:", ", ".join(missing_tools))
        print("Please install via conda: conda install -c bioconda", " ".join(missing_tools))
        return False
    return True

def run_sra_download(args, config):
    """Download data from NCBI SRA"""
    try:
        from src.sra_download import run_sra_down
        logger = setup_logging(config.get('logging', {}).get('log_file'))
        run_sra_down(
            sra_acc_file=args.sra_file,
            out_dir=args.output,
            sra_folder=os.path.join(args.output, 'sra_files'),
            path_of_fastq_dump=config['paths'].get('fastq_dump', 'parallel-fastq-dump'),
            delete_fastq=config['options'].get('delete_fastq_files', True),
            logger_name="picota_sra"
        )
        logger.info("SRA download completed successfully")
    except Exception as e:
        print(f"Error in SRA download: {e}")
        sys.exit(1)

def run_assembly(args, config):
    """Run genome assembly"""
    try:
        from src.assembly import run_assembly_pipeline
        logger = setup_logging(config.get('logging', {}).get('log_file'))
        run_assembly_pipeline(
            fastq_path=args.fastq,
            output_dir=args.output,
            threads=args.threads,
            config=config
        )
        logger.info("Assembly completed successfully")
    except Exception as e:
        print(f"Error in assembly: {e}")
        sys.exit(1)

def run_analysis(args, config):
    """Run cycle detection and scoring analysis"""
    try:
        from src.cycle_finder import GraphWork
        from src.scoringv3_blast import run_scoring_pipeline
        from src.analyze_blocks import analyze_gfa_file

        logger = setup_logging(config.get('logging', {}).get('log_file'))

        # Validate input files
        if not os.path.exists(args.gfa):
            print(f"Error: GFA file not found: {args.gfa}")
            sys.exit(1)

        # Create output directory
        os.makedirs(args.output, exist_ok=True)

        # Run analysis
        logger.info(f"Starting analysis on {args.gfa}")

        # Parse GFA and find cycles
        cycles = analyze_gfa_file(args.gfa, config)

        if not cycles:
            logger.warning("No cycles detected in assembly graph")
            return

        # Score cycles
        results = run_scoring_pipeline(
            cycle_data=cycles,
            output_dir=args.output,
            config=config
        )

        # Save results
        results_file = os.path.join(args.output, 'transposon_candidates.tsv')
        results.to_csv(results_file, sep='\t', index=False)

        logger.info(f"Analysis completed. Results saved to {results_file}")
        print(f"Found {len(results)} potential composite transposons")

    except Exception as e:
        print(f"Error in analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def run_scoring(args, config):
    """Run scoring on existing cycle data"""
    try:
        from src.scoringv3_blast import run_scoring_pipeline
        logger = setup_logging(config.get('logging', {}).get('log_file'))

        if not os.path.exists(args.cycle_folder):
            print(f"Error: Cycle folder not found: {args.cycle_folder}")
            sys.exit(1)

        results = run_scoring_pipeline(
            cycle_folder=args.cycle_folder,
            output_dir=args.output,
            config=config
        )

        results_file = os.path.join(args.output, 'transposon_candidates.tsv')
        results.to_csv(results_file, sep='\t', index=False)

        logger.info(f"Scoring completed. Results saved to {results_file}")

    except Exception as e:
        print(f"Error in scoring: {e}")
        sys.exit(1)

def run_all_pipeline(args, config):
    """Run complete pipeline from SRA to results"""
    try:
        logger = setup_logging(config.get('logging', {}).get('log_file'))
        logger.info("Starting complete PICOTA pipeline")

        # Step 1: SRA download
        if args.sra:
            logger.info("Step 1: Downloading SRA data")
            run_sra_download(args, config)

        # Step 2: Assembly
        fastq_path = args.fastq or os.path.join(args.output, 'fastq')
        if os.path.exists(fastq_path):
            logger.info("Step 2: Running genome assembly")
            assembly_args = argparse.Namespace(
                fastq=fastq_path,
                output=os.path.join(args.output, 'assembly'),
                threads=args.threads
            )
            run_assembly(assembly_args, config)

        # Step 3: Analysis
        gfa_path = args.gfa or os.path.join(args.output, 'assembly', 'assembly.gfa')
        if os.path.exists(gfa_path):
            logger.info("Step 3: Running transposon analysis")
            analysis_args = argparse.Namespace(
                gfa=gfa_path,
                output=os.path.join(args.output, 'results'),
                threads=args.threads
            )
            run_analysis(analysis_args, config)

        logger.info("Complete pipeline finished successfully")

    except Exception as e:
        print(f"Error in complete pipeline: {e}")
        sys.exit(1)

def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description="PICOTA - Pipeline for Identification of Composite Transposons",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download and analyze from SRA
  python picota.py all --sra SRR11362851 --output results/

  # Analyze existing GFA file
  python picota.py analysis --gfa assembly.gfa --output results/

  # Download SRA data only
  python picota.py sra_download --sra SRR11362851 --output data/

  # Run assembly on FASTQ files
  python picota.py assembly --fastq reads.fq --output assembly/ --threads 8

  # Score existing cycles
  python picota.py scoring --cycle_folder cycles/ --output results/
        """
    )

    parser.add_argument('--config', '-c',
                       default='picota/config.yaml',
                       help='Configuration file path (default: picota/config.yaml)')

    parser.add_argument('--threads', '-t', type=int, default=4,
                       help='Number of threads to use (default: 4)')

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # SRA download command
    sra_parser = subparsers.add_parser('sra_download', help='Download data from NCBI SRA')
    sra_parser.add_argument('--sra', required=True, help='SRA accession number')
    sra_parser.add_argument('--output', '-o', required=True, help='Output directory')

    # Assembly command
    assembly_parser = subparsers.add_parser('assembly', help='Run genome assembly')
    assembly_parser.add_argument('--fastq', required=True, help='Input FASTQ file/directory')
    assembly_parser.add_argument('--output', '-o', required=True, help='Output directory')

    # Analysis command
    analysis_parser = subparsers.add_parser('analysis', help='Detect and score transposons')
    analysis_parser.add_argument('--gfa', required=True, help='Input GFA file')
    analysis_parser.add_argument('--output', '-o', required=True, help='Output directory')
    analysis_parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads')

    # Scoring command
    scoring_parser = subparsers.add_parser('scoring', help='Score existing cycles')
    scoring_parser.add_argument('--cycle_folder', required=True, help='Cycle data folder')
    scoring_parser.add_argument('--output', '-o', required=True, help='Output directory')

    # Complete pipeline
    all_parser = subparsers.add_parser('all', help='Run complete pipeline')
    all_parser.add_argument('--sra', help='SRA accession number')
    all_parser.add_argument('--fastq', help='Input FASTQ file/directory')
    all_parser.add_argument('--gfa', help='Input GFA file')
    all_parser.add_argument('--output', '-o', required=True, help='Output directory')

    # Parse arguments
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Load configuration
    config = load_config(args.config)

    # Validate environment
    if not validate_environment():
        print("Warning: Some tools are missing. Pipeline may fail.")

    # Execute command
    try:
        if args.command == 'sra_download':
            run_sra_download(args, config)
        elif args.command == 'assembly':
            run_assembly(args, config)
        elif args.command == 'analysis':
            run_analysis(args, config)
        elif args.command == 'scoring':
            run_scoring(args, config)
        elif args.command == 'all':
            run_all_pipeline(args, config)

    except KeyboardInterrupt:
        print("\nInterrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
