#!/usr/bin/env python3
"""
PICOTA - Pipeline for Identification of Composite Transposons from Assembly graphs

A bioinformatics pipeline for detecting composite transposons in bacterial genomes
using assembly graph analysis and homology-based scoring.

Every subcommand is a thin front end over the stage functions in
picota_testv3.py, so the CLI and the config-driven pipeline run the same code.

Author: Recep Canaltinbag
Version: 1.1.0
License: MIT
"""

import argparse
import os
import shutil
import sys

# The stage modules import as "src.X" and "picota_testv3", both relative to
# this directory, whether run from a checkout or from an installed package.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import picota_testv3 as pipeline  # noqa: E402
from src.config_loader import load_config  # noqa: E402
from src.logger_setup import setup_logger_from_config  # noqa: E402

__version__ = "1.1.0"

def required_tools(command, cfg):
    """The executables `command` will call, as named in the config."""
    download = [cfg.paths.fastq_dump]
    assembler = (cfg.paths.path_of_megahit if cfg.paths.assembler_type == 'megahit'
                 else cfg.paths.path_of_spades)
    assembly = ([] if cfg.paths.skip_filtering else [cfg.paths.path_of_fastp]) + [assembler]
    scoring = [cfg.options.path_of_prodigal, cfg.options.path_of_blastn,
               cfg.options.path_of_blastp, cfg.options.path_of_makeblastdb]
    return {
        'sra_download': download,
        'assembly':     assembly,
        'analysis':     scoring,
        'scoring':      scoring,
        'all':          download + assembly + scoring,
    }[command]


def missing_tools(command, cfg):
    """Return the tools `command` needs that are not on PATH."""
    return [t for t in required_tools(command, cfg) if t and not shutil.which(t)]


def prepare(args):
    """Load the config, apply command-line overrides and start the logger."""
    if not os.path.exists(args.config):
        sys.exit(f"Error: configuration file not found: {args.config}")
    cfg = load_config(args.config)
    cfg.paths.outdir = args.output
    if args.threads is not None:
        cfg.paths.assembly_threads = args.threads
        cfg.options.mapping_threads = args.threads
        os.environ['PICOTA_BLAST_THREADS'] = str(args.threads)
    os.makedirs(args.output, exist_ok=True)
    pipeline.logger = setup_logger_from_config(cfg)
    return cfg


def cmd_sra_download(args, cfg):
    files = pipeline.run_sra_download(
        args.sra, args.output, os.path.join(args.output, 'sra', args.sra),
        cfg.paths.fastq_dump, cfg.logging.logger_name)
    if not files:
        sys.exit(f"Error: no FASTQ files were produced for {args.sra}")
    return files


def cmd_assembly(args, cfg):
    name = args.name or os.path.basename(args.fastq[0]).split('_')[0].split('.')[0]
    gfa_files = pipeline.run_assembly(name, args.fastq, args.output, cfg,
                                      cfg.logging.logger_name)
    if not gfa_files:
        sys.exit(f"Error: assembly produced no GFA in {args.output}")
    return gfa_files


def cmd_scoring(args, cfg):
    if not os.path.exists(args.cycle_folder):
        sys.exit(f"Error: cycle file or folder not found: {args.cycle_folder}")
    name = os.path.splitext(os.path.basename(os.path.normpath(args.cycle_folder)))[0]
    return pipeline.run_scoring(name, args.cycle_folder, args.output, cfg)


def cmd_analysis(args, cfg):
    """Cycle detection, scoring and transposon/cargo split on one GFA."""
    if not os.path.exists(args.gfa):
        sys.exit(f"Error: GFA file not found: {args.gfa}")
    name = os.path.splitext(os.path.basename(args.gfa))[0]
    cycle_dir = os.path.join(args.output, 'cycles')
    scoring_dir = os.path.join(args.output, 'scoring')
    annot_dir = os.path.join(args.output, 'annot')
    for d in (cycle_dir, scoring_dir, annot_dir):
        os.makedirs(d, exist_ok=True)

    cycle_file = os.path.join(cycle_dir, f"{name}.fasta")
    pipeline.run_cycle_analysis(name, args.gfa, cycle_file, cfg)
    if not os.path.exists(cycle_file):
        pipeline.logger.warning("No cycles detected in the assembly graph")
        return None

    final_tab = pipeline.run_scoring(name, cycle_file, scoring_dir, cfg)
    if final_tab and os.path.exists(final_tab):
        pipeline.split_cycles_from_picota(final_tab, cycle_file, annot_dir,
                                          cfg.options.split_min_score)
    return final_tab


def cmd_all(args, cfg):
    """Run every stage, from whichever input was given."""
    if args.gfa:
        return cmd_analysis(argparse.Namespace(gfa=args.gfa, output=args.output), cfg)
    if args.fastq:
        asm_out = os.path.join(args.output, 'assembly')
        gfa_files = cmd_assembly(
            argparse.Namespace(fastq=args.fastq, output=asm_out, name=None), cfg)
        return cmd_analysis(argparse.Namespace(gfa=gfa_files[0], output=args.output), cfg)
    if args.sra:
        pipeline.process_accession(args.sra, args.long, cfg)
        return None
    for short_id, long_id in pipeline.load_sra_pairs(cfg.paths.sra_id_file):
        pipeline.process_accession(short_id, long_id, cfg)
    return None


COMMANDS = {
    'sra_download': cmd_sra_download,
    'assembly':     cmd_assembly,
    'analysis':     cmd_analysis,
    'scoring':      cmd_scoring,
    'all':          cmd_all,
}


def build_parser():
    parser = argparse.ArgumentParser(
        prog='picota',
        description="PICOTA - Pipeline for Identification of Composite Transposons",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  picota all          --sra SRR11362851 --output results/
  picota analysis     --gfa assembly.gfa --output results/
  picota sra_download --sra SRR11362851 --output data/
  picota assembly     --fastq reads_1.fastq reads_2.fastq --output assembly/ -t 8
  picota scoring      --cycle_folder cycles/ --output results/
        """
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    # Shared by every subcommand, so they may follow the subcommand name.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--config', '-c', default='picota/config.yaml',
                        help='Configuration file (default: picota/config.yaml)')
    common.add_argument('--threads', '-t', type=int, default=None,
                        help='Threads for assembly, mapping and BLAST '
                             '(default: the values in the config)')
    common.add_argument('--output', '-o', required=True, help='Output directory')

    sub = parser.add_subparsers(dest='command', metavar='COMMAND')

    p = sub.add_parser('sra_download', parents=[common], help='Download reads from NCBI SRA')
    p.add_argument('--sra', required=True, help='SRA run accession')

    p = sub.add_parser('assembly', parents=[common], help='Assemble reads into a GFA graph')
    p.add_argument('--fastq', required=True, nargs='+', help='FASTQ file(s), paired files together')
    p.add_argument('--name', help='Sample name (default: taken from the first FASTQ)')

    p = sub.add_parser('analysis', parents=[common],
                       help='Detect, score and split composite transposons in a GFA')
    p.add_argument('--gfa', required=True, help='Input GFA file')

    p = sub.add_parser('scoring', parents=[common], help='Score existing cycles')
    p.add_argument('--cycle_folder', required=True,
                   help='Cycle FASTA file, or a folder of *.fasta files')

    p = sub.add_parser('all', parents=[common],
                       help='Run every stage (no input: the samples in sra_id_file)')
    src = p.add_mutually_exclusive_group()
    src.add_argument('--sra', help='Short-read SRA run accession')
    src.add_argument('--fastq', nargs='+', help='FASTQ file(s); skips the download')
    src.add_argument('--gfa', help='GFA file; skips download and assembly')
    p.add_argument('--long', help='Long-read SRA accession for validation (with --sra)')
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 1

    cfg = prepare(args)
    missing = missing_tools(args.command, cfg)
    if missing:
        print("Warning: tools not found on PATH: " + ", ".join(missing), file=sys.stderr)

    try:
        COMMANDS[args.command](args, cfg)
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        return 130
    return 0


if __name__ == '__main__':
    sys.exit(main())
