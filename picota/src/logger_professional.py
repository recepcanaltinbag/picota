"""
Professional logging and reporting module for PICOTA pipeline

Provides structured logging, progress tracking, and formatted output for
bioinformatics analysis workflows.
"""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional
import json


class PICOTAFormatter(logging.Formatter):
    """Custom formatter for PICOTA logs with color support and structure"""

    COLORS = {
        'DEBUG': '\033[36m',    # Cyan
        'INFO': '\033[32m',     # Green
        'WARNING': '\033[33m',  # Yellow
        'ERROR': '\033[31m',    # Red
        'CRITICAL': '\033[41m', # Red background
    }
    RESET = '\033[0m'

    def format(self, record):
        """Format log record with color"""
        log_color = self.COLORS.get(record.levelname, self.RESET)
        timestamp = datetime.fromtimestamp(record.created).strftime('%Y-%m-%d %H:%M:%S')

        # Format message
        if record.levelname == 'INFO':
            formatted = f"{log_color}[{timestamp}] [{record.name}] {record.getMessage()}{self.RESET}"
        else:
            formatted = f"{log_color}[{timestamp}] [{record.levelname}] [{record.name}] {record.getMessage()}{self.RESET}"

        return formatted


class PICOTALogger:
    """Professional logger for PICOTA with file and console output"""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if not hasattr(self, '_initialized'):
            self.logger = logging.getLogger('PICOTA')
            self.logger.setLevel(logging.DEBUG)
            self._initialized = True
            self.log_file = None

    @staticmethod
    def setup(log_file: Optional[str] = None, level: str = 'INFO') -> logging.Logger:
        """Setup logger with console and file handlers"""
        logger_instance = PICOTALogger()
        logger = logger_instance.logger

        # Clear existing handlers
        logger.handlers = []

        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, level.upper()))
        console_handler.setFormatter(PICOTAFormatter())
        logger.addHandler(console_handler)

        # File handler
        if log_file:
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(logging.DEBUG)
            file_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
            logger_instance.log_file = log_file

        return logger

    @staticmethod
    def get_logger(name: str = 'PICOTA') -> logging.Logger:
        """Get logger instance"""
        return logging.getLogger(name)


class AnalysisProgress:
    """Track and report analysis progress"""

    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.steps = {}
        self.current_step = None
        self.start_time = datetime.now()

    def add_step(self, step_id: str, description: str, substeps: int = 0):
        """Add an analysis step"""
        self.steps[step_id] = {
            'description': description,
            'substeps': substeps,
            'completed': 0,
            'status': 'pending'
        }

    def start_step(self, step_id: str):
        """Mark step as started"""
        if step_id in self.steps:
            self.current_step = step_id
            self.steps[step_id]['status'] = 'in_progress'
            self.logger.info(f"▶ Starting: {self.steps[step_id]['description']}")

    def complete_step(self, step_id: str, message: str = ""):
        """Mark step as completed"""
        if step_id in self.steps:
            self.steps[step_id]['status'] = 'completed'
            self.steps[step_id]['completed'] += 1
            msg = f"✓ Completed: {self.steps[step_id]['description']}"
            if message:
                msg += f" - {message}"
            self.logger.info(msg)
            self.current_step = None

    def report_summary(self):
        """Report overall progress summary"""
        elapsed = (datetime.now() - self.start_time).total_seconds()
        self.logger.info("=" * 70)
        self.logger.info("ANALYSIS SUMMARY")
        self.logger.info("=" * 70)

        for step_id, step_info in self.steps.items():
            status_symbol = "✓" if step_info['status'] == 'completed' else "✗"
            self.logger.info(f"{status_symbol} {step_info['description']}")

        self.logger.info("=" * 70)
        self.logger.info(f"Total elapsed time: {elapsed:.2f} seconds")


class ResultsFormatter:
    """Format and export analysis results"""

    @staticmethod
    def format_ct_metadata(ct_data: dict) -> dict:
        """Format composite transposon metadata for output"""
        return {
            'CT_ID': ct_data.get('ct_id', 'CT000'),
            'Category': ct_data.get('category', 'Unknown'),  # 'Novel' or 'Known'
            'CT_Length_bp': ct_data.get('ct_length', 0),
            'IS_Family': ct_data.get('is_family', 'Unknown'),
            'IS_Group': ct_data.get('is_group', 'Unknown'),
            'IS_Length_bp': ct_data.get('is_length', 0),
            'Antibiotic_Classes': '; '.join(ct_data.get('antibiotic_classes', [])),
            'Resistance_Genes': '; '.join(ct_data.get('resistance_genes', [])),
            'Score': f"{ct_data.get('score', 0):.2f}",
            'Organism': ct_data.get('organism', 'Unknown'),
            'Novel': ct_data.get('novel', False),
        }

    @staticmethod
    def to_csv(results: list, output_file: str) -> None:
        """Export results to CSV format"""
        import csv
        if not results:
            return

        fieldnames = list(results[0].keys())

        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

    @staticmethod
    def to_TSV(results: list, output_file: str) -> None:
        """Export results to TSV format"""
        import csv
        if not results:
            return

        fieldnames = list(results[0].keys())

        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(results)

    @staticmethod
    def to_json(results: list, output_file: str) -> None:
        """Export results to JSON format"""
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)

    @staticmethod
    def print_summary(results: list, logger: logging.Logger):
        """Print summary statistics"""
        if not results:
            logger.warning("No results to summarize")
            return

        logger.info("=" * 70)
        logger.info("RESULTS SUMMARY")
        logger.info("=" * 70)

        # Count by category
        novel_count = sum(1 for r in results if r.get('Category') == 'Novel')
        known_count = sum(1 for r in results if r.get('Category') == 'Known')

        logger.info(f"Total Composite Transposons Detected: {len(results)}")
        logger.info(f"  ├─ Novel: {novel_count}")
        logger.info(f"  └─ Known: {known_count}")

        # Antibiotic resistance summary
        all_classes = set()
        for r in results:
            classes = r.get('Antibiotic_Classes', '').split('; ')
            all_classes.update(filter(None, classes))

        if all_classes:
            logger.info(f"Antibiotic Resistance Classes Found: {', '.join(sorted(all_classes))}")

        # Top scoring CTs
        top_cts = sorted(results, key=lambda x: float(x.get('Score', 0)), reverse=True)[:5]
        if top_cts:
            logger.info("Top Scoring Composite Transposons:")
            for i, ct in enumerate(top_cts, 1):
                logger.info(f"  {i}. {ct.get('CT_ID')} - Score: {ct.get('Score')}")

        logger.info("=" * 70)
