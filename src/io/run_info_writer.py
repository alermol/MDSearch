"""Run information file output operations."""

import datetime
import platform
import subprocess
from pathlib import Path
from typing import List

from ..core.vcf_parser import VCFData
from ..utils.memory_monitor import MemoryMonitor
from ..version import __version__ as mdsearch_version

__all__ = ["RunInfoWriter"]


class RunInfoWriter:
    """Handles run information file output."""

    def __init__(self, memory_monitor: MemoryMonitor):
        """Initialize run info writer with memory monitor.

        Args:
            memory_monitor: MemoryMonitor instance for tracking memory usage

        Example:
            >>> from src.utils.memory_monitor import MemoryMonitor
            >>> from src.utils.logging_setup import setup_logger
            >>> logger = setup_logger("test")
            >>> monitor = MemoryMonitor(logger)
            >>> writer = RunInfoWriter(monitor)
        """
        self.memory_monitor = memory_monitor

    def write_run_info(
        self, 
        output_prefix: Path, 
        vcf_data: VCFData, 
        snp_sets: List[List[str]],
        config_data: dict
    ) -> None:
        """Write comprehensive run information to a file.

        Args:
            output_prefix: Path to output directory
            vcf_data: Parsed VCF data containing sample and SNP information
            snp_sets: List of SNP sets that were generated
            config_data: Dictionary containing configuration information

        Example:
            >>> from pathlib import Path
            >>> writer = RunInfoWriter(memory_monitor)
            >>> config = {"ploidy": 2, "min_distance": 3, "convert_het": False}
            >>> writer.write_run_info(Path("output"), vcf_data, snp_sets, config)
            >>> # Creates output/run_info.txt with comprehensive run information
        """
        run_info_path = output_prefix / "run_info.txt"

        # Get git commit hash if available
        try:
            git_commit = (
                subprocess.run(
                    ["git", "rev-parse", "--short", "HEAD"],
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.DEVNULL,
                    text=True,
                ).stdout.strip()
                or "unknown"
            )
        except Exception:
            git_commit = "unknown"

        # Get current timestamp
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC")

        # Get comprehensive memory information (after run completion)
        memory_summary = self.memory_monitor.get_memory_summary()
        threshold_info = self.memory_monitor.get_threshold_info()
        estimated_matrix_mb = self.memory_monitor.estimate_matrix_memory_mb(
            len(vcf_data.snp_genotypes), len(vcf_data.headers.samples)
        )

        # Prepare run info content
        lines = [
            "MDSearch Run Information",
            "=======================",
            "",
            f"Version: {mdsearch_version}",
            f"Git commit: {git_commit}",
            f"Python version: {platform.python_version()}",
            f"Platform: {platform.platform()}",
            "",
            f"Run timestamp: {timestamp}",
            "",
            "System Memory Information:",
            f"  Current process memory: {memory_summary['current_mb']:.1f} MB",
            f"  Peak process memory: {memory_summary['peak_mb']:.1f} MB",
            f"  Available system memory: {memory_summary['available_mb']:.1f} MB",
            f"  Total system memory: {memory_summary['total_mb']:.1f} MB",
            f"  Memory warning threshold: {memory_summary['warning_threshold_mb']:.1f} MB ({threshold_info['warning_percent']:.0f}%)",
            f"  Critical threshold: {memory_summary['critical_threshold_mb']:.1f} MB ({threshold_info['critical_percent']:.0f}%)",
            f"  Estimated genotype matrix memory: {estimated_matrix_mb:.1f} MB",
            "",
            "Input Configuration:",
            f"  Input VCF: {config_data['input_vcf']}",
            f"  Input format: {config_data['input_format']}",
            f"  Output directory: {config_data['output_prefix']}",
            f"  Output format: {config_data['output_format']}",
            f"  Ploidy: {config_data['ploidy']}",
            f"  Minimum distance: {config_data['min_distance']}",
            f"  Convert heterozygous: {config_data['convert_het']}",
            f"  Target number of sets: {config_data['n_sets']}",
            f"  Weight entropy: {config_data['weight_entropy']}",
            f"  Weight MAF: {config_data['weight_maf']}",
            f"  Verbose: {config_data['verbose']}",
            f"  Log level: {config_data['log_level'] or 'default'}",
            f"  Log format: {config_data['log_format']}",
            "",
            "Input Data Summary:",
            f"  Number of samples: {len(vcf_data.headers.samples)}",
            f"  Number of SNPs: {len(vcf_data.snp_genotypes)}",
            f"  Chromosomes: {', '.join(sorted(vcf_data.headers.contigs))}",
            "",
            "Output Summary:",
            f"  Number of SNP sets found: {len(snp_sets)}",
        ]

        # Add details for each SNP set
        for i, snp_set in enumerate(snp_sets, 1):
            lines.append(f"  Set {i}: {len(snp_set)} SNPs")

        # Write to file
        with open(run_info_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

        return run_info_path
