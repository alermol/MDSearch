"""Main application coordinator for MDSearch."""

import logging
from pathlib import Path
from typing import Optional, Union, Callable
from dataclasses import dataclass
import time

from .core import VCFParser, DistanceCalculator, SNPSelector
from .io import (
    VCFWriter,
    WriteConfig,
    SummaryWriter,
    RunInfoWriter,
    StructureInfoWriter,
)
from .utils import MemoryMonitor, setup_logger
from .utils import ensure_variant_index

__all__ = ["MDSearchConfig", "MDSearchApp"]


@dataclass
class MDSearchConfig:
    """Configuration for MDSearch application."""

    input_vcf: Path
    output_prefix: Path
    ploidy: int = 2
    min_distance: int = 1
    convert_het: bool = False
    n_sets: Union[int, str] = 1
    verbose: bool = True
    log_level: Optional[str] = None
    log_format: str = "text"
    weight_entropy: float = 0.5
    weight_maf: float = 0.5
    input_format: str = "auto"
    output_format: str = "v"


class MDSearchApp:
    """Main application coordinator with separated concerns."""

    def __init__(
        self,
        config: MDSearchConfig,
        shutdown_checker: Optional[Callable[[], bool]] = None,
    ):
        """Initialize MDSearch application with configuration."""
        self.config = config
        self.shutdown_checker = shutdown_checker
        self.logger = setup_logger(
            "mdsearch", config.log_level, config.log_format, config.verbose
        )
        self.memory_monitor = MemoryMonitor(self.logger)

        # Initialize components
        self.vcf_parser = VCFParser(self.memory_monitor, self.logger, shutdown_checker)
        self.distance_calc = DistanceCalculator(
            self.memory_monitor, self.logger, shutdown_checker=shutdown_checker
        )
        self.snp_selector = SNPSelector(
            self.distance_calc,
            self.memory_monitor,
            self.logger,
            shutdown_checker,
            weight_entropy=self.config.weight_entropy,
            weight_maf=self.config.weight_maf,
        )
        self.vcf_writer = VCFWriter()
        self.summary_writer = SummaryWriter(self.distance_calc)
        self.run_info_writer = RunInfoWriter(self.memory_monitor)
        self.structure_info_writer = StructureInfoWriter(self.config.output_format)

        # Initial memory check
        self.memory_monitor.check_memory_and_warn("initialization")

    def run(self) -> None:
        """Execute the complete MDSearch pipeline."""
        start_time = time.time()

        # Ensure index exists for compressed inputs
        try:
            ensure_variant_index(
                self.config.input_vcf, self.config.input_format, self.logger
            )
        except Exception as e:
            if self.logger.isEnabledFor(logging.ERROR):
                self.logger.error(f"Failed to ensure variant index: {e}")
            raise

        # Parse and validate VCF
        vcf_data = self.vcf_parser.parse_and_validate(
            self.config.input_vcf, self.config.ploidy, self.config.convert_het
        )

        # Search for optimal SNP sets
        snp_sets = self.snp_selector.search_optimal_sets(
            vcf_data,
            self.config.min_distance,
            self.config.n_sets,
        )

        # Write output files
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Writing selected SNPs in VCF...")

        write_config = WriteConfig(
            ploidy=self.config.ploidy,
            convert_het=self.config.convert_het,
            output_format=self.config.output_format,
        )
        self.vcf_writer.write_snp_sets(
            self.config.input_vcf, self.config.output_prefix, snp_sets, write_config
        )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Done")

        # Write summary
        summary_path = self.config.output_prefix / "summary.tsv"
        self.summary_writer.write_summary(
            snp_sets, self.config.output_prefix, summary_path, vcf_data
        )
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Summary TSV written: {summary_path}")

        # Write run information
        config_data = {
            "input_vcf": self.config.input_vcf,
            "input_format": self.config.input_format,
            "output_prefix": self.config.output_prefix,
            "output_format": self.config.output_format,
            "ploidy": self.config.ploidy,
            "min_distance": self.config.min_distance,
            "convert_het": self.config.convert_het,
            "n_sets": self.config.n_sets,
            "weight_entropy": self.config.weight_entropy,
            "weight_maf": self.config.weight_maf,
            "verbose": self.config.verbose,
            "log_level": self.config.log_level,
            "log_format": self.config.log_format,
        }
        run_info_path = self.run_info_writer.write_run_info(
            self.config.output_prefix, vcf_data, snp_sets, config_data, start_time
        )
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Run information written to: {run_info_path}")

        # Write output directory structure information
        structure_info_path = self.structure_info_writer.write_structure_info(
            self.config.output_prefix, snp_sets
        )
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Output structure information written to: {structure_info_path}"
            )

        # Find best SNP set and copy to best_set.vcf
        if snp_sets:
            best_set_index, best_snp_set, best_entropy = (
                self.summary_writer.find_best_snp_set(snp_sets, vcf_data)
            )

            if self.logger.isEnabledFor(logging.INFO):
                self.logger.info(
                    f"Best SNP set identified: set #{best_set_index} "
                    f"with Shannon entropy {best_entropy:.3f}"
                )

            # Copy best set to best_set.vcf
            self.vcf_writer.copy_snp_set_to_best_set(
                self.config.input_vcf,
                self.config.output_prefix,
                best_snp_set,
                write_config,
            )

            if self.logger.isEnabledFor(logging.INFO):
                best_set_path = self.config.output_prefix / "best_set.vcf"
                self.logger.info(f"Best SNP set copied to: {best_set_path}")

        # Final memory usage summary
        self.memory_monitor.check_memory_and_warn("processing complete")
        if self.logger.isEnabledFor(logging.INFO):
            final_memory = self.memory_monitor.get_memory_usage_mb()
            self.logger.info(f"Final memory usage: {final_memory:.1f}MB")
