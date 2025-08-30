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
        self.vcf_writer = VCFWriter(self.logger)
        self.summary_writer = SummaryWriter(self.distance_calc)
        self.run_info_writer = RunInfoWriter(self.memory_monitor)
        self.structure_info_writer = StructureInfoWriter(self.config.output_format)

        self.memory_monitor.check_memory_and_warn("initialization")

    def run(self) -> None:
        """Execute the complete MDSearch pipeline.
        
        Performs the full SNP selection workflow:
        1. Index verification
        2. VCF parsing and validation
        3. SNP set generation
        4. Output file creation
        """
        start_time = time.time()

        try:
            index_start = time.time()
            ensure_variant_index(
                self.config.input_vcf, self.config.input_format, self.logger
            )
            index_time = time.time() - index_start
            if self.logger.isEnabledFor(logging.DEBUG):
                self.logger.debug(f"Index verification completed in {index_time:.2f}s")
        except Exception as e:
            if self.logger.isEnabledFor(logging.ERROR):
                self.logger.error(f"Failed to ensure variant index: {e}")
            raise

        parse_start = time.time()
        vcf_data = self.vcf_parser.parse_and_validate(
            self.config.input_vcf, self.config.ploidy, self.config.convert_het
        )
        parse_time = time.time() - parse_start
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"VCF parsing completed in {parse_time:.2f}s")

        search_start = time.time()
        snp_sets = self.snp_selector.search_optimal_sets(
            vcf_data,
            self.config.min_distance,
            self.config.n_sets,
        )
        search_time = time.time() - search_start
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"SNP set search completed in {search_time:.2f}s")

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Writing selected SNPs in VCF...")

        write_start = time.time()
        write_config = WriteConfig(
            ploidy=self.config.ploidy,
            convert_het=self.config.convert_het,
            output_format=self.config.output_format,
        )
        self.vcf_writer.write_snp_sets(
            self.config.input_vcf, self.config.output_prefix, snp_sets, write_config
        )
        write_time = time.time() - write_start
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"VCF writing completed in {write_time:.2f}s")

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Done")

        summary_start = time.time()
        summary_path = self.config.output_prefix / "summary.tsv"
        self.summary_writer.write_summary(
            snp_sets, self.config.output_prefix, summary_path, vcf_data
        )
        summary_time = time.time() - summary_start
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Summary TSV written: {summary_path}")
        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug(f"Summary generation completed in {summary_time:.2f}s")

        run_info_start = time.time()
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
        run_info_time = time.time() - run_info_start
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Run information written to: {run_info_path}")
        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug(f"Run info writing completed in {run_info_time:.2f}s")

        structure_start = time.time()
        structure_info_path = self.structure_info_writer.write_structure_info(
            self.config.output_prefix, snp_sets
        )
        structure_time = time.time() - structure_start
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Output structure information written to: {structure_info_path}"
            )
        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug(
                f"Structure info writing completed in {structure_time:.2f}s"
            )

        if snp_sets:
            best_set_start = time.time()
            best_set_index, best_snp_set, best_entropy = (
                self.summary_writer.find_best_snp_set(snp_sets, vcf_data)
            )
            best_set_time = time.time() - best_set_start

            if self.logger.isEnabledFor(logging.INFO):
                self.logger.info(
                    f"Best SNP set identified: set #{best_set_index} "
                    f"with Shannon entropy {best_entropy:.3f}"
                )

            copy_start = time.time()
            self.vcf_writer.copy_snp_set_to_best_set(
                self.config.input_vcf,
                self.config.output_prefix,
                best_snp_set,
                write_config,
            )
            copy_time = time.time() - copy_start

            if self.logger.isEnabledFor(logging.INFO):
                best_set_path = self.config.output_prefix / "best_set.vcf"
                self.logger.info(f"Best SNP set copied to: {best_set_path}")

            if self.logger.isEnabledFor(logging.DEBUG):
                self.logger.debug(
                    f"Best set identification: {best_set_time:.2f}s, copying: {copy_time:.2f}s"
                )

        self.memory_monitor.check_memory_and_warn("processing complete")
        if self.logger.isEnabledFor(logging.INFO):
            final_memory = self.memory_monitor.get_memory_usage_mb()
            self.logger.info(f"Final memory usage: {final_memory:.1f}MB")

        total_time = time.time() - start_time
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Total execution time: {total_time:.2f}s")

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug("Performance breakdown:")
            self.logger.debug(f"  Index verification: {index_time:.2f}s")
            self.logger.debug(f"  VCF parsing: {parse_time:.2f}s")
            self.logger.debug(f"  SNP set search: {search_time:.2f}s")
            self.logger.debug(f"  VCF writing: {write_time:.2f}s")
            self.logger.debug(f"  Summary generation: {summary_time:.2f}s")
            self.logger.debug(f"  Run info writing: {run_info_time:.2f}s")
            self.logger.debug(f"  Structure info: {structure_time:.2f}s")
            if snp_sets:
                self.logger.debug(
                    f"  Best set operations: {best_set_time + copy_time:.2f}s"
                )
            self.logger.debug(f"  Total: {total_time:.2f}s")
