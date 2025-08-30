"""Main application coordinator for MDSearch."""

import logging
from pathlib import Path
from typing import Optional, Union, Callable
from dataclasses import dataclass

from .core import VCFParser, DistanceCalculator, SNPSelector
from .io import VCFWriter, WriteConfig, SummaryWriter
from .utils import MemoryMonitor, setup_logger
from .utils import ensure_variant_index

__all__ = ["MDSearchConfig", "MDSearchApp"]


@dataclass
class MDSearchConfig:
    """Configuration for MDSearch application.

    Attributes:
        input_vcf: Path to input VCF file
        output_prefix: Path to output folder (will be created if absent)
        ploidy: Ploidy level (e.g., 2 for diploid)
        min_distance: Minimum Hamming distance required
        convert_het: Whether to convert heterozygous calls to missing
        n_sets: Number of alternative SNP sets to generate (0 for unlimited)
        verbose: Whether to enable verbose logging
        log_level: Logging level override
        log_format: Logging format (text or json)
        input_format: Input format identifier (auto, v, z, u, b)
        output_format: Output format identifier (v, z, u, b)

    Example:
        >>> from pathlib import Path
        >>> config = MDSearchConfig(
        ...     input_vcf=Path("sample.vcf"),
        ...     output_prefix=Path("output"),
        ...     ploidy=2,
        ...     min_distance=3,
        ...     n_sets=2,
        ...     # lazy_loading removed
        ... )
        >>> print(f"Input: {config.input_vcf}, Output folder: {config.output_prefix}")
        Input: sample.vcf, Output folder: output
    """

    input_vcf: Path
    output_prefix: Path
    ploidy: int = 2
    min_distance: int = 1
    convert_het: bool = False
    n_sets: Union[int, str] = 1
    verbose: bool = True
    log_level: Optional[str] = None
    log_format: str = "text"
    # Lazy loading fields removed - no longer needed
    input_format: str = "auto"  # one of: auto|v|z|u|b (bcftools letters)
    output_format: str = "v"  # one of: v|z|u|b (bcftools letters)


class MDSearchApp:
    """Main application coordinator with separated concerns."""

    def __init__(
        self,
        config: MDSearchConfig,
        shutdown_checker: Optional[Callable[[], bool]] = None,
    ):
        """Initialize MDSearch application with configuration.

        Args:
            config: Application configuration
            shutdown_checker: Optional function to check if shutdown was requested

        Example:
            >>> from src.app import MDSearchApp, MDSearchConfig
            >>> from pathlib import Path
            >>>
            >>> config = MDSearchConfig(
            ...     input_vcf=Path("input.vcf"),
            ...     output_prefix=Path("output"),
            ...     ploidy=2,
            ...     min_distance=3
            ... )
            >>> app = MDSearchApp(config)
            >>> print(f"Initialized with {config.ploidy}-ploidy analysis")
            Initialized with 2-ploidy analysis
        """
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
            self.distance_calc, self.memory_monitor, self.logger, shutdown_checker
        )
        self.vcf_writer = VCFWriter()
        self.summary_writer = SummaryWriter(self.distance_calc)

        # Initial memory check
        self.memory_monitor.check_memory_and_warn("initialization")

    def run(self) -> None:
        """Execute the complete MDSearch pipeline.

        This method orchestrates the entire workflow:
        1. Index verification for compressed inputs
        2. VCF parsing and validation
        3. SNP set selection and optimization
        4. Output file generation
        5. Summary statistics (if requested)

        Example:
            >>> app = MDSearchApp(config)
            >>> app.run()
            >>> # Will process VCF, select SNPs, and generate output files
        """
        # 1. Ensure index exists for compressed inputs (VCF.gz/BCF)
        try:
            ensure_variant_index(
                self.config.input_vcf, self.config.input_format, self.logger
            )
        except Exception as e:
            if self.logger.isEnabledFor(logging.ERROR):
                self.logger.error(f"Failed to ensure variant index: {e}")
            raise

        # 2. Parse and validate VCF
        vcf_data = self.vcf_parser.parse_and_validate(
            self.config.input_vcf, self.config.ploidy, self.config.convert_het
        )

        # 3. Search for optimal SNP sets
        snp_sets = self.snp_selector.search_optimal_sets(
            vcf_data,
            self.config.min_distance,
            self.config.n_sets,
        )

        # 5. Write output files
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

        # 6. Write summary (always created)
        summary_path = self.config.output_prefix / "summary.tsv"
        self.summary_writer.write_summary(
            snp_sets, self.config.output_prefix, summary_path, vcf_data
        )
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Summary TSV written: {summary_path}")

        # 7. Find best SNP set based on Shannon entropy and copy to best_set.vcf
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

        # Final memory usage summary and cache statistics
        self.memory_monitor.check_memory_and_warn("processing complete")
        if self.logger.isEnabledFor(logging.INFO):
            final_memory = self.memory_monitor.get_memory_usage_mb()
            self.logger.info(f"Final memory usage: {final_memory:.1f}MB")


# Lazy loading cache statistics removed - no longer needed
