"""Main application coordinator for MDSearch."""

import logging
from pathlib import Path
from typing import Optional
from dataclasses import dataclass, field

from .core import VCFParser, DistanceCalculator, SNPSelector
from .core.snp_selector import OverlapConstraints
from .io import VCFWriter, WriteConfig, SummaryWriter
from .utils import MemoryMonitor, setup_logger


@dataclass
class MDSearchConfig:
    """Configuration for MDSearch application."""

    input_vcf: Path
    output_prefix: Path
    ploidy: int = 2
    max_snps: int = 0
    min_distance: int = 1
    convert_het: bool = False
    n_sets: int = 1
    overlap_constraints: OverlapConstraints = field(default_factory=OverlapConstraints)
    verbose: bool = True
    log_level: Optional[str] = None
    log_format: str = "text"
    summary_tsv: Optional[Path] = None


class MDSearchApp:
    """Main application coordinator with separated concerns."""

    def __init__(self, config: MDSearchConfig):
        self.config = config
        self.logger = setup_logger(
            "mdsearch", config.log_level, config.log_format, config.verbose
        )
        self.memory_monitor = MemoryMonitor(self.logger)

        # Initialize components
        self.vcf_parser = VCFParser(self.memory_monitor, self.logger)
        self.distance_calc = DistanceCalculator(self.memory_monitor, self.logger)
        self.snp_selector = SNPSelector(
            self.distance_calc, self.memory_monitor, self.logger
        )
        self.vcf_writer = VCFWriter()
        self.summary_writer = SummaryWriter(self.distance_calc)

        # Initial memory check
        self.memory_monitor.check_memory_and_warn("initialization")

    def run(self) -> None:
        """Execute the complete MDSearch pipeline."""
        # 1. Parse and validate VCF
        vcf_data = self.vcf_parser.parse_and_validate(
            self.config.input_vcf, self.config.ploidy, self.config.convert_het
        )

        # 2. Search for optimal SNP sets
        snp_sets = self.snp_selector.search_optimal_sets(
            vcf_data,
            self.config.min_distance,
            self.config.n_sets,
            self.config.max_snps,
            self.config.overlap_constraints,
        )

        # 3. Log selected SNP details when not using TSV summary
        if not self.config.summary_tsv:
            for si, s in enumerate(snp_sets, start=1):
                min_d = self.distance_calc.calc_distance_for_snp_ids(
                    s, vcf_data.snp_genotypes
                )
                snp_ids_str = ",".join(sorted(s))
                self.logger.info(
                    f"Set {si}: {len(s)} SNPs, min_distance={int(min_d)}, "
                    f"SNP_IDs=[{snp_ids_str}]"
                )

        # 4. Write output files
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Writing selected SNPs in VCF...")

        write_config = WriteConfig(
            ploidy=self.config.ploidy, convert_het=self.config.convert_het
        )
        self.vcf_writer.write_snp_sets(
            self.config.input_vcf, self.config.output_prefix, snp_sets, write_config
        )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Done")

        # 5. Write summary if requested
        if self.config.summary_tsv:
            self.summary_writer.write_summary(
                snp_sets, self.config.output_prefix, self.config.summary_tsv, vcf_data
            )
            if self.logger.isEnabledFor(logging.INFO):
                self.logger.info(f"Summary TSV written: {self.config.summary_tsv}")

        # Final memory usage summary
        self.memory_monitor.check_memory_and_warn("processing complete")
        if self.logger.isEnabledFor(logging.INFO):
            final_memory = self.memory_monitor.get_memory_usage_mb()
            self.logger.info(f"Final memory usage: {final_memory:.1f}MB")
