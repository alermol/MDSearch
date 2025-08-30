"""TSV summary file output operations."""

from pathlib import Path
from typing import List, Tuple
from dataclasses import dataclass
import logging

from ..core.distance_calculator import DistanceCalculator
from ..core.vcf_parser import VCFData
from ..core.genotype_utils import calculate_shannon_entropy
from collections import Counter

__all__ = ["SetStatistics", "SummaryWriter"]


@dataclass
class SetStatistics:
    """Statistics for a SNP set.

    Attributes:
        set_index: 1-based index of the SNP set
        output_vcf: Name of the output VCF file
        num_snps: Number of SNPs in the set
        min_distance: Minimum Hamming distance between samples
        shannon_entropy: Shannon entropy of chromosome distribution

    Example:
        >>> stats = SetStatistics(
        ...     set_index=1, output_vcf="output_1.vcf",
        ...     num_snps=5, min_distance=3.0, shannon_entropy=1.585
        ... )
        >>> print(f"Set {stats.set_index}: {stats.num_snps} SNPs, distance {stats.min_distance}, entropy {stats.shannon_entropy}")
        Set 1: 5 SNPs, distance 3.0, entropy 1.585
    """

    set_index: int
    output_vcf: str
    num_snps: int
    min_distance: float
    shannon_entropy: float


class SummaryWriter:
    """Handles TSV summary file output."""

    def __init__(self, distance_calc: DistanceCalculator):
        """Initialize summary writer with distance calculator.

        Args:
            distance_calc: DistanceCalculator instance for computing distances

        Example:
            >>> from src.core.distance_calculator import DistanceCalculator
            >>> from src.utils.memory_monitor import MemoryMonitor
            >>> from src.utils.logging_setup import setup_logger
            >>> logger = setup_logger("summary_writer")
            >>> memory_monitor = MemoryMonitor(logger)
            >>> distance_calc = DistanceCalculator(memory_monitor, logger)
            >>> writer = SummaryWriter(distance_calc)
        """
        self.distance_calc = distance_calc
        self.logger = logging.getLogger(__name__)

    def find_best_snp_set(
        self,
        snp_sets: List[List[str]],
        vcf_data: VCFData,
    ) -> Tuple[int, List[str], float]:
        """Find the best SNP set based on Shannon entropy.

        The best set is determined by the highest Shannon entropy value,
        which indicates the most balanced distribution of SNPs across chromosomes.

        Args:
            snp_sets: List of SNP sets, each containing SNP IDs
            vcf_data: VCF data containing SNP information

        Returns:
            Tuple of (best_set_index, best_snp_set, best_entropy)
            - best_set_index: 1-based index of the best set
            - best_snp_set: List of SNP IDs in the best set
            - best_entropy: Shannon entropy value of the best set

        Example:
            >>> from src.io.summary_writer import SummaryWriter
            >>> writer = SummaryWriter(distance_calc)
            >>> snp_sets = [["rs1", "rs2"], ["rs3", "rs4"]]
            >>> best_idx, best_set, best_entropy = writer.find_best_snp_set(snp_sets, vcf_data)
            >>> print(f"Best set: {best_idx}, Entropy: {best_entropy:.3f}")
            Best set: 2, Entropy: 1.585
        """
        if not snp_sets:
            raise ValueError("No SNP sets provided")

        best_set_index = 1
        best_snp_set = snp_sets[0]
        best_entropy = 0.0

        for si, s in enumerate(snp_sets, start=1):
            chromosome_counts: Counter[str] = Counter()

            for snp_id in s:
                if snp_id in vcf_data.snp_genotypes:
                    chromosome_counts[vcf_data.snp_genotypes[snp_id].chromosome] += 1

            for chrom in vcf_data.headers.contigs:
                if chrom not in chromosome_counts:
                    chromosome_counts[chrom] = 0

            current_entropy = calculate_shannon_entropy(chromosome_counts)

            if current_entropy > best_entropy:
                best_entropy = current_entropy
                best_snp_set = s
                best_set_index = si

        return best_set_index, best_snp_set, best_entropy

    def write_summary(
        self,
        snp_sets: List[List[str]],
        output_prefix: Path,
        output_path: Path,
        vcf_data: VCFData,
    ) -> None:
        """Write per-set summary statistics to TSV.

        Args:
            snp_sets: List of SNP sets, each containing SNP IDs
            output_prefix: Path to output folder containing VCF files in 'mdss' subdirectory
            output_path: Path for the summary TSV file
            vcf_data: VCF data containing SNP information

        Example:
            >>> from pathlib import Path
            >>> from src.io.summary_writer import SummaryWriter
            >>>
            >>> writer = SummaryWriter(distance_calc)
            >>> snp_sets = [["rs1", "rs2"], ["rs3", "rs4"]]
            >>>
            >>> writer.write_summary(
            ...     snp_sets, Path("output"), Path("summary.tsv"), vcf_data
            ... )
            >>> # Creates summary.tsv with statistics for each set
            >>>
            >>> # Example output format:
            >>> # set_index	output_vcf	num_snps	min_distance	shannon_entropy
            >>> # 1	output_1.vcf	2	3	1.000
            >>> # 2	output_2.vcf	2	4	1.585
        """
        if self.logger and self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Generating summary statistics for {len(snp_sets)} SNP set(s)"
            )
            self.logger.info(f"Summary will be written to: {output_path}")

        header = [
            "set_index",
            "output_vcf",
            "num_snps",
            "min_distance",
            "shannon_entropy",
        ]

        lines = ["\t".join(header)]

        for si, s in enumerate(snp_sets, start=1):
            if self.logger and self.logger.isEnabledFor(logging.DEBUG):
                self.logger.debug(f"Processing SNP set #{si} with {len(s)} SNPs")

            out_vcf = f"minimal_set_{si}.vcf"
            min_d = self.distance_calc.calc_distance_for_snp_ids(
                s, vcf_data.snp_genotypes
            )

            chromosome_counts: Counter[str] = Counter()

            for snp_id in s:
                if snp_id in vcf_data.snp_genotypes:
                    chromosome_counts[vcf_data.snp_genotypes[snp_id].chromosome] += 1

            for chrom in vcf_data.headers.contigs:
                if chrom not in chromosome_counts:
                    chromosome_counts[chrom] = 0

            shannon_entropy = calculate_shannon_entropy(chromosome_counts)

            if self.logger and self.logger.isEnabledFor(logging.DEBUG):
                self.logger.debug(
                    f"Set #{si}: distance={min_d}, entropy={shannon_entropy:.3f}, "
                    f"chromosomes={dict(chromosome_counts)}"
                )

            row = [
                str(si),
                out_vcf,
                str(len(s)),
                str(int(min_d)),
                f"{shannon_entropy:.3f}",
            ]
            lines.append("\t".join(row))

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

        if self.logger and self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Summary statistics written successfully to: {output_path}"
            )
