"""TSV summary file output operations."""

from pathlib import Path
from typing import List, Union
from dataclasses import dataclass

from ..core.distance_calculator import DistanceCalculator
from ..core.vcf_parser import VCFData, LazyVCFData

__all__ = ["SetStatistics", "SummaryWriter"]


@dataclass
class SetStatistics:
    """Statistics for a SNP set.

    Attributes:
        set_index: 1-based index of the SNP set
        output_vcf: Name of the output VCF file
        num_snps: Number of SNPs in the set
        min_distance: Minimum Hamming distance between samples
        snp_ids: Comma-separated list of SNP IDs

    Example:
        >>> stats = SetStatistics(
        ...     set_index=1, output_vcf="output_1.vcf",
        ...     num_snps=5, min_distance=3.0, snp_ids="rs1,rs2,rs3,rs4,rs5"
        ... )
        >>> print(f"Set {stats.set_index}: {stats.num_snps} SNPs, distance {stats.min_distance}")
        Set 1: 5 SNPs, distance 3.0
    """

    set_index: int
    output_vcf: str
    num_snps: int
    min_distance: float
    snp_ids: str


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

    def write_summary(
        self,
        snp_sets: List[List[str]],
        output_prefix: Path,
        output_path: Path,
        vcf_data: Union[VCFData, LazyVCFData],
    ) -> None:
        """Write per-set summary statistics to TSV.

        Args:
            snp_sets: List of SNP sets, each containing SNP IDs
            output_prefix: Prefix for output VCF files
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
            >>> # set_index	output_vcf	num_snps	min_distance	snp_ids
            >>> # 1	output_1.vcf	2	3	rs1,rs2
            >>> # 2	output_2.vcf	2	4	rs3,rs4
        """
        header = [
            "set_index",
            "output_vcf",
            "num_snps",
            "min_distance",
            "snp_ids",
        ]

        lines = ["\t".join(header)]

        for si, s in enumerate(snp_sets, start=1):
            out_vcf = f"{output_prefix}_{si}.vcf"
            min_d = self.distance_calc.calc_distance_for_snp_ids(
                s, vcf_data.snp_genotypes
            )
            snp_ids_str = ",".join(sorted(s))
            row = [str(si), out_vcf, str(len(s)), str(int(min_d)), snp_ids_str]
            lines.append("\t".join(row))

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
