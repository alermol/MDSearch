"""Output directory structure information file operations."""

import datetime
from pathlib import Path
from typing import List

__all__ = ["StructureInfoWriter"]


class StructureInfoWriter:
    """Handles output directory structure information file output."""

    def __init__(self, output_format: str = "v"):
        """Initialize structure info writer with output format.

        Args:
            output_format: Output format identifier (v, z, u, b)

        Example:
            >>> writer = StructureInfoWriter(output_format="z")
            >>> print(f"Writer initialized for format: {writer.output_format}")
            Writer initialized for format: z
        """
        self.output_format = output_format

    def write_structure_info(
        self, output_prefix: Path, snp_sets: List[List[str]]
    ) -> Path:
        """Write information about the output directory structure to a file.

        Args:
            output_prefix: Path to output directory
            snp_sets: List of SNP sets that were generated

        Returns:
            Path to the created structure info file

        Example:
            >>> from pathlib import Path
            >>> writer = StructureInfoWriter(output_format="v")
            >>> snp_sets = [["rs1", "rs2"], ["rs3", "rs4"]]
            >>> path = writer.write_structure_info(Path("output"), snp_sets)
            >>> print(f"Structure info written to: {path}")
            Structure info written to: output/output_structure.txt
        """
        structure_info_path = output_prefix / "output_structure.txt"

        # Get current timestamp
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC")

        # Prepare structure info content
        lines = [
            "MDSearch Output Directory Structure",
            "===================================",
            "",
            f"Generated: {timestamp}",
            "",
            "Directory Overview:",
            f"  Output root: {output_prefix}",
            f"  Total SNP sets generated: {len(snp_sets)}",
            "",
            "File Structure:",
            "",
            "mdss/ (SNP Set VCF Files)",
            "  This subdirectory contains individual VCF files for each SNP set:",
        ]

        # Add details for each SNP set
        for i, snp_set in enumerate(snp_sets, 1):
            file_extension = self._get_vcf_extension()
            filename = f"minimal_set_{i}.{file_extension}"
            lines.append(f"  • {filename} - {len(snp_set)} SNPs")

        lines.extend(
            [
                "",
                "Root Directory Files:",
                "  • summary.tsv - Summary statistics for all SNP sets",
                "  • best_set.vcf - Copy of the best SNP set (highest Shannon entropy)",
                "  • run_info.txt - Complete run information and memory usage",
                "  • output_structure.txt - This file (directory structure information)",
                "",
                "File Descriptions:",
                "",
                "summary.tsv:",
                "  Tab-separated values file with columns:",
                "  - set_index: Sequential number of the SNP set (1, 2, 3...)",
                "  - output_vcf: Filename of the corresponding VCF file",
                "  - num_snps: Number of SNPs in the set",
                "  - min_distance: Minimum Hamming distance achieved between samples",
                "  - shannon_entropy: Shannon entropy of chromosome distribution",
                "",
                "best_set.vcf:",
                "  Contains the SNP set with the highest Shannon entropy, representing the most balanced chromosome distribution among all generated sets.",
                "",
                "run_info.txt:",
                "  Comprehensive information about the analysis run including:",
                "  - Version and configuration details",
                "  - Input data summary",
                "  - System memory usage (current, peak, available, thresholds)",
                "  - Output statistics",
                "",
                "VCF File Contents:",
                "  Each VCF file contains:",
                "  - Only the SNPs from the corresponding set",
                "  - All original samples from the input VCF",
                "  - Preserved VCF format and metadata",
                "  - Genotype data for selected SNPs only",
                "",
                "Usage Notes:",
                "  • All VCF files maintain the same format as the input file",
                "  • SNP sets are mutually exclusive (no overlapping SNPs)",
                "  • Each set provides the minimum number of SNPs needed for the specified Hamming distance requirement",
                "  • The best_set.vcf is automatically selected based on chromosome distribution balance (Shannon entropy)",
            ]
        )

        # Write to file
        with open(structure_info_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

        return structure_info_path

    def _get_vcf_extension(self) -> str:
        """Get the appropriate file extension based on output format.

        Returns:
            File extension string (vcf, vcf.gz, bcf, etc.)

        Example:
            >>> writer = StructureInfoWriter(output_format="z")
            >>> ext = writer._get_vcf_extension()
            >>> print(f"VCF extension: {ext}")
            VCF extension: vcf.gz
        """
        format_map = {"v": "vcf", "z": "vcf.gz", "u": "bcf", "b": "bcf"}
        return format_map.get(self.output_format, "vcf")
