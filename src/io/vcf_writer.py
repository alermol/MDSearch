"""VCF file output operations."""

import logging
from pathlib import Path
from typing import List, Any
from dataclasses import dataclass

from ..core.genotype_utils import is_het
import pysam

__all__ = ["WriteConfig", "VCFWriter"]


@dataclass
class WriteConfig:
    """Configuration for VCF writing.

    Attributes:
        ploidy: Ploidy level (e.g., 2 for diploid)
        convert_het: Whether to convert heterozygous calls to missing
        output_format: Output format identifier (v, z, u, b)

    Example:
        >>> config = WriteConfig(ploidy=2, convert_het=True, output_format="z")
        >>> print(f"Output: {config.output_format}, Convert het: {config.convert_het}")
        Output: z, Convert het: True
    """

    ploidy: int
    convert_het: bool
    output_format: str = "v"  # one of: v|z|u|b (bcftools letters)


class VCFWriter:
    """Handles VCF file output operations."""

    def __init__(self, logger: logging.Logger = None):
        """Initialize VCF writer with optional logger."""
        self.logger = logger

    def write_snp_sets(
        self,
        input_vcf: Path,
        output_prefix: Path,
        snp_sets: List[List[str]],
        config: WriteConfig,
    ) -> None:
        """Write each SNP set to output files.

        - Uses pysam VariantFile for all formats (vcf, vcf.gz, bcf).
        - Creates output folder structure with 'mdss' subdirectory.

        Args:
            input_vcf: Path to input VCF file
            output_prefix: Path to output folder (will be created if absent)
            snp_sets: List of SNP sets, each containing SNP IDs
            config: Write configuration specifying format and options

        Example:
            >>> from pathlib import Path
            >>> from src.io.vcf_writer import VCFWriter, WriteConfig
            >>>
            >>> writer = VCFWriter()
            >>> config = WriteConfig(ploidy=2, convert_het=False, output_format="v")
            >>> snp_sets = [["rs1", "rs2"], ["rs3", "rs4"]]
            >>>
            >>> writer.write_snp_sets(
            ...     Path("input.vcf"), Path("output"), snp_sets, config
            ... )
            >>> # Creates output/mdss/output_1.vcf and output/mdss/output_2.vcf
        """
        if self.logger and self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Writing {len(snp_sets)} SNP set(s) to output directory: {output_prefix}"
            )
            self.logger.info(
                f"Output format: {config.output_format}, Ploidy: {config.ploidy}"
            )
            if config.convert_het:
                self.logger.info("Converting heterozygous calls to missing values")

        # Create output directory structure
        mdss_dir = output_prefix / "mdss"
        mdss_dir.mkdir(parents=True, exist_ok=True)

        if self.logger and self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug(f"Created output directory: {mdss_dir}")

        # Map bcftools letters to pysam modes
        mode_map = {"v": "w", "z": "wz", "u": "wb0", "b": "wb"}
        if config.output_format not in mode_map:
            raise ValueError(f"Unsupported output_format: {config.output_format}")
        out_mode: Any = mode_map[config.output_format]

        with pysam.VariantFile(str(input_vcf)) as invcf:
            for si, s in enumerate(snp_sets, start=1):
                if self.logger and self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(f"Processing SNP set #{si} with {len(s)} SNPs")

                # Gather contigs used by this set
                contigs_needed = set()
                for rec in invcf.fetch():
                    if rec.id and rec.id in s:
                        contigs_needed.add(rec.chrom)

                # Prepare header: copy and add missing contigs and GT format
                header = invcf.header.copy()
                for chrom in sorted(contigs_needed):
                    if chrom not in header.contigs:
                        header.contigs.add(chrom)
                if "GT" not in header.formats:
                    header.formats.add("GT", 1, "String", "Genotype")

                if config.output_format == "v":
                    output_file = mdss_dir / f"minimal_set_{si}.vcf"
                    if str(input_vcf).endswith(".vcf"):
                        with open(input_vcf) as infh, open(output_file, "w") as outfh:
                            for vcf_line in infh:
                                if vcf_line.startswith("#"):
                                    outfh.write(vcf_line)
                                else:
                                    parts = vcf_line.strip().split("\t")
                                    if len(parts) > 2 and parts[2] in s:
                                        outfh.write(vcf_line)
                    else:
                        with pysam.VariantFile(
                            str(output_file), mode=out_mode, header=header
                        ) as outvcf:
                            for rec in invcf.fetch():
                                if rec.id and rec.id in s:
                                    outvcf.write(rec)
                else:
                    if config.output_format == "z":
                        output_file = mdss_dir / f"minimal_set_{si}.vcf.gz"
                    elif config.output_format == "u":
                        output_file = mdss_dir / f"minimal_set_{si}.bcf"
                    elif config.output_format == "b":
                        output_file = mdss_dir / f"minimal_set_{si}.bcf"
                    else:
                        output_file = (
                            mdss_dir / f"minimal_set_{si}.{config.output_format}"
                        )

                    with pysam.VariantFile(
                        str(output_file), mode=out_mode, header=header
                    ) as outvcf:
                        for rec in invcf.fetch():
                            if rec.id and rec.id in s:
                                outvcf.write(rec)

                if self.logger and self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(f"Written SNP set #{si} to: {output_file}")

        if self.logger and self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Successfully wrote all {len(snp_sets)} SNP sets to {mdss_dir}"
            )

    def _convert_het_in_record(
        self, rec: pysam.VariantRecord, config: WriteConfig
    ) -> None:
        """Convert heterozygous GTs to missing for a record in-place.

        Args:
            rec: pysam VariantRecord to modify
            config: Write configuration specifying ploidy

        Example:
            >>> # This method is called internally by write_snp_sets
            >>> # when convert_het=True in the config
        """
        if config.ploidy and config.ploidy > 1:
            missing_tuple = tuple([None] * config.ploidy)
        else:
            missing_tuple = (None,)

        for sample in rec.header.samples:
            sdata = rec.samples[sample]
            gt = sdata.get("GT")
            if gt is None:
                continue
            sep_detect = "|" if any("|" in str(x) for x in (sdata.get("GT"),)) else "/"
            gt_str = sep_detect.join(["." if g is None else str(g) for g in gt])
            if is_het(gt_str):
                sdata["GT"] = missing_tuple
                is_phased = "|" in gt_str
                if is_phased and len(missing_tuple) > 1:
                    try:
                        sdata.phased = True
                    except Exception:
                        pass

    def _write_line_with_het_conversion(
        self, line: List[str], outfh: Any, config: WriteConfig
    ) -> None:
        """Write VCF line with heterozygous calls converted to missing (text path).

        Args:
            line: VCF line parts as list of strings
            outfh: Output file handle
            config: Write configuration specifying format and options

        Example:
            >>> # This method is called internally by write_snp_sets
            >>> # when processing text VCF files with convert_het=True
        """
        format_field = line[8] if len(line) > 8 else "GT"
        keys = format_field.split(":") if format_field else []
        gt_index = keys.index("GT") if "GT" in keys else -1

        first_gt = None
        if gt_index >= 0 and len(line) > 9:
            first_parts = line[9].split(":")
            if gt_index < len(first_parts):
                first_gt = first_parts[gt_index]
        sep = "/" if (first_gt and "/" in first_gt) else "|"
        missing_gt = (
            sep.join(["."] * config.ploidy)
            if (config.ploidy and config.ploidy > 1)
            else "."
        )

        new_samples: List[str] = []
        for sample_field in line[9:]:
            if gt_index == -1:
                new_samples.append(sample_field)
                continue
            parts = sample_field.split(":")
            gt_val = parts[gt_index] if gt_index < len(parts) else "."
            if is_het(gt_val):
                if gt_index < len(parts):
                    parts[gt_index] = missing_gt
                    new_samples.append(":".join(parts))
                else:
                    new_samples.append(sample_field)
            else:
                new_samples.append(sample_field)

        line = line[:9] + new_samples
        outfh.write("\t".join(line) + "\n")

    def copy_snp_set_to_best_set(
        self,
        input_vcf: Path,
        output_prefix: Path,
        snp_set: List[str],
        config: WriteConfig,
    ) -> None:
        """Copy a specific SNP set to best_set.vcf in the output directory.

        Args:
            input_vcf: Path to input VCF file
            output_prefix: Path to output folder (will be created if absent)
            snp_set: List of SNP IDs to include in the best set
            config: Write configuration specifying format and options

        Example:
            >>> from pathlib import Path
            >>> from src.io.vcf_writer import VCFWriter, WriteConfig
            >>>
            >>> writer = VCFWriter()
            >>> config = WriteConfig(ploidy=2, convert_het=False, output_format="v")
            >>> snp_set = ["rs1", "rs2", "rs3"]
            >>>
            >>> writer.copy_snp_set_to_best_set(
            ...     Path("input.vcf"), Path("output"), snp_set, config
            ... )
            >>> # Creates output/best_set.vcf with the specified SNPs
        """
        # Create output directory if it doesn't exist
        output_prefix.mkdir(parents=True, exist_ok=True)

        # Map bcftools letters to pysam modes
        mode_map = {"v": "w", "z": "wz", "u": "wb0", "b": "wb"}
        if config.output_format not in mode_map:
            raise ValueError(f"Unsupported output_format: {config.output_format}")
        out_mode: Any = mode_map[config.output_format]

        with pysam.VariantFile(str(input_vcf)) as invcf:
            # Gather contigs used by this set
            contigs_needed = set()
            for rec in invcf.fetch():
                if rec.id and rec.id in snp_set:
                    contigs_needed.add(rec.chrom)

            # Prepare header: copy and add missing contigs and GT format
            header = invcf.header.copy()
            for chrom in sorted(contigs_needed):
                if chrom not in header.contigs:
                    header.contigs.add(chrom)
            if "GT" not in header.formats:
                header.formats.add("GT", 1, "String", "Genotype")

            # Output file name based on format
            if config.output_format == "v":
                output_file = output_prefix / "best_set.vcf"
                if str(input_vcf).endswith(".vcf"):
                    with open(input_vcf) as infh, open(output_file, "w") as outfh:
                        for vcf_line in infh:
                            if vcf_line.startswith("#"):
                                outfh.write(vcf_line)
                            else:
                                parts = vcf_line.strip().split("\t")
                                if len(parts) > 2 and parts[2] in snp_set:
                                    if config.convert_het:
                                        self._write_line_with_het_conversion(
                                            parts, outfh, config
                                        )
                                    else:
                                        outfh.write(vcf_line)
                else:
                    with open(output_file, "w") as outfh:
                        for line in str(header).split("\n"):
                            if line.strip():
                                outfh.write(line + "\n")

                        for rec in invcf.fetch():
                            if not rec.id or rec.id not in snp_set:
                                continue
                            if config.convert_het:
                                self._convert_het_in_record(rec, config)
                            outfh.write(str(rec))
            else:
                if config.output_format == "z":
                    output_file = output_prefix / "best_set.vcf.gz"
                elif config.output_format == "u":
                    output_file = output_prefix / "best_set.vcf.bz2"
                elif config.output_format == "b":
                    output_file = output_prefix / "best_set.bcf"
                else:
                    output_file = output_prefix / "best_set.bcf"

                with pysam.VariantFile(
                    str(output_file), out_mode, header=header
                ) as outvcf:
                    for rec in invcf.fetch():
                        if not rec.id or rec.id not in snp_set:
                            continue
                        if config.convert_het:
                            self._convert_het_in_record(rec, config)
                        outvcf.write(rec)
