"""VCF file output operations."""

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

    def write_snp_sets(
        self,
        input_vcf: Path,
        output_prefix: Path,
        snp_sets: List[List[str]],
        config: WriteConfig,
    ) -> None:
        """Write each SNP set to output files.

        - Uses pysam VariantFile for all formats (vcf, vcf.gz, bcf).

        Args:
            input_vcf: Path to input VCF file
            output_prefix: Prefix for output file names
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
            >>> # Creates output_1.vcf and output_2.vcf
        """
        # Map bcftools letters to pysam modes
        mode_map = {"v": "w", "z": "wz", "u": "wb0", "b": "wb"}
        if config.output_format not in mode_map:
            raise ValueError(f"Unsupported output_format: {config.output_format}")
        out_mode: Any = mode_map[config.output_format]

        with pysam.VariantFile(str(input_vcf)) as invcf:
            for si, s in enumerate(snp_sets, start=1):
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

                # Output file name based on format
                if config.output_format == "v":
                    output_file = f"{output_prefix}_{si}.vcf"
                    # Text VCF path: preserve original FORMAT and sample subfields when input is text.
                    if str(input_vcf).endswith(".vcf"):
                        with open(input_vcf) as infh, open(output_file, "w") as outfh:
                            for vcf_line in infh:
                                if vcf_line.startswith("#"):
                                    outfh.write(vcf_line)
                                else:
                                    parts = vcf_line.rstrip("\n").split("\t")
                                    if len(parts) < 3:
                                        continue
                                    vid = parts[2]
                                    if not vid or vid not in s:
                                        continue
                                    if config.convert_het and len(parts) >= 9:
                                        self._write_line_with_het_conversion(
                                            parts, outfh, config
                                        )
                                    else:
                                        outfh.write(vcf_line)
                    else:
                        # Compressed inputs: synthesize simple GT-only lines using pysam header
                        with open(output_file, "w") as outfh:
                            outfh.write(str(header))
                            for rec in invcf.fetch():
                                if not rec.id or rec.id not in s:
                                    continue
                                chrom = rec.chrom
                                pos = rec.pos
                                vid = rec.id
                                ref = rec.ref
                                alt = rec.alts[0] if rec.alts else "."
                                qual = "."
                                flt = "PASS"
                                info = "."
                                fmt = "GT"
                                sample_strs: List[str] = []
                                for sample in invcf.header.samples:
                                    data = rec.samples[sample]
                                    gt_tuple = data.get("GT")
                                    phased = getattr(data, "phased", False)
                                    sep = "|" if phased else "/"
                                    if gt_tuple is None or all(
                                        g is None for g in gt_tuple
                                    ):
                                        gt_s = sep.join(["."] * max(1, config.ploidy))
                                    else:
                                        parts = [
                                            "." if g is None else str(g)
                                            for g in gt_tuple
                                        ]
                                        gt_s = sep.join(parts)
                                    if config.convert_het and is_het(gt_s):
                                        gt_s = sep.join(["."] * max(1, config.ploidy))
                                    sample_strs.append(gt_s)
                                line = "\t".join(
                                    [
                                        str(chrom),
                                        str(pos),
                                        str(vid),
                                        str(ref),
                                        str(alt),
                                        qual,
                                        flt,
                                        info,
                                        fmt,
                                    ]
                                    + sample_strs
                                )
                                outfh.write(line + "\n")
                else:
                    if config.output_format == "z":
                        output_file = f"{output_prefix}_{si}.vcf.gz"
                    elif config.output_format == "u":
                        output_file = f"{output_prefix}_{si}.bcf"
                    else:
                        output_file = f"{output_prefix}_{si}.bcf"
                    with pysam.VariantFile(
                        output_file, out_mode, header=header
                    ) as outvcf:
                        for rec in invcf.fetch():
                            if not rec.id or rec.id not in s:
                                continue
                            if config.convert_het:
                                self._convert_het_in_record(rec, config)
                            outvcf.write(rec)

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
        # Compose missing GT according to ploidy
        if config.ploidy and config.ploidy > 1:
            missing_tuple = tuple([None] * config.ploidy)
        else:
            missing_tuple = (None,)

        for sample in rec.header.samples:
            sdata = rec.samples[sample]
            gt = sdata.get("GT")
            if gt is None:
                continue
            # Determine heterozygosity by string form for reliability
            sep_detect = "|" if any("|" in str(x) for x in (sdata.get("GT"),)) else "/"
            gt_str = sep_detect.join(["." if g is None else str(g) for g in gt])
            if is_het(gt_str):
                sdata["GT"] = missing_tuple
                # Preserve phasing for missing if original looked phased
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
            config: Write configuration specifying ploidy

        Example:
            >>> # This method is called internally by write_snp_sets
            >>> # when processing text VCF files with convert_het=True
        """
        format_field = line[8] if len(line) > 8 else "GT"
        keys = format_field.split(":") if format_field else []
        gt_index = keys.index("GT") if "GT" in keys else -1

        # Determine separator from first sample's GT
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
