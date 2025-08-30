"""VCF file parsing and validation."""

import sys
import logging
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional, Callable
from dataclasses import dataclass

from ..utils.memory_monitor import MemoryMonitor
from .genotype_utils import (
    gt_to_value,
    calculate_maf,
    calculate_snp_information_entropy,
)
import pysam

__all__ = [
    "VCFHeaders",
    "SNPData",
    "VCFData",
    "VCFParser",
]


@dataclass
class VCFHeaders:
    """VCF header information."""

    samples: List[str]
    fileformat_found: bool
    contigs: List[str]


@dataclass
class SNPData:
    """Data for a single SNP."""

    genotypes: List[float]
    sample_fields: List[str]
    maf: float
    chromosome: str


@dataclass
class VCFData:
    """Complete VCF data structure."""

    headers: VCFHeaders
    snp_genotypes: Dict[str, SNPData]
    snp_maf_cache: Dict[str, float]
    snp_entropy_cache: Dict[str, float]


class VCFParser:
    """Handles VCF file parsing and validation."""

    def __init__(
        self,
        memory_monitor: MemoryMonitor,
        logger: logging.Logger,
        shutdown_checker: Optional[Callable[[], bool]] = None,
    ):
        """Initialize VCF parser with memory monitoring and logging."""
        self.memory_monitor = memory_monitor
        self.logger = logger
        self.shutdown_checker = shutdown_checker

    def parse_and_validate(
        self, vcf_path: Path, ploidy: int, convert_het: bool
    ) -> VCFData:
        """Parse VCF file and return structured data (loads all SNPs into memory).
        
        Args:
            vcf_path: Path to VCF file to parse
            ploidy: Ploidy level for genotype interpretation
            convert_het: Whether to convert heterozygous calls to missing values
            
        Returns:
            VCFData object containing parsed headers, genotypes, and cached scores
            
        Raises:
            SystemExit: If VCF validation fails (missing headers, malformed data)
        """
        headers = self._validate_headers(vcf_path)

        snp_genotypes, snp_maf_cache, snp_entropy_cache = self._validate_variants(
            vcf_path, headers, ploidy, convert_het
        )

        self.memory_monitor.check_memory_and_warn("VCF parsing")
        self.memory_monitor.warn_for_large_dataset(
            len(snp_genotypes), len(headers.samples)
        )

        return VCFData(
            headers=headers,
            snp_genotypes=snp_genotypes,
            snp_maf_cache=snp_maf_cache,
            snp_entropy_cache=snp_entropy_cache,
        )

    def _validate_headers(self, vcf_path: Path) -> VCFHeaders:
        """Validate VCF headers and extract sample information."""
        try:
            with pysam.VariantFile(str(vcf_path)) as vf:
                samples = list(vf.header.samples)
                fileformat_found = True if getattr(vf.header, "version", None) else True

                contigs = [str(contig) for contig in vf.header.contigs]

                if len(samples) < 2:
                    sys.exit(
                        f"ERROR: Insufficient samples for distance calculation. "
                        f"Found {len(samples)} samples, need at least 2."
                    )

                return VCFHeaders(
                    samples=samples, fileformat_found=fileformat_found, contigs=contigs
                )
        except Exception as e:
            sys.exit(f"ERROR: Failed to read VCF/BCF headers via pysam: {e}")

    def _validate_variants(
        self, vcf_path: Path, headers: VCFHeaders, ploidy: int, convert_het: bool
    ) -> Tuple[Dict[str, SNPData], Dict[str, float], Dict[str, float]]:
        """Parse and validate variant lines."""
        seen_ids: Set[str] = set()
        data_line_count = 0
        snp_genotypes: Dict[str, SNPData] = {}
        snp_maf_cache: Dict[str, float] = {}
        snp_entropy_cache: Dict[str, float] = {}

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Parsing VCF variants with {len(headers.samples)} samples, ploidy={ploidy}"
            )
            if convert_het:
                self.logger.info("Converting heterozygous calls to missing values")

        with pysam.VariantFile(str(vcf_path)) as vf:
            line_number = 0
            for rec in vf:
                line_number += 1
                data_line_count += 1

                if self.shutdown_checker and self.shutdown_checker():
                    if self.logger.isEnabledFor(logging.INFO):
                        self.logger.info(
                            "Graceful shutdown requested. Stopping VCF parsing."
                        )
                    break

                snp_id = rec.id or ""
                if (not snp_id) or (snp_id == "."):
                    sys.exit(
                        f"ERROR: Line {line_number}: Missing or placeholder SNP ID. "
                        "Ensure IDs are present and non-'.'."
                    )

                if snp_id in seen_ids:
                    sys.exit(
                        f"ERROR: Line {line_number}: Duplicate SNP ID '{snp_id}'. "
                        "Ensure SNP IDs are unique."
                    )
                seen_ids.add(snp_id)

                if rec.alts and len(rec.alts) > 1:
                    sys.exit(
                        f"ERROR: Line {line_number}: Multiallelic site detected (ALT='{','.join(rec.alts)}'). "
                        "Filter or split multiallelic sites before processing."
                    )

                if len(rec.samples) != len(headers.samples):
                    sys.exit(
                        f"ERROR: Line {line_number}: Invalid number of samples. "
                        f"Expected {len(headers.samples)}, found {len(rec.samples)}."
                    )

                geno: List[float] = []
                sample_fields_list: List[str] = []
                for sample_idx, sample in enumerate(headers.samples):
                    data = rec.samples[sample]
                    gt_tuple = data.get("GT")
                    if gt_tuple is None or all(g is None for g in gt_tuple):
                        gt = "."
                    else:
                        gt_parts = ["." if g is None else str(g) for g in gt_tuple]
                        gt = "/".join(gt_parts)

                    try:
                        alleles = [
                            int(x)
                            for x in gt.replace("|", "/").split("/")
                            if x not in ("", ".")
                        ]
                    except ValueError:
                        alleles = []
                    if any(a > 1 for a in alleles):
                        sys.exit(
                            f"ERROR: Line {line_number}, Sample {sample_idx + 1}: "
                            f"Multiallelic genotype detected (GT='{gt}'). "
                            "Filter or split multiallelic sites before processing."
                        )

                    sample_fields_list.append(gt)
                    geno.append(gt_to_value(gt, ploidy, convert_het))

                maf = calculate_maf(geno, ploidy)
                entropy = calculate_snp_information_entropy(geno)
                snp_data = SNPData(
                    genotypes=geno,
                    sample_fields=sample_fields_list,
                    maf=maf,
                    chromosome=rec.chrom,
                )
                snp_genotypes[snp_id] = snp_data
                snp_maf_cache[snp_id] = maf
                snp_entropy_cache[snp_id] = entropy

                if data_line_count % 10000 == 0 and self.logger.isEnabledFor(
                    logging.INFO
                ):
                    self.logger.info(f"Processed {data_line_count} variants...")

        if data_line_count == 0:
            sys.exit(
                "ERROR: No data lines found in VCF. "
                "Ensure VCF contains variant records."
            )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"VCF validation complete: {data_line_count} variants, "
                f"{len(headers.samples)} samples processed successfully."
            )

        if self.logger.isEnabledFor(logging.DEBUG):
            maf_values = list(snp_maf_cache.values())
            entropy_values = list(snp_entropy_cache.values())
            chrom_counts = {}
            for snp_data in snp_genotypes.values():
                chrom_counts[snp_data.chromosome] = (
                    chrom_counts.get(snp_data.chromosome, 0) + 1
                )

            self.logger.debug("VCF statistics:")
            self.logger.debug(
                f"  Chromosomes: {len(chrom_counts)} ({', '.join(sorted(chrom_counts.keys()))})"
            )
            self.logger.debug(
                f"  MAF range: {min(maf_values):.3f} - {max(maf_values):.3f}"
            )
            self.logger.debug(
                f"  Entropy range: {min(entropy_values):.3f} - {max(entropy_values):.3f}"
            )
            self.logger.debug(f"  Average MAF: {sum(maf_values) / len(maf_values):.3f}")
            self.logger.debug(
                f"  Average entropy: {sum(entropy_values) / len(entropy_values):.3f}"
            )

        return snp_genotypes, snp_maf_cache, snp_entropy_cache
