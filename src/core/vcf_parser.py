"""VCF file parsing and validation."""

import sys
import logging
from pathlib import Path
from typing import List, Dict, Set, Union, Any
from dataclasses import dataclass

from ..utils.memory_monitor import MemoryMonitor
from .genotype_utils import extract_gt, gt_to_value, calculate_maf


@dataclass
class VCFHeaders:
    """VCF header information."""
    samples: List[str]
    fileformat_found: bool


@dataclass 
class SNPData:
    """Data for a single SNP."""
    genotypes: List[float]
    sample_fields: List[str]
    maf: float


@dataclass
class VCFData:
    """Complete VCF data structure."""
    headers: VCFHeaders
    snp_genotypes: Dict[str, SNPData]
    snp_maf_cache: Dict[str, float]
    

class VCFParser:
    """Handles VCF file parsing and validation."""
    
    def __init__(self, memory_monitor: MemoryMonitor, logger: logging.Logger):
        self.memory_monitor = memory_monitor
        self.logger = logger
    
    def parse_and_validate(self, vcf_path: Path, ploidy: int, convert_het: bool) -> VCFData:
        """Parse VCF file and return structured data."""
        # First pass: validate headers and count samples
        headers = self._validate_headers(vcf_path)
        
        # Second pass: parse and validate data lines
        snp_genotypes, snp_maf_cache = self._validate_variants(
            vcf_path, headers, ploidy, convert_het
        )
        
        # Check memory usage after VCF parsing and warn for large datasets
        self.memory_monitor.check_memory_and_warn("VCF parsing")
        self.memory_monitor.warn_for_large_dataset(
            len(snp_genotypes), len(headers.samples)
        )
        
        return VCFData(
            headers=headers,
            snp_genotypes=snp_genotypes,
            snp_maf_cache=snp_maf_cache
        )
    
    def _validate_headers(self, vcf_path: Path) -> VCFHeaders:
        """Validate VCF headers and extract sample information."""
        header_line = None
        samples = []
        fileformat_found = False
        line_number = 0

        with open(vcf_path) as vcf:
            for line in vcf:
                line_number += 1
                line = line.strip()

                if not line:
                    continue

                if line.startswith("##fileformat=VCF"):
                    fileformat_found = True
                elif line.startswith("#CHROM"):
                    header_line = line
                    parts = line.split("\t")
                    if len(parts) < 9:
                        sys.exit(
                            f"ERROR: Line {line_number}: Invalid header format. "
                            f"Expected at least 9 columns, found {len(parts)}."
                        )
                    samples = parts[9:]
                    break
                elif line.startswith("#"):
                    continue  # Skip other header lines
                else:
                    sys.exit(
                        f"ERROR: Line {line_number}: Found data line before required #CHROM header. "
                        "VCF must have proper header structure."
                    )

        # Validate headers
        if not fileformat_found:
            sys.exit(
                "ERROR: Missing required VCF header '##fileformat=VCFv4.x'. "
                "Ensure input is a valid VCF file."
            )

        if header_line is None:
            sys.exit(
                "ERROR: Missing required #CHROM header line. "
                "VCF must have column headers."
            )

        # Validate minimum sample count
        if len(samples) < 2:
            sys.exit(
                f"ERROR: Insufficient samples for distance calculation. "
                f"Found {len(samples)} samples, need at least 2."
            )
            
        return VCFHeaders(samples=samples, fileformat_found=fileformat_found)
    
    def _validate_variants(
        self, 
        vcf_path: Path, 
        headers: VCFHeaders, 
        ploidy: int, 
        convert_het: bool
    ) -> tuple[Dict[str, SNPData], Dict[str, float]]:
        """Parse and validate variant lines."""
        seen_ids: Set[str] = set()
        data_line_count = 0
        snp_genotypes: Dict[str, SNPData] = {}
        snp_maf_cache: Dict[str, float] = {}

        with open(vcf_path) as vcf:
            line_number = 0
            for line in vcf:
                line_number += 1
                line = line.strip()

                if line.startswith("#") or not line:
                    continue

                data_line_count += 1
                parts = line.split("\t")

                # Validate line format
                expected_columns = 9 + len(headers.samples)
                if len(parts) != expected_columns:
                    sys.exit(
                        f"ERROR: Line {line_number}: Invalid number of columns. "
                        f"Expected {expected_columns}, found {len(parts)}."
                    )

                # Validate required columns exist
                if len(parts) < 5:
                    sys.exit(
                        f"ERROR: Line {line_number}: Missing required VCF columns "
                        "(CHROM, POS, ID, REF, ALT)."
                    )

                # Validate SNP ID
                snp_id = parts[2]
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

                # ALT-based multiallelic detection
                alt_field = parts[4]
                if "," in alt_field:
                    sys.exit(
                        f"ERROR: Line {line_number}: Multiallelic site detected (ALT='{alt_field}'). "
                        "Filter or split multiallelic sites before processing."
                    )

                # Validate FORMAT and sample fields
                format_field = parts[8] if len(parts) > 8 else "GT"
                sample_fields = parts[9:]

                geno: List[float] = []
                for sample_idx, sample_field in enumerate(sample_fields):
                    gt = extract_gt(format_field, sample_field)

                    # Detect multiallelic allele indices in GT
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

                    geno.append(gt_to_value(gt, ploidy, convert_het))

                # Store validated data
                maf = calculate_maf(geno, ploidy)
                snp_data = SNPData(
                    genotypes=geno,
                    sample_fields=sample_fields,
                    maf=maf
                )
                snp_genotypes[snp_id] = snp_data
                snp_maf_cache[snp_id] = maf

        # Final validation
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
            
        return snp_genotypes, snp_maf_cache
