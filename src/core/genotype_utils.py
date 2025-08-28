"""Genotype conversion and utility functions."""

from typing import List, Optional
import numpy as np


def extract_gt(format_field: str, sample_field: str) -> str:
    """Extract GT subfield from sample data."""
    keys = format_field.split(":") if format_field else []
    if "GT" not in keys:
        return "."
    gt_index = keys.index("GT")
    parts = sample_field.split(":")
    return parts[gt_index] if gt_index < len(parts) else "."


def gt_to_value(
    gt: str, ploidy: Optional[int], convert_het: Optional[bool]
) -> float:
    """Convert GT string to numeric value based on ploidy and het conversion."""
    if "." in gt:
        return np.nan
    tokens = [t for t in gt.replace("|", "/").split("/") if t != ""]
    try:
        alleles = [int(x) for x in tokens]
    except ValueError:
        return np.nan
    if len(alleles) == 1 and ploidy == 1:
        return float(alleles[0])
    count_alt = sum(1 for x in alleles if x == 1)
    if count_alt == 0:
        return 0.0
    if count_alt == ploidy:
        return 1.0
    return np.nan if convert_het else count_alt / float(ploidy)


def calculate_maf(geno: List[float], ploidy: Optional[int]) -> float:
    """Compute minor allele frequency for a SNP given numeric genotypes."""
    allele0 = 0.0
    allele1 = 0.0
    valid_genotypes = 0
    for i in geno:
        if isinstance(i, float) and np.isnan(i):
            continue
        valid_genotypes += 1
        if i == 0:
            allele0 += ploidy
        elif i == 1:
            allele1 += ploidy
        else:
            allele1 += i * ploidy
            allele0 += ploidy - (i * ploidy)
    total_alleles = valid_genotypes * ploidy
    if total_alleles == 0:
        return 0.0
    return min(allele0 / total_alleles, allele1 / total_alleles)


def is_het(gt: str) -> bool:
    """Return True if GT subfield represents a heterozygous call."""
    if not gt or gt == ".":
        return False
    alleles = [a for a in gt.replace("|", "/").split("/") if a != ""]
    return len(set(alleles)) > 1
