"""Genotype conversion and utility functions."""

from typing import List, Optional
import numpy as np

__all__ = ["extract_gt", "gt_to_value", "calculate_maf", "is_het"]


def extract_gt(format_field: str, sample_field: str) -> str:
    """Extract GT subfield from sample data.

    Args:
        format_field: FORMAT field from VCF (e.g., "GT:DP:GQ")
        sample_field: Sample data field (e.g., "0/1:10:99")

    Returns:
        GT value as string, or "." if GT not found

    Example:
        >>> extract_gt("GT:DP:GQ", "0/1:10:99")
        '0/1'
        >>> extract_gt("DP:GQ", "10:99")  # No GT field
        '.'
        >>> extract_gt("GT", "0/0")
        '0/0'
    """
    keys = format_field.split(":") if format_field else []
    if "GT" not in keys:
        return "."
    gt_index = keys.index("GT")
    parts = sample_field.split(":")
    return parts[gt_index] if gt_index < len(parts) else "."


def gt_to_value(gt: str, ploidy: Optional[int], convert_het: Optional[bool]) -> float:
    """Convert GT string to numeric value based on ploidy and het conversion.

    Args:
        gt: Genotype string (e.g., "0/1", "1/1", "0/0")
        ploidy: Ploidy level (e.g., 2 for diploid)
        convert_het: Whether to convert heterozygous calls to missing

    Returns:
        Numeric genotype value (0.0, 1.0, or NaN for missing/het)

    Example:
        >>> gt_to_value("0/0", 2, False)  # Homozygous reference
        0.0
        >>> gt_to_value("1/1", 2, False)  # Homozygous alternate
        1.0
        >>> gt_to_value("0/1", 2, False)  # Heterozygous
        0.5
        >>> gt_to_value("0/1", 2, True)   # Het converted to missing
        nan
        >>> gt_to_value(".", 2, False)    # Missing
        nan
        >>> gt_to_value("0", 1, False)    # Haploid
        0.0
    """
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
    if ploidy is not None and count_alt == ploidy:
        return 1.0
    if ploidy is None:
        return np.nan
    return np.nan if convert_het else count_alt / float(ploidy)


def calculate_maf(geno: List[float], ploidy: Optional[int]) -> float:
    """Compute minor allele frequency for a SNP given numeric genotypes.

    Args:
        geno: List of numeric genotype values
        ploidy: Ploidy level (e.g., 2 for diploid)

    Returns:
        Minor allele frequency (0.0 to 0.5)

    Example:
        >>> # Diploid genotypes: 0=ref/ref, 1=alt/alt, 0.5=ref/alt
        >>> genotypes = [0.0, 0.5, 1.0, 0.0, 0.5]
        >>> maf = calculate_maf(genotypes, ploidy=2)
        >>> print(f"Minor allele frequency: {maf:.2f}")
        Minor allele frequency: 0.30

        >>> # Haploid genotypes: 0=ref, 1=alt
        >>> hap_genotypes = [0, 1, 0, 1, 1]
        >>> maf = calculate_maf(hap_genotypes, ploidy=1)
        >>> print(f"Minor allele frequency: {maf:.2f}")
        Minor allele frequency: 0.40
    """
    if ploidy is None:
        return 0.0

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
    """Return True if GT subfield represents a heterozygous call.

    Args:
        gt: Genotype string to check

    Returns:
        True if heterozygous, False otherwise

    Example:
        >>> is_het("0/1")    # Heterozygous
        True
        >>> is_het("0|1")    # Phased heterozygous
        True
        >>> is_het("0/0")    # Homozygous reference
        False
        >>> is_het("1/1")    # Homozygous alternate
        False
        >>> is_het(".")      # Missing
        False
        >>> is_het("")       # Empty
        False
    """
    if not gt or gt == ".":
        return False
    alleles = [a for a in gt.replace("|", "/").split("/") if a != ""]
    return len(set(alleles)) > 1
