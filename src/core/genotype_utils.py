"""Genotype conversion and utility functions."""

from typing import List, Optional, Dict
import numpy as np
import math

__all__ = [
    "extract_gt",
    "gt_to_value",
    "calculate_maf",
    "is_het",
    "calculate_shannon_entropy",
    "calculate_snp_information_entropy",
    "calculate_snp_entropy_score",
]


def extract_gt(format_field: str, sample_field: str) -> str:
    """Extract GT subfield from sample field based on format specification."""
    if "GT" not in format_field:
        return "."
    gt_index = format_field.split(":").index("GT")
    parts = sample_field.split(":")
    return parts[gt_index] if gt_index < len(parts) else "."


def gt_to_value(gt: str, ploidy: Optional[int], convert_het: Optional[bool]) -> float:
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
    if ploidy is not None and count_alt == ploidy:
        return 1.0
    if ploidy is None:
        return np.nan
    return np.nan if convert_het else count_alt / float(ploidy)


def calculate_maf(geno: List[float], ploidy: Optional[int]) -> float:
    """Calculate Minor Allele Frequency (MAF) from genotype values."""
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
    """Return True if GT subfield represents a heterozygous call."""
    if not gt or gt == ".":
        return False
    alleles = [a for a in gt.replace("|", "/").split("/") if a != ""]
    return len(set(alleles)) > 1


def calculate_shannon_entropy(chromosome_counts: Dict[str, int]) -> float:
    """Calculate Shannon entropy for chromosome distribution."""
    if not chromosome_counts:
        return 0.0

    total_snps = sum(chromosome_counts.values())
    if total_snps == 0:
        return 0.0

    entropy = 0.0
    for count in chromosome_counts.values():
        if count > 0:
            probability = count / total_snps
            entropy -= probability * math.log2(probability)

    return entropy


def calculate_snp_information_entropy(genotypes: List[float]) -> float:
    """Calculate Shannon entropy for SNP information content based on genotype distribution."""
    if not genotypes:
        return 0.0

    # Count occurrences of each genotype category
    genotype_counts = {
        "ref": 0,  # 0.0 (homozygous reference)
        "het": 0,  # 0.5 (heterozygous)
        "alt": 0,  # 1.0 (homozygous alternate)
        "missing": 0,  # NaN (missing data)
    }

    total_valid = 0
    for gt in genotypes:
        if isinstance(gt, float) and np.isnan(gt):
            genotype_counts["missing"] += 1
        elif gt == 0.0:
            genotype_counts["ref"] += 1
        elif gt == 0.5:
            genotype_counts["het"] += 1
        elif gt == 1.0:
            genotype_counts["alt"] += 1
        else:
            # Handle any other unexpected values as missing
            genotype_counts["missing"] += 1

        total_valid += 1

    if total_valid == 0:
        return 0.0

    # Calculate entropy
    entropy = 0.0
    for count in genotype_counts.values():
        if count > 0:
            probability = count / total_valid
            entropy -= probability * math.log2(probability)

    return entropy


def calculate_snp_entropy_score(
    genotypes: List[float],
    maf: float,
    weight_entropy: float = 0.7,
    weight_maf: float = 0.3,
    entropy: Optional[float] = None,
) -> float:
    """Calculate combined SNP scoring using Shannon entropy and MAF."""
    if not genotypes:
        return 0.0

        # Calculate entropy if not provided
    if entropy is None:
        entropy = calculate_snp_information_entropy(genotypes)

    # Normalize entropy to [0, 1] range (max entropy is log2(4) = 2.0)
    max_entropy = 2.0
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0.0

    # Normalize MAF to [0, 1] range (max MAF is 0.5)
    max_maf = 0.5
    normalized_maf = maf / max_maf if max_maf > 0 else 0.0

    # Calculate combined score
    score = weight_entropy * normalized_entropy + weight_maf * normalized_maf

    return score
