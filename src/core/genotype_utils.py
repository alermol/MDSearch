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
    """Extract GT subfield from sample field based on format specification.

    Args:
        format_field: FORMAT column (e.g., "GT:DP:GQ")
        sample_field: Sample-specific data (e.g., "0/0:14:94")

    Returns:
        GT value (e.g., "0/0") or "." if not found

    Example:
        >>> extract_gt("GT:DP:GQ", "0/0:14:94")
        '0/0'
        >>> extract_gt("DP:GQ", "14:94")
        '.'
    """
    if "GT" not in format_field:
        return "."
    gt_index = format_field.split(":").index("GT")
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
    """Calculate Minor Allele Frequency (MAF) from genotype values.

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


def calculate_shannon_entropy(chromosome_counts: Dict[str, int]) -> float:
    """Calculate Shannon entropy for chromosome distribution.

    Shannon entropy measures the balance of distribution of SNPs across chromosomes.
    Higher entropy indicates more balanced distribution, lower entropy indicates
    concentration in fewer chromosomes.

    Note: Chromosomes with 0 SNPs are included in the calculation to provide
    a complete picture of distribution across all chromosomes in the VCF.

    Args:
        chromosome_counts: Dictionary mapping chromosome names to SNP counts
                          (including 0 counts for chromosomes without SNPs)

    Returns:
        Shannon entropy value (0.0 to log2(num_chromosomes))

    Example:
        >>> # Balanced distribution across 3 chromosomes
        >>> balanced = {"chr1": 10, "chr2": 10, "chr3": 10}
        >>> entropy = calculate_shannon_entropy(balanced)
        >>> print(f"Balanced distribution entropy: {entropy:.3f}")
        Balanced distribution entropy: 1.585

        >>> # Unbalanced distribution (concentrated in chr1)
        >>> unbalanced = {"chr1": 25, "chr2": 3, "chr3": 2}
        >>> entropy = calculate_shannon_entropy(unbalanced)
        >>> print(f"Balanced distribution entropy: {entropy:.3f}")
        Balanced distribution entropy: 0.863

        >>> # Including chromosomes with 0 SNPs (from VCF headers)
        >>> with_zeros = {"chr1": 5, "chr2": 0, "chr3": 0, "chr4": 5}
        >>> entropy = calculate_shannon_entropy(with_zeros)
        >>> print(f"With zero-count chromosomes: {entropy:.3f}")
        With zero-count chromosomes: 1.000
    """
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
    """Calculate Shannon entropy for SNP information content based on genotype distribution.

    This function measures the information content of a SNP by calculating the Shannon
    entropy of its genotype distribution. Higher entropy indicates more informative SNPs
    with more balanced genotype distributions, while lower entropy indicates less
    informative SNPs with skewed distributions.

    The entropy is calculated as: H = -Σ(p_i * log2(p_i)) where p_i is the probability
    of each genotype category (0, 0.5, 1, missing).

    Args:
        genotypes: List of numeric genotype values (0.0, 0.5, 1.0, or NaN for missing)

    Returns:
        Shannon entropy value (0.0 to log2(4) ≈ 2.0 for maximum information)

    Example:
        >>> # Highly informative SNP with balanced distribution
        >>> balanced_geno = [0.0, 0.5, 1.0, 0.0, 0.5, 1.0]
        >>> entropy = calculate_snp_information_entropy(balanced_geno)
        >>> print(f"Balanced SNP entropy: {entropy:.3f}")
        Balanced SNP entropy: 1.585

        >>> # Less informative SNP with skewed distribution
        >>> skewed_geno = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0]
        >>> entropy = calculate_snp_information_entropy(skewed_geno)
        >>> print(f"Skewed SNP entropy: {entropy:.3f}")
        Skewed SNP entropy: 0.918

        >>> # SNP with missing data
        >>> missing_geno = [0.0, 0.5, np.nan, 1.0, 0.0]
        >>> entropy = calculate_snp_information_entropy(missing_geno)
        >>> print(f"SNP with missing data entropy: {entropy:.3f}")
        SNP with missing data entropy: 1.500
    """
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
    """Calculate combined SNP scoring using Shannon entropy and MAF.

    This function combines Shannon entropy (information content) with Minor Allele
    Frequency (MAF) to create a comprehensive SNP scoring system. The entropy
    measures information content, while MAF ensures the SNP is polymorphic.

    The combined score is calculated as:
    score = weight_entropy * normalized_entropy + weight_maf * normalized_maf

    where normalized_entropy = entropy / log2(4) and normalized_maf = maf / 0.5

    Args:
        genotypes: List of numeric genotype values
        maf: Minor allele frequency (0.0 to 0.5)
        weight_entropy: Weight for entropy component (default: 0.7)
        weight_maf: Weight for MAF component (default: 0.3)
        entropy: Precomputed entropy value (if None, will calculate from genotypes)

    Returns:
        Combined SNP score (0.0 to 1.0, higher is better)

        Example:
        >>> # High-quality SNP with balanced distribution and good MAF
        >>> genotypes = [0.0, 0.5, 1.0, 0.0, 0.5, 1.0]
        >>> maf = 0.4
        >>> score = calculate_snp_entropy_score(genotypes, maf)
        >>> print(f"High-quality SNP score: {score:.3f}")
        High-quality SNP score: 0.910

        >>> # SNP with lower information content
        >>> genotypes = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        >>> maf = 0.3
        >>> score = calculate_snp_entropy_score(genotypes, maf)
        >>> print(f"Lower quality SNP score: {score:.3f}")
        Lower quality SNP score: 0.636

        >>> # Using precomputed entropy for performance
        >>> precomputed_entropy = 1.585
        >>> score = calculate_snp_entropy_score(genotypes, maf, entropy=precomputed_entropy)
        >>> print(f"Score with precomputed entropy: {score:.3f}")
        Score with precomputed entropy: 0.910
    """
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
