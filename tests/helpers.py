import os
import shutil
from pathlib import Path
from typing import List, Dict
import numpy as np


def _artifacts_enabled() -> bool:
    value = os.getenv("MDSEARCH_SAVE_VCFS", "0")
    return value.lower() in ("1", "true", "yes", "on")


def _project_root() -> Path:
    # helpers.py resides in tests/, go one level up
    return Path(__file__).resolve().parents[1]


def _artifacts_dir(subdir: str = "") -> Path:
    d = _project_root() / "saved_vcfs"
    if subdir:
        d = d / subdir
    return d


def save_artifact(src: Path, subdir: str = ""):
    if not _artifacts_enabled():
        return
    src = Path(src)
    dst_dir = _artifacts_dir(subdir)
    dst_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst_dir / src.name)


def write_vcf(path: Path, samples: List[str], variants: List[Dict]):
    """
    Write a minimal VCF with provided variants.

    Each variant dict must contain keys:
      - chrom (str)
      - pos (int)
      - id (str)
      - ref (str)
      - alt (str)
      - genotypes (List[str]) aligned to samples order
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=mdsearch-tests\n")
        f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n"
        )
        for v in variants:
            line = [
                v.get("chrom", "1"),
                str(v.get("pos", 1)),
                v["id"],
                v.get("ref", "A"),
                v.get("alt", "T"),
                ".",
                "PASS",
                ".",
                "GT",
            ]
            line += v["genotypes"]
            f.write("\t".join(line) + "\n")
    # Optionally save generated VCF (original/expected) as artifact
    save_artifact(path)


def genotype_to_value(geno: str, ploidy: int, convert_het: bool):
    if "." in geno:
        return np.nan
    digits = [int(x) for x in geno.replace("|", "/").split("/")]
    # Handle haploid representation like "0" or "1"
    if len(digits) == 1 and ploidy == 1:
        return float(digits[0])
    count_alt = sum(1 for x in digits if x == 1)
    if count_alt == 0:
        return 0.0
    if count_alt == ploidy:
        return 1.0
    return np.nan if convert_het else count_alt / float(ploidy)


def parse_matrix_for_distance(
    vcf_path: Path, ploidy: int, convert_het: bool
) -> List[np.ndarray]:
    """Return list of arrays (variants x samples) with numeric genotypes like MDSearch uses.

    Supports both simple FORMAT with GT only and multi-field FORMAT like GT:DP:GQ by
    extracting only the GT subfield per sample.
    """
    matrices: List[np.ndarray] = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            format_field = parts[8] if len(parts) > 8 else "GT"
            keys = format_field.split(":") if format_field else []
            gt_index = keys.index("GT") if "GT" in keys else -1
            sample_fields = parts[9:]
            gt_values = []
            for sf in sample_fields:
                if gt_index == -1:
                    gt = sf
                else:
                    sub = sf.split(":")
                    gt = sub[gt_index] if gt_index < len(sub) else "."
                gt_values.append(genotype_to_value(gt, ploidy, convert_het))
            matrices.append(np.array(gt_values, dtype=float))
    return matrices


def calculate_min_hamming_distance(
    vcf_path: Path, ploidy: int, convert_het: bool
) -> int:
    matrices = parse_matrix_for_distance(vcf_path, ploidy, convert_het)
    if not matrices:
        return 0
    snps_array = np.vstack(matrices)  # shape: (num_snps, num_samples)
    n_samples = snps_array.shape[1]
    distances = np.zeros((n_samples, n_samples), dtype=float)
    for i in range(n_samples):
        for j in range(n_samples):
            col_i = snps_array[:, i]
            col_j = snps_array[:, j]
            valid_mask = ~np.isnan(col_i) & ~np.isnan(col_j)
            distance = np.nansum(col_i[valid_mask] != col_j[valid_mask])
            distances[i, j] = distance
    tri = distances[np.triu_indices(n_samples, k=1)]
    return int(np.min(tri)) if tri.size else 0


def get_snp_ids(vcf_path: Path) -> List[str]:
    ids: List[str] = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ids.append(line.rstrip().split("\t")[2])
    return ids


def assert_discriminative(
    vcf_path: Path, ploidy: int, min_dist: int, convert_het: bool
):
    observed = calculate_min_hamming_distance(vcf_path, ploidy, convert_het)
    assert (
        observed >= min_dist
    ), f"Min Hamming distance {observed} < required {min_dist}"


def save_out_prefix_vcfs(out_prefix: Path, subdir: str = ""):
    """Copy mdsearch outputs like <prefix>_1.vcf, <prefix>_2.vcf and summary.tsv into saved_vcfs/ if enabled."""
    if not _artifacts_enabled():
        return
    out_prefix = Path(out_prefix)
    out_dir = out_prefix.parent
    stem = out_prefix.name
    for f in out_dir.glob(f"{stem}_*.vcf"):
        save_artifact(f, subdir=subdir)
    for f in out_dir.glob(f"summary.tsv"):
        save_artifact(f, subdir=subdir)
