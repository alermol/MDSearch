from pathlib import Path
import sys
import logging
import time
import json
from typing import List, Dict, Set, Optional, Union, Any

import numpy as np
from itertools import combinations
from math import floor


def parser_resolve_path(path: str) -> Path:
    """Resolve CLI-provided path string to an absolute Path."""
    return Path(path).resolve()


class MDSearch:
    """Search for minimal SNP sets that discriminate samples.

    The algorithm builds a primary set by greedily maximizing MAF, then
    prunes it by deterministic elimination while maintaining a minimal
    Hamming distance between all sample pairs. Multiple alternative sets
    can be generated with optional overlap constraints.
    """

    def __init__(
        self,
        in_vcf: Path,
        out_vcf: Path,
        ploidy: Optional[int] = None,
        max_snps: Optional[int] = None,
        min_dist: Optional[int] = None,
        convert_het: Optional[bool] = None,
        n_sets: Optional[int] = None,
        overlap_max_number: int = -1,
        overlap_max_fraction: float = -1.0,
        verbose: bool = True,
        log_level: Optional[str] = None,
        log_format: str = "text",
        summary_tsv: Optional[Path] = None,
    ) -> None:
        self.in_vcf = in_vcf
        self.out_vcf = out_vcf
        self.ploidy = ploidy
        self.max_snps = max_snps
        self.min_dist = min_dist
        self.convert_het = convert_het
        self.n_sets = n_sets
        self.overlap_max_number = overlap_max_number
        self.overlap_max_fraction = overlap_max_fraction
        self.verbose = verbose
        self.summary_tsv = Path(summary_tsv) if summary_tsv else None
        # Configure logger
        self.logger = logging.getLogger("mdsearch")
        # Avoid duplicate handlers if multiple instances
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            if log_format == "json":

                class JsonFormatter(logging.Formatter):
                    def format(self, record: logging.LogRecord) -> str:  # type: ignore[override]
                        payload: Dict[str, Any] = {
                            "level": record.levelname,
                            "message": record.getMessage(),
                            "time": time.strftime(
                                "%Y-%m-%dT%H:%M:%SZ", time.gmtime(record.created)
                            ),
                        }
                        # If message is JSON-like string, keep it as message
                        return json.dumps(payload, ensure_ascii=False)

                handler.setFormatter(JsonFormatter())
            else:
                handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
            self.logger.addHandler(handler)
            self.logger.propagate = False
        # Determine level
        level_name = (log_level or ("INFO" if verbose else "WARNING")).upper()
        level = getattr(logging, level_name, logging.INFO)
        self.logger.setLevel(level)

        # Validate mutually exclusive overlap constraints
        if (self.overlap_max_number is not None and self.overlap_max_number >= 0) and (
            self.overlap_max_fraction is not None and self.overlap_max_fraction >= 0
        ):
            sys.exit(
                "Specify only one of -oMx (max overlap number) or -oMf "
                "(max overlap fraction), not both."
            )

        # Validate CLI numeric constraints early
        if self.max_snps is not None and self.max_snps < 0:
            sys.exit("-ts (total SNPs) must be >= 0")
        if self.n_sets is None or self.n_sets < 1:
            sys.exit("-ns (number of sets) must be >= 1")

        # calculate target number of genotypes and create list containing genotype for each SNP
        self.snp_genotypes: Dict[str, List[Union[List[float], List[str]]]] = {}
        self.snp_maf_cache: Dict[str, float] = {}
        self._parse_and_validate_vcf()
        self.main()

    def _parse_and_validate_vcf(self) -> None:
        """Parse and validate VCF format with comprehensive error checking."""
        # First pass: validate headers and count samples
        header_line = None
        samples = []
        fileformat_found = False
        line_number = 0

        with open(self.in_vcf) as vcf:
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

        # Second pass: parse and validate data lines
        seen_ids: set[str] = set()
        data_line_count = 0

        with open(self.in_vcf) as vcf:
            line_number = 0
            for line in vcf:
                line_number += 1
                line = line.strip()

                if line.startswith("#") or not line:
                    continue

                data_line_count += 1
                parts = line.split("\t")

                # Validate line format
                expected_columns = 9 + len(samples)
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

                geno: list[float] = []
                for sample_idx, sample_field in enumerate(sample_fields):
                    gt = self._extract_gt(format_field, sample_field)

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

                    geno.append(self._gt_to_value(gt, self.ploidy, self.convert_het))

                # Store validated data
                self.snp_genotypes[snp_id] = [geno, sample_fields]
                # Cache MAF for this SNP
                self.snp_maf_cache[snp_id] = MDSearch._calculate_maf(
                    self.snp_genotypes[snp_id][0], self.ploidy
                )

        # Final validation
        if data_line_count == 0:
            sys.exit(
                "ERROR: No data lines found in VCF. "
                "Ensure VCF contains variant records."
            )

        if self.verbose:
            self.logger.info(
                f"VCF validation complete: {data_line_count} variants, "
                f"{len(samples)} samples processed successfully."
            )

    def _extract_gt(self, format_field: str, sample_field: str) -> str:
        """Extract GT subfield from sample data."""
        keys = format_field.split(":") if format_field else []
        if "GT" not in keys:
            return "."
        gt_index = keys.index("GT")
        parts = sample_field.split(":")
        return parts[gt_index] if gt_index < len(parts) else "."

    def _gt_to_value(
        self, gt: str, ploidy: Optional[int], convert_het: Optional[bool]
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

    @staticmethod
    def _calculate_maf(geno: List[float], ploidy: Optional[int]) -> float:
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

    def _calc_min_dist(self, snps: List[List[float]]) -> float:
        """Return minimal pairwise Hamming distance across samples for given SNPs.

        Optimized: compute only upper-triangle pairwise distances (j > i) and avoid
        self-pairs. Uses vectorized comparisons per anchor column.
        """
        if self.verbose:
            self.logger.info(
                f"Calculate pairwise distance based on {len(snps)} SNPs..."
            )
        snps_array = np.array([i for i in snps])  # shape: (num_snps, num_samples)
        num_samples = snps_array.shape[1]
        pairwise_distances: list[float] = []
        if num_samples <= 1:
            if self.verbose:
                self.logger.info("Distance between samples (min/med/avg/max): 0/0/0/0")
            return 0.0
        for i in range(num_samples - 1):
            col_i = snps_array[:, i][:, None]  # (num_snps, 1)
            rest = snps_array[:, i + 1 :]  # (num_snps, num_samples - i - 1)
            if rest.shape[1] == 0:
                continue
            valid_mask = (~np.isnan(col_i)) & (~np.isnan(rest))
            diffs = (col_i != rest) & valid_mask
            dists = diffs.sum(axis=0).astype(float)  # per pair distances vs column i
            pairwise_distances.extend(dists.tolist())
        res = np.array(pairwise_distances, dtype=float)
        if self.verbose:
            self.logger.info(
                "Distance between samples (min/med/avg/max): "
                f"{np.min(res)}/{np.median(res)}/{round(float(np.mean(res)), 1)}/{np.max(res)}"
            )
        return float(np.min(res))

    def select_first_snp(self, excluded: Optional[Set[str]] = None) -> str:
        """Select SNP with highest MAF not in excluded; tie-break by SNP ID."""
        # select SNP with max MAF; deterministic tie-breaker by SNP ID asc
        excluded = excluded or set()
        candidates = []
        for sid, sg in self.snp_genotypes.items():
            if sid in excluded:
                continue
            maf = self.snp_maf_cache.get(sid, 0.0)
            candidates.append((sid, maf))
        if not candidates:
            sys.exit("No SNPs available for selection after exclusions.")
        candidates.sort(key=lambda x: (-x[1], x[0]))
        return candidates[0][0]

    def is_het(self, gt: str) -> bool:
        """Return True if GT subfield represents a heterozygous call."""
        if not gt or gt == ".":
            return False
        alleles = [a for a in gt.replace("|", "/").split("/") if a != ""]
        return len(set(alleles)) > 1

    def _calc_min_dist_for_set_ids(self, snp_ids: List[str]) -> float:
        """Helper computing minimal distance for a list of SNP IDs."""
        if not snp_ids:
            return 0.0
        genos = [self.snp_genotypes[sid][0] for sid in snp_ids]
        return self._calc_min_dist(genos)

    class BuildError(Exception):
        pass

    def _build_primary_set(self, excluded: Optional[Set[str]] = None) -> List[str]:
        """Greedily build initial SNP set maximizing MAF under exclusions."""
        excluded = excluded or set()
        current_snp_set: List[str] = []
        current_snps_geno: List[List[float]] = []

        # choose first snp
        current_snp = self.select_first_snp(excluded=excluded)
        current_snps_geno.append(self.snp_genotypes[current_snp][0])
        current_snp_set.append(current_snp)

        if self.verbose:
            self.logger.info("Primary SNP selection...")

        # identify primary set of SNPs deterministically with tie-breakers
        while self._calc_min_dist(current_snps_geno) < self.min_dist:
            if self.verbose:
                self.logger.info(
                    f"Current SNP set contains {len(current_snp_set)} SNPs..."
                )
            parent_nodes_info = []
            for s, g in self.snp_genotypes.items():
                if (s in current_snp_set) or (s in excluded):
                    continue
                maf = self.snp_maf_cache.get(s, 0.0)
                parent_nodes_info.append((s, maf))
            if not parent_nodes_info:
                raise MDSearch.BuildError(
                    "Not enough polymorphic SNP to discriminate samples."
                )
            # Sort by MAF desc, SNP ID asc
            parent_nodes_info.sort(key=lambda x: (-x[1], x[0]))
            current_snp = parent_nodes_info[0][0]
            current_snps_geno.append(self.snp_genotypes[current_snp][0])
            current_snp_set.append(current_snp)

        if self.verbose:
            self.logger.info(
                f"After 1st step {len(current_snp_set)} primary SNP selected"
            )
        return current_snp_set

    def _deterministic_eliminate(self, snp_set: List[str]) -> List[str]:
        """Greedy backward elimination preserving minimal distance constraint.

        Optimized to avoid repeated distance recomputation by precomputing per-SNP
        contributions to pairwise distances across the upper triangle and updating
        totals incrementally when SNPs are removed.
        """
        if self.verbose:
            self.logger.info("Backward one-by-one elimination (deterministic)...")

        optimized_ids: List[str] = list(snp_set)
        # Build matrix (rows=SNPs, cols=samples) in current order
        snps_matrix = np.array([self.snp_genotypes[sid][0] for sid in optimized_ids])
        num_snps, num_samples = snps_matrix.shape

        # Precompute per-row contributions to upper-triangle pairwise distances
        pair_contribs: list[np.ndarray] = []
        for r in range(num_snps):
            row = snps_matrix[r, :]
            row_contrib_parts: list[np.ndarray] = []
            for i in range(num_samples - 1):
                vi = row[i]
                rest = row[i + 1 :]
                vi_col = np.full(rest.shape, vi)
                valid = (~np.isnan(vi_col)) & (~np.isnan(rest))
                diffs = (vi_col != rest) & valid
                row_contrib_parts.append(diffs.astype(np.int32))
            if row_contrib_parts:
                row_contrib = np.concatenate(row_contrib_parts)
            else:
                row_contrib = np.zeros(0, dtype=np.int32)
            pair_contribs.append(row_contrib)
        if pair_contribs:
            contrib_matrix = np.vstack(pair_contribs)  # shape: (num_snps, num_pairs)
            current_pair_dists = contrib_matrix.sum(axis=0).astype(np.int32)
        else:
            contrib_matrix = np.zeros((0, 0), dtype=np.int32)
            current_pair_dists = np.zeros(0, dtype=np.int32)

        # Greedy removal in stable order (SNP ID asc)
        for sid in sorted(list(optimized_ids)):
            # Skip if sid already removed
            if sid not in optimized_ids:
                continue
            idx = optimized_ids.index(sid)
            # Compute distances if this SNP is removed
            candidate_pair_dists = current_pair_dists - contrib_matrix[idx]
            if (
                candidate_pair_dists.size == 0
                or candidate_pair_dists.min() >= self.min_dist
            ):
                # Accept removal: update state
                optimized_ids.pop(idx)
                current_pair_dists = candidate_pair_dists
                if contrib_matrix.size:
                    contrib_matrix = np.delete(contrib_matrix, idx, axis=0)
        return optimized_ids

    def optimal_snp_set_search(self) -> List[List[str]]:
        """Find up to n_sets minimal discriminating SNP lists under constraints."""
        # Base minimal set
        try:
            base_primary = self._build_primary_set(excluded=set())
        except MDSearch.BuildError:
            sys.exit("Not enough polymorphic SNP to discriminate samples. Exit.")
        base_minimal = self._deterministic_eliminate(base_primary)

        # Enumerate alternatives by excluding each SNP from the base minimal set
        unique_sets = []
        seen = set()

        def add_set(s: List[str]) -> None:
            key = tuple(sorted(s))
            if key not in seen:
                seen.add(key)
                unique_sets.append(list(key))

        add_set(base_minimal)

        # Determine allowed maximum overlap (default = no cap)
        allowed_max_num = (
            self.overlap_max_number
            if (self.overlap_max_number is not None and self.overlap_max_number >= 0)
            else len(base_minimal)
        )
        allowed_max_frac = (
            floor(self.overlap_max_fraction * len(base_minimal))
            if (
                self.overlap_max_fraction is not None and self.overlap_max_fraction >= 0
            )
            else len(base_minimal)
        )
        allowed_max = min(allowed_max_num, allowed_max_frac)

        # If we need to cap overlap, exclude combinations of base SNPs to achieve it
        exclude_size = max(1, len(base_minimal) - allowed_max)
        for excl in combinations(sorted(base_minimal), exclude_size):
            try:
                alt_primary = self._build_primary_set(excluded=set(excl))
            except MDSearch.BuildError:
                continue
            alt_minimal = self._deterministic_eliminate(alt_primary)
            overlap = len(set(alt_minimal).intersection(base_minimal))
            if len(alt_minimal) == len(base_minimal) and (overlap <= allowed_max):
                add_set(alt_minimal)
            if len(unique_sets) >= self.n_sets:
                break
        # Fallback: if no cap requested (allowed_max equals base size), use single exclusions
        if len(unique_sets) < self.n_sets and allowed_max >= len(base_minimal):
            for sid in sorted(base_minimal):
                try:
                    alt_primary = self._build_primary_set(excluded={sid})
                except MDSearch.BuildError:
                    continue
                alt_minimal = self._deterministic_eliminate(alt_primary)
                overlap = len(set(alt_minimal).intersection(base_minimal))
                if len(alt_minimal) == len(base_minimal):
                    add_set(alt_minimal)
                if len(unique_sets) >= self.n_sets:
                    break

        if self.verbose:
            self.logger.info(f"{len(unique_sets)} discriminating SNP sets selected.")

        # Optionally expand by PIC to reach max_snps
        best_snp_sets_final = []
        for si, s in enumerate(unique_sets[: self.n_sets], start=1):
            orig_snp_number = len(s)
            if self.max_snps > len(s):
                snp_maf = {
                    sid: self.snp_maf_cache[sid]
                    for sid in self.snp_genotypes
                    if sid not in s
                }
                snp_pic = sorted(
                    [
                        (
                            sid,
                            1 - ((maf**2) + ((1 - maf) ** 2)),
                        )
                        for sid, maf in snp_maf.items()
                    ],
                    key=lambda x: x[1],
                )[::-1]
                n_snps_to_add = self.max_snps - len(s)
                s = list(s) + [i[0] for i in snp_pic[:n_snps_to_add]]
                if self.verbose:
                    self.logger.info(
                        f"{n_snps_to_add} SNPs added to set {si} "
                        f"(original set contains {orig_snp_number} SNPs, "
                        f"total number of SNPs: {len(s)})."
                    )
                best_snp_sets_final.append(s)
            else:
                best_snp_sets_final.append(list(s))
                if self.verbose:
                    self.logger.info(
                        f"Addition of SNPs to discriminating set {si} "
                        f"is not required (total number of SNPs: {len(s)})."
                    )
        return best_snp_sets_final

    def write_vcf(self, snp_sets_list: List[List[str]]) -> None:
        """Write each SNP set to a separate VCF suffixed by index starting at 1."""
        for si, s in enumerate(snp_sets_list, start=1):
            with open(self.in_vcf) as invcf, open(
                f"{self.out_vcf}_{si}.vcf", "w"
            ) as outvcf:
                for vcf_line in invcf:
                    if vcf_line.startswith("#"):
                        outvcf.write(vcf_line)
                    else:
                        line = vcf_line.strip().split("\t")
                        if (line[2] in s) and self.convert_het:
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
                                sep.join(["."] * self.ploidy)
                                if (self.ploidy and self.ploidy > 1)
                                else "."
                            )
                            new_samples = []
                            for sample_field in line[9:]:
                                if gt_index == -1:
                                    new_samples.append(sample_field)
                                    continue
                                parts = sample_field.split(":")
                                gt_val = (
                                    parts[gt_index] if gt_index < len(parts) else "."
                                )
                                if self.is_het(gt_val):
                                    if gt_index < len(parts):
                                        parts[gt_index] = missing_gt
                                        new_samples.append(":".join(parts))
                                    else:
                                        new_samples.append(sample_field)
                                else:
                                    new_samples.append(sample_field)
                            line = line[:9] + new_samples
                            outvcf.write("\t".join(line) + "\n")
                        elif (line[2] in s) and (not self.convert_het):
                            outvcf.write(vcf_line)
                        else:
                            continue

    def write_summary_tsv(self, snp_sets_list: List[List[str]]) -> None:
        """Emit a TSV summary with per-set stats if a summary path is provided."""
        if not self.summary_tsv:
            return
        header = [
            "set_index",
            "output_vcf",
            "num_snps",
            "min_distance",
            "snp_ids",
        ]
        lines = ["\t".join(header)]
        for si, s in enumerate(snp_sets_list, start=1):
            out_vcf = f"{self.out_vcf}_{si}.vcf"
            min_d = self._calc_min_dist_for_set_ids(s)
            snp_ids_str = ",".join(sorted(s))
            row = [str(si), out_vcf, str(len(s)), str(int(min_d)), snp_ids_str]
            lines.append("\t".join(row))
        self.summary_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(self.summary_tsv, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        if self.verbose:
            self.logger.info(f"Summary TSV written: {self.summary_tsv}")

    def main(self) -> None:
        """Entrypoint: search optimal SNP sets and write them to VCF files."""
        selected_snps = self.optimal_snp_set_search()

        # Log selected SNP details at INFO level when not using TSV summary
        if not self.summary_tsv:
            for si, s in enumerate(selected_snps, start=1):
                min_d = self._calc_min_dist_for_set_ids(s)
                snp_ids_str = ",".join(sorted(s))
                self.logger.info(
                    f"Set {si}: {len(s)} SNPs, min_distance={int(min_d)}, "
                    f"SNP_IDs=[{snp_ids_str}]"
                )

        if self.verbose:
            self.logger.info("Writing selected SNPs in VCF...")
        self.write_vcf(selected_snps)
        if self.verbose:
            self.logger.info("Done")
        # Optional summary TSV
        self.write_summary_tsv(selected_snps)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Select minimal discriminatory SNP sets from a VCF given a minimal "
            "pairwise Hamming distance; supports multiple alternative sets and "
            "optional constraints."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Notes: Input VCF must be bi-allelic with SNP IDs present. "
            "Multiallelics are rejected via ALT or GT."
        ),
    )

    parser.add_argument(
        "ivcf",
        help="Input VCF file (bi-allelic, with SNP IDs)",
        type=parser_resolve_path,
        metavar="IVCF",
    )
    parser.add_argument(
        "ovcf_prefix",
        help="Prefix for output VCF(s)",
        type=parser_resolve_path,
        metavar="OVCF_PREFIX",
    )
    # Legacy flags -s/-t/-c were removed; no-ops kept out of argparse to avoid confusion
    parser.add_argument(
        "-pl", help="VCF ploidy (default: 2)", default=2, type=int, metavar="PLOIDY"
    )
    parser.add_argument(
        "-ts",
        help="Total SNPs in output set (0 = keep minimal; default: 0)",
        default=0,
        type=int,
        metavar="TOTAL_SNP",
    )
    parser.add_argument(
        "-md",
        help="Minimal Hamming distance between samples (default: 1)",
        default=1,
        type=int,
        metavar="MIN_DIST",
    )
    parser.add_argument(
        "-ch",
        help="Convert heterozygous calls into NA (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-ns",
        help="Number of distinct SNP sets in output (default: 1)",
        default=1,
        type=int,
        metavar="N_SETS",
    )
    parser.add_argument(
        "-oMx",
        help=(
            "Maximum overlap count allowed with the base minimal set for "
            "alternative sets (-1 = unlimited; default: -1)"
        ),
        default=-1,
        type=int,
        metavar="OVERLAP_MAX_N",
    )
    parser.add_argument(
        "-oMf",
        help=(
            "Maximum overlap fraction allowed with the base minimal set for "
            "alternative sets (-1 = unlimited; default: -1.0)"
        ),
        default=-1.0,
        type=float,
        metavar="OVERLAP_MAX_FRAC",
    )
    parser.add_argument(
        "--quiet",
        help="Suppress progress output (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--log-level",
        help=(
            "Logging level (DEBUG, INFO, WARNING, ERROR); default depends on --quiet"
        ),
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default=None,
    )
    parser.add_argument(
        "--log-format",
        help="Logging format: text or json (default: text)",
        choices=["text", "json"],
        default="text",
    )
    parser.add_argument(
        "--summary-tsv",
        help=(
            "Write per-set summary TSV to this path (columns: set_index, output_vcf, "
            "num_snps, min_distance, snp_ids)."
        ),
        type=parser_resolve_path,
        default=None,
        metavar="SUMMARY_TSV",
    )

    args = parser.parse_args()

    MDSearch(
        in_vcf=args.ivcf,
        out_vcf=args.ovcf_prefix,
        ploidy=args.pl,
        max_snps=args.ts,
        min_dist=args.md,
        convert_het=args.ch,
        n_sets=args.ns,
        overlap_max_number=args.oMx,
        overlap_max_fraction=args.oMf,
        verbose=not args.quiet,
        log_level=args.log_level,
        log_format=args.log_format,
        summary_tsv=args.summary_tsv,
    )
