"""SNP selection and optimization algorithms."""

import sys
import logging
from typing import List, Set, Optional, Union
from dataclasses import dataclass
from itertools import combinations
from math import floor

import numpy as np
from numpy.typing import NDArray

from .distance_calculator import DistanceCalculator
from .vcf_parser import VCFData, LazyVCFData
from ..utils.memory_monitor import MemoryMonitor

# Type alias for VCF data that supports both regular and lazy loading
VCFDataType = Union[VCFData, LazyVCFData]

__all__ = ["OverlapConstraints", "BuildError", "SNPSelector", "VCFDataType"]


@dataclass
class OverlapConstraints:
    """Overlap constraints for alternative set generation."""

    max_number: int = -1
    max_fraction: float = -1.0


class BuildError(Exception):
    """Exception raised when SNP set building fails."""

    pass


class SNPSelector:
    """SNP selection and optimization algorithms."""

    def __init__(
        self,
        distance_calculator: DistanceCalculator,
        memory_monitor: MemoryMonitor,
        logger: logging.Logger,
    ):
        self.distance_calc = distance_calculator
        self.memory_monitor = memory_monitor
        self.logger = logger

    def select_first_snp(
        self, vcf_data: VCFDataType, excluded: Optional[Set[str]] = None
    ) -> str:
        """Select SNP with highest MAF not in excluded; tie-break by SNP ID."""
        excluded = excluded or set()
        candidates = []

        for sid, snp_data in vcf_data.snp_genotypes.items():
            if sid in excluded:
                continue
            maf = vcf_data.snp_maf_cache.get(sid, 0.0)
            candidates.append((sid, maf))

        if not candidates:
            sys.exit("No SNPs available for selection after exclusions.")

        # Sort by MAF desc, SNP ID asc
        candidates.sort(key=lambda x: (-x[1], x[0]))
        return candidates[0][0]

    def build_primary_set(
        self, vcf_data: VCFDataType, min_distance: int, excluded: Optional[Set[str]] = None
    ) -> List[str]:
        """Greedily build initial SNP set maximizing MAF under exclusions."""
        excluded = excluded or set()
        current_snp_set: List[str] = []
        current_snps_geno: List[List[float]] = []

        # Choose first SNP
        current_snp = self.select_first_snp(vcf_data, excluded=excluded)
        current_snps_geno.append(vcf_data.snp_genotypes[current_snp].genotypes)
        current_snp_set.append(current_snp)

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Primary SNP selection...")

        # Identify primary set of SNPs deterministically with tie-breakers
        while self.distance_calc.calc_min_distance(current_snps_geno) < min_distance:
            if self.logger.isEnabledFor(logging.INFO):
                self.logger.info(
                    f"Current SNP set contains {len(current_snp_set)} SNPs..."
                )
            parent_nodes_info = []
            for sid, snp_data in vcf_data.snp_genotypes.items():
                if (sid in current_snp_set) or (sid in excluded):
                    continue
                maf = vcf_data.snp_maf_cache.get(sid, 0.0)
                parent_nodes_info.append((sid, maf))

            if not parent_nodes_info:
                raise BuildError("Not enough polymorphic SNP to discriminate samples.")

            # Sort by MAF desc, SNP ID asc
            parent_nodes_info.sort(key=lambda x: (-x[1], x[0]))
            current_snp = parent_nodes_info[0][0]
            current_snps_geno.append(vcf_data.snp_genotypes[current_snp].genotypes)
            current_snp_set.append(current_snp)

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"After 1st step {len(current_snp_set)} primary SNP selected"
            )
        return current_snp_set

    def deterministic_eliminate(
        self, snp_set: List[str], vcf_data: VCFDataType, min_distance: int
    ) -> List[str]:
        """Greedy backward elimination preserving minimal distance constraint.

        Optimized to avoid repeated distance recomputation by precomputing per-SNP
        contributions to pairwise distances across the upper triangle and updating
        totals incrementally when SNPs are removed.
        """
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Backward one-by-one elimination...")

        # Monitor memory before creating large matrices
        self.memory_monitor.check_memory_and_warn("elimination algorithm start")

        optimized_ids: List[str] = list(snp_set)
        # Build matrix (rows=SNPs, cols=samples) in current order
        snps_matrix = np.array(
            [vcf_data.snp_genotypes[sid].genotypes for sid in optimized_ids]
        )
        num_snps, num_samples = snps_matrix.shape

        # Precompute per-row contributions to upper-triangle pairwise distances
        pair_contribs: List[NDArray[np.int32]] = []
        for r in range(num_snps):
            row = snps_matrix[r, :]
            row_contrib_parts: List[NDArray[np.int32]] = []
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
                or candidate_pair_dists.min() >= min_distance
            ):
                # Accept removal: update state
                optimized_ids.pop(idx)
                current_pair_dists = candidate_pair_dists
                if contrib_matrix.size:
                    contrib_matrix = np.delete(contrib_matrix, idx, axis=0)

        # Monitor memory after elimination algorithm and force cleanup
        self.memory_monitor.check_memory_and_warn("elimination algorithm complete")
        self.memory_monitor.force_garbage_collection()

        return optimized_ids

    def search_optimal_sets(
        self,
        vcf_data: VCFDataType,
        min_distance: int,
        n_sets: int,
        max_snps: int,
        overlap_constraints: OverlapConstraints,
    ) -> List[List[str]]:
        """Find up to n_sets minimal discriminating SNP lists under constraints."""
        # Base minimal set
        try:
            base_primary = self.build_primary_set(
                vcf_data, min_distance, excluded=set()
            )
        except BuildError:
            sys.exit("Not enough polymorphic SNP to discriminate samples. Exit.")
        base_minimal = self.deterministic_eliminate(
            base_primary, vcf_data, min_distance
        )

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
            overlap_constraints.max_number
            if (
                overlap_constraints.max_number is not None
                and overlap_constraints.max_number >= 0
            )
            else len(base_minimal)
        )
        allowed_max_frac = (
            floor(overlap_constraints.max_fraction * len(base_minimal))
            if (
                overlap_constraints.max_fraction is not None
                and overlap_constraints.max_fraction >= 0
            )
            else len(base_minimal)
        )
        allowed_max = min(allowed_max_num, allowed_max_frac)

        # If we need to cap overlap, exclude combinations of base SNPs to achieve it
        exclude_size = max(1, len(base_minimal) - allowed_max)
        for excl in combinations(sorted(base_minimal), exclude_size):
            try:
                alt_primary = self.build_primary_set(
                    vcf_data, min_distance, excluded=set(excl)
                )
            except BuildError:
                continue
            alt_minimal = self.deterministic_eliminate(
                alt_primary, vcf_data, min_distance
            )
            overlap = len(set(alt_minimal).intersection(base_minimal))
            if len(alt_minimal) == len(base_minimal) and (overlap <= allowed_max):
                add_set(alt_minimal)
            if len(unique_sets) >= n_sets:
                break

        # Fallback: if no cap requested (allowed_max equals base size), use single exclusions
        if len(unique_sets) < n_sets and allowed_max >= len(base_minimal):
            for sid in sorted(base_minimal):
                try:
                    alt_primary = self.build_primary_set(
                        vcf_data, min_distance, excluded={sid}
                    )
                except BuildError:
                    continue
                alt_minimal = self.deterministic_eliminate(
                    alt_primary, vcf_data, min_distance
                )
                overlap = len(set(alt_minimal).intersection(base_minimal))
                if len(alt_minimal) == len(base_minimal):
                    add_set(alt_minimal)
                if len(unique_sets) >= n_sets:
                    break

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"{len(unique_sets)} discriminating SNP sets selected.")

        # Optionally expand by PIC to reach max_snps
        best_snp_sets_final = []
        for si, s in enumerate(unique_sets[:n_sets], start=1):
            orig_snp_number = len(s)
            if max_snps > len(s):
                snp_maf = {
                    sid: vcf_data.snp_maf_cache[sid]
                    for sid in vcf_data.snp_genotypes
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
                n_snps_to_add = max_snps - len(s)
                s = list(s) + [i[0] for i in snp_pic[:n_snps_to_add]]
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        f"{n_snps_to_add} SNPs added to set {si} "
                        f"(original set contains {orig_snp_number} SNPs, "
                        f"total number of SNPs: {len(s)})."
                    )
                best_snp_sets_final.append(s)
            else:
                best_snp_sets_final.append(list(s))
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        f"Addition of SNPs to discriminating set {si} "
                        f"is not required (total number of SNPs: {len(s)})."
                    )
        return best_snp_sets_final
