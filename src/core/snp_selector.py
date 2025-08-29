"""SNP selection and optimization algorithms."""

import sys
import logging
from typing import List, Set, Optional

import numpy as np
from numpy.typing import NDArray

from .distance_calculator import DistanceCalculator
from .vcf_parser import VCFData
from ..utils.memory_monitor import MemoryMonitor

# VCFDataType alias removed - no longer needed without lazy loading

__all__ = ["BuildError", "SNPSelector"]


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
        """Initialize SNP selector with required components.

        Args:
            distance_calculator: DistanceCalculator instance for computing distances
            memory_monitor: MemoryMonitor instance for tracking memory usage
            logger: Logger instance for output

        Example:
            >>> from src.core.distance_calculator import DistanceCalculator
            >>> from src.utils.memory_monitor import MemoryMonitor
            >>> from src.utils.logging_setup import setup_logger
            >>> logger = setup_logger("snp_selector")
            >>> memory_monitor = MemoryMonitor(logger)
            >>> distance_calc = DistanceCalculator(memory_monitor, logger)
            >>> selector = SNPSelector(distance_calc, memory_monitor, logger)
        """
        self.distance_calc = distance_calculator
        self.memory_monitor = memory_monitor
        self.logger = logger

    def select_first_discriminatory_snp(
        self, vcf_data: VCFData, excluded: Optional[Set[str]] = None
    ) -> str:
        """Select SNP with highest MAF not in excluded; tie-break by SNP ID.

        Args:
            vcf_data: VCF data containing SNP information
            excluded: Set of SNP IDs to exclude from selection

        Returns:
            SNP ID with highest MAF (or lexicographically first if tied)

        Raises:
            BuildError: If no SNPs available after exclusions

        Example:
            >>> selector = SNPSelector(distance_calc, memory_monitor, logger)
            >>> excluded_snps = {"rs999"}  # Exclude specific SNP
            >>> first_snp = selector.select_first_snp(vcf_data, excluded=excluded_snps)
            >>> print(f"Selected first SNP: {first_snp}")
            Selected first SNP: rs123
        """
        excluded = excluded or set()
        candidates = []

        for sid, _ in vcf_data.snp_genotypes.items():
            if sid in excluded:
                continue
            maf = vcf_data.snp_maf_cache.get(sid, 0.0)
            candidates.append((sid, maf))

        if not candidates:
            raise BuildError("No SNPs available for selection after exclusions.")

        # Sort by MAF desc, SNP ID asc
        candidates.sort(key=lambda x: (-x[1], x[0]))
        return candidates[0][0]

    def build_discriminatory_set(
        self,
        vcf_data: VCFData,
        min_distance: int,
        excluded: Optional[Set[str]] = None,
    ) -> List[str]:
        """Greedily build initial SNP set maximizing MAF under exclusions.

        Args:
            vcf_data: VCF data containing SNP information
            min_distance: Minimum Hamming distance required between samples
            excluded: Set of SNP IDs to exclude from selection
            log_start: Whether to log start of primary selection
            log_distances: Whether to log distance calculations

        Returns:
            List of SNP IDs forming primary set

        Example:
            >>> selector = SNPSelector(distance_calc, memory_monitor, logger)
            >>> primary_set = selector.build_discriminatory_set(
            ...     vcf_data, min_distance=3, excluded={"rs999"}
            ... )
            >>> print(f"Primary set contains {len(primary_set)} SNPs")
            >>> print(f"SNPs: {', '.join(primary_set)}")
            Primary set contains 5 SNPs
            SNPs: rs123, rs456, rs789, rs101, rs202
        """
        excluded = excluded or set()
        current_snp_set: List[str] = []
        current_snps_geno: List[List[float]] = []

        # Choose first SNP
        current_snp = self.select_first_discriminatory_snp(vcf_data, excluded=excluded)
        current_snps_geno.append(vcf_data.snp_genotypes[current_snp].genotypes)
        current_snp_set.append(current_snp)

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Discriminatory SNP selection...")

        # Identify primary set of SNPs deterministically with tie-breakers
        while (
            self.distance_calc.calc_min_distance(current_snps_geno, log=False)
            < min_distance
        ):
            parent_nodes_info = []
            for sid, _ in vcf_data.snp_genotypes.items():
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
        self,
        snp_set: List[str],
        vcf_data: VCFData,
        min_distance: int,
        log_start: bool = True,
    ) -> List[str]:
        """Greedy backward elimination preserving minimal distance constraint.

        Optimized to avoid repeated distance recomputation by precomputing per-SNP
        contributions to pairwise distances across the upper triangle and updating
        totals incrementally when SNPs are removed.

        Args:
            snp_set: List of SNP IDs to optimize
            vcf_data: VCF data containing SNP information
            min_distance: Minimum Hamming distance to maintain
            log_start: Whether to log start of elimination process

        Returns:
            Optimized list of SNP IDs with minimal distance preserved

        Example:
            >>> selector = SNPSelector(distance_calc, memory_monitor, logger)
            >>> primary_set = ["rs1", "rs2", "rs3", "rs4", "rs5"]
            >>> optimized_set = selector.deterministic_eliminate(
            ...     primary_set, vcf_data, min_distance=3
            ... )
            >>> print(f"Optimized from {len(primary_set)} to {len(optimized_set)} SNPs")
            >>> print(f"Removed: {set(primary_set) - set(optimized_set)}")
            Optimized from 5 to 3 SNPs
            Removed: {'rs2', 'rs4'}
        """
        if log_start and self.logger.isEnabledFor(logging.INFO):
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
        vcf_data: VCFData,
        min_distance: int,
        n_sets: int,
    ) -> List[List[str]]:
        """Find all perfectly orthogonal (disjoint) minimal discriminating SNP sets.

        This implementation enumerates disjoint minimal sets across the entire VCF.
        The n_sets parameter is ignored: all disjoint sets found are returned.
        """
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Searching for disjoint minimal SNP sets...")

        # Enumerate orthogonal sets using global exclusions to enforce disjointness
        disjoint_sets: List[List[str]] = []
        excluded_global: Set[str] = set()

        alt_index = 0
        while n_sets > len(disjoint_sets):
            try:
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        f"Attempting to build discriminatory set #{alt_index + 1}"
                    )
                selected_set = self.build_discriminatory_set(
                    vcf_data,
                    min_distance,
                    excluded=set(excluded_global),
                )
            except BuildError:
                sys.exit(
                    f"Not enough polymorphic SNP to select #{alt_index + 1} discriminatory set. Exit."
                )
            selected_set_minimal = self.deterministic_eliminate(
                selected_set, vcf_data, min_distance, log_start=False
            )

            if not selected_set_minimal:
                break
            # Ensure strict disjointness from all previously selected sets
            if all(
                len(set(selected_set_minimal).intersection(set(s))) == 0
                for s in disjoint_sets
            ):
                disjoint_sets.append(sorted(selected_set_minimal))
                excluded_global.update(selected_set_minimal)
                alt_index += 1
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        f"Added alternative set #{alt_index} with {len(selected_set_minimal)} SNP(s)"
                    )
            else:
                excluded_global.update(selected_set_minimal)
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        "Candidate alternative set overlapped with existing; skipping"
                    )
                continue

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"{len(disjoint_sets)} orthogonal discriminating SNP set(s) selected."
            )
        return disjoint_sets
