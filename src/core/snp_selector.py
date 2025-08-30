"""SNP selection and optimization algorithms."""

import logging
from typing import List, Set, Optional, Union, Callable

import numpy as np
from numpy.typing import NDArray

from .distance_calculator import DistanceCalculator
from .vcf_parser import VCFData
from .genotype_utils import calculate_snp_entropy_score
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
        shutdown_checker: Optional[Callable[[], bool]] = None,
    ):
        """Initialize SNP selector with required components.

        Args:
            distance_calculator: DistanceCalculator instance for computing distances
            memory_monitor: MemoryMonitor instance for tracking memory usage
            logger: Logger instance for output
            shutdown_checker: Optional function to check if shutdown was requested

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
        self.shutdown_checker = shutdown_checker

    def select_first_discriminatory_snp(
        self, vcf_data: VCFData, excluded: Optional[Set[str]] = None
    ) -> str:
        """Select SNP with highest entropy score not in excluded; tie-break by SNP ID.

        Args:
            vcf_data: VCF data containing SNP information
            excluded: Set of SNP IDs to exclude from selection

        Returns:
            SNP ID with highest entropy score (or lexicographically first if tied)

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

        for sid, snp_data in vcf_data.snp_genotypes.items():
            if sid in excluded:
                continue
            maf = vcf_data.snp_maf_cache.get(sid, 0.0)
            entropy = vcf_data.snp_entropy_cache.get(sid, 0.0)
            # Use cached entropy for scoring
            entropy_score = calculate_snp_entropy_score(
                snp_data.genotypes, maf, entropy=entropy
            )
            candidates.append((sid, entropy_score))

        if not candidates:
            raise BuildError("No SNPs available for selection after exclusions.")

        # Sort by entropy score desc, SNP ID asc
        candidates.sort(key=lambda x: (-x[1], x[0]))
        return candidates[0][0]

    def build_discriminatory_set(
        self,
        vcf_data: VCFData,
        min_distance: int,
        excluded: Optional[Set[str]] = None,
    ) -> List[str]:
        """Greedily build initial SNP set maximizing entropy score under exclusions.

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
            # Check for graceful shutdown request
            if self.shutdown_checker and self.shutdown_checker():
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        "Graceful shutdown requested. Stopping discriminatory set building."
                    )
                break

            parent_nodes_info = []
            for sid, snp_data in vcf_data.snp_genotypes.items():
                if (sid in current_snp_set) or (sid in excluded):
                    continue
                maf = vcf_data.snp_maf_cache.get(sid, 0.0)
                entropy = vcf_data.snp_entropy_cache.get(sid, 0.0)
                # Use cached entropy for scoring
                entropy_score = calculate_snp_entropy_score(
                    snp_data.genotypes, maf, entropy=entropy
                )
                parent_nodes_info.append((sid, entropy_score))

            if not parent_nodes_info:
                raise BuildError("Not enough polymorphic SNP to discriminate samples.")

            # Sort by entropy score desc, SNP ID asc
            parent_nodes_info.sort(key=lambda x: (-x[1], x[0]))
            current_snp = parent_nodes_info[0][0]
            current_snps_geno.append(vcf_data.snp_genotypes[current_snp].genotypes)
            current_snp_set.append(current_snp)

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"After 1st step {len(current_snp_set)} primary SNP selected"
            )
        return current_snp_set

    def get_snp_entropy_scores(
        self, vcf_data: VCFData, excluded: Optional[Set[str]] = None
    ) -> List[tuple[str, float, float, float]]:
        """Get entropy scores for all SNPs in the dataset.

        This method calculates and returns entropy scores for all SNPs, which can be
        useful for analysis, debugging, and understanding the information content
        distribution across the dataset.

        Args:
            vcf_data: VCF data containing SNP information
            excluded: Set of SNP IDs to exclude from scoring

        Returns:
            List of tuples containing (SNP_ID, entropy_score, maf, raw_entropy)

        Example:
            >>> selector = SNPSelector(distance_calc, memory_monitor, logger)
            >>> scores = selector.get_snp_entropy_scores(vcf_data, excluded={"rs999"})
            >>> for snp_id, score, maf, entropy in sorted(scores, key=lambda x: -x[1])[:5]:
            ...     print(f"{snp_id}: score={score:.3f}, maf={maf:.3f}, entropy={entropy:.3f}")
            rs123: score=0.910, maf=0.400, entropy=1.585
            rs456: score=0.863, maf=0.350, entropy=1.500
        """
        excluded = excluded or set()
        scores = []

        for sid, snp_data in vcf_data.snp_genotypes.items():
            if sid in excluded:
                continue

            maf = vcf_data.snp_maf_cache.get(sid, 0.0)
            raw_entropy = vcf_data.snp_entropy_cache.get(sid, 0.0)
            entropy_score = calculate_snp_entropy_score(
                snp_data.genotypes, maf, entropy=raw_entropy
            )

            scores.append((sid, entropy_score, maf, raw_entropy))

        return scores

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
        n_sets: Union[int, str],
    ) -> List[List[str]]:
        """Find all perfectly orthogonal (disjoint) minimal discriminating SNP sets.

        This implementation enumerates disjoint minimal sets across the entire VCF.
        When n_sets is 0 (unlimited mode), it finds all possible disjoint sets.
        When n_sets > 0, it finds up to that many sets.
        """
        # Convert n_sets to int for internal use
        n_sets_int = int(n_sets) if isinstance(n_sets, str) else n_sets

        if self.logger.isEnabledFor(logging.INFO):
            if n_sets_int == 0:
                self.logger.info(
                    "Searching for all possible disjoint minimal SNP sets (unlimited mode)..."
                )
            else:
                self.logger.info(
                    f"Searching for up to {n_sets_int} disjoint minimal SNP sets..."
                )

        # Enumerate orthogonal sets using global exclusions to enforce disjointness
        disjoint_sets: List[List[str]] = []
        excluded_global: Set[str] = set()

        alt_index = 0

        # Continue searching until either we reach n_sets (if specified) or no more sets can be found
        while n_sets_int == 0 or n_sets_int > len(disjoint_sets):
            # Check for graceful shutdown request
            if self.shutdown_checker and self.shutdown_checker():
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        "Graceful shutdown requested. Stopping SNP set search."
                    )
                break

            # Report available SNPs for this set generation
            available_snps = len(vcf_data.snp_genotypes) - len(excluded_global)
            if self.logger.isEnabledFor(logging.INFO):
                self.logger.info(
                    f"Attempting to build discriminatory set #{alt_index + 1} "
                    f"(SNPs available: {available_snps})"
                )
            try:
                selected_set = self.build_discriminatory_set(
                    vcf_data,
                    min_distance,
                    excluded=set(excluded_global),
                )
            except BuildError:
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info("No more disjoint SNP sets can be found")
                break
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
                remaining_snps = len(vcf_data.snp_genotypes) - len(excluded_global)
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        f"Added alternative set #{alt_index} with {len(selected_set_minimal)} SNP(s) "
                        f"(SNPs remaining: {remaining_snps})"
                    )
            else:
                excluded_global.update(selected_set_minimal)
                remaining_snps = len(vcf_data.snp_genotypes) - len(excluded_global)
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        f"Candidate alternative set overlapped with existing; skipping "
                        f"(SNPs remaining: {remaining_snps})"
                    )
                continue

        if self.logger.isEnabledFor(logging.INFO):
            if n_sets_int == 0:
                self.logger.info(
                    f"Found {len(disjoint_sets)} orthogonal discriminating SNP set(s) in unlimited mode."
                )
            else:
                self.logger.info(
                    f"{len(disjoint_sets)} orthogonal discriminating SNP set(s) selected."
                )
        return disjoint_sets
