"""SNP selection and optimization algorithms."""

import logging
from typing import List, Set, Optional, Union, Callable

import numpy as np
from numpy.typing import NDArray

from .distance_calculator import DistanceCalculator
from .vcf_parser import VCFData
from .genotype_utils import calculate_snp_entropy_score
from ..utils.memory_monitor import MemoryMonitor

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
        weight_entropy: float = 0.5,
        weight_maf: float = 0.5,
    ):
        """Initialize SNP selector with required components."""
        self.distance_calc = distance_calculator
        self.memory_monitor = memory_monitor
        self.logger = logger
        self.shutdown_checker = shutdown_checker
        self.weight_entropy = weight_entropy
        self.weight_maf = weight_maf

    def select_first_discriminatory_snp(
        self, vcf_data: VCFData, excluded: Optional[Set[str]] = None
    ) -> str:
        """Select SNP with highest entropy score not in excluded; tie-break by SNP ID.
        
        Args:
            vcf_data: VCF data containing SNP genotypes and cached scores
            excluded: Set of SNP IDs to exclude from selection
            
        Returns:
            SNP ID with highest entropy score
            
        Raises:
            BuildError: If no SNPs available after exclusions
        """
        excluded = excluded or set()
        candidates = []

        for sid, snp_data in vcf_data.snp_genotypes.items():
            if sid in excluded:
                continue
            maf = vcf_data.snp_maf_cache.get(sid, 0.0)
            entropy = vcf_data.snp_entropy_cache.get(sid, 0.0)
            entropy_score = calculate_snp_entropy_score(
                snp_data.genotypes,
                maf,
                self.weight_entropy,
                self.weight_maf,
                entropy=entropy,
            )
            candidates.append((sid, entropy_score))

        if not candidates:
            raise BuildError("No SNPs available for selection after exclusions.")

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
            vcf_data: VCF data containing SNP genotypes and cached scores
            min_distance: Minimum Hamming distance requirement between samples
            excluded: Set of SNP IDs to exclude from selection
            
        Returns:
            List of SNP IDs forming the discriminatory set
            
        Raises:
            BuildError: If not enough SNPs available to meet distance requirement
        """
        excluded = excluded or set()
        current_snp_set: List[str] = []
        current_snps_geno: List[List[float]] = []

        current_snp = self.select_first_discriminatory_snp(vcf_data, excluded=excluded)
        current_snps_geno.append(vcf_data.snp_genotypes[current_snp].genotypes)
        current_snp_set.append(current_snp)

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Discriminatory SNP selection...")

        while (
            self.distance_calc.calc_min_distance(current_snps_geno, log=False)
            < min_distance
        ):
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
                entropy_score = calculate_snp_entropy_score(
                    snp_data.genotypes,
                    maf,
                    self.weight_entropy,
                    self.weight_maf,
                    entropy=entropy,
                )
                parent_nodes_info.append((sid, entropy_score))

            if not parent_nodes_info:
                raise BuildError("Not enough polymorphic SNP to discriminate samples.")

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
        
        Args:
            vcf_data: VCF data containing SNP genotypes and cached scores
            excluded: Set of SNP IDs to exclude from scoring
            
        Returns:
            List of tuples: (snp_id, entropy_score, maf, raw_entropy)
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
        
        Args:
            snp_set: Initial list of SNP IDs to optimize
            vcf_data: VCF data containing SNP genotypes
            min_distance: Minimum Hamming distance requirement to maintain
            log_start: Whether to log the start of elimination process
            
        Returns:
            Optimized list of SNP IDs maintaining minimum distance
        """
        if log_start and self.logger.isEnabledFor(logging.INFO):
            self.logger.info("Backward one-by-one elimination...")
            self.logger.info(
                f"Starting with {len(snp_set)} SNPs, target distance: {min_distance}"
            )

        self.memory_monitor.check_memory_and_warn("elimination algorithm start")

        optimized_ids: List[str] = list(snp_set)
        snps_matrix = np.array(
            [vcf_data.snp_genotypes[sid].genotypes for sid in optimized_ids]
        )
        num_snps, num_samples = snps_matrix.shape

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug(
                f"Elimination matrix: {num_snps} SNPs × {num_samples} samples"
            )
            self.logger.debug(f"Initial SNP set: {optimized_ids}")

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
            contrib_matrix = np.vstack(pair_contribs)
            current_pair_dists = contrib_matrix.sum(axis=0).astype(np.int32)
        else:
            contrib_matrix = np.zeros((0, 0), dtype=np.int32)
            current_pair_dists = np.zeros(0, dtype=np.int32)

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug(f"Contribution matrix shape: {contrib_matrix.shape}")
            self.logger.debug(
                f"Initial pairwise distances: min={current_pair_dists.min() if current_pair_dists.size > 0 else 'N/A'}"
            )

        removed_count = 0
        for sid in sorted(list(optimized_ids)):
            if sid not in optimized_ids:
                continue
            idx = optimized_ids.index(sid)
            candidate_pair_dists = current_pair_dists - contrib_matrix[idx]
            if (
                candidate_pair_dists.size == 0
                or candidate_pair_dists.min() >= min_distance
            ):
                optimized_ids.pop(idx)
                current_pair_dists = candidate_pair_dists
                if contrib_matrix.size:
                    contrib_matrix = np.delete(contrib_matrix, idx, axis=0)

                removed_count += 1
                if self.logger.isEnabledFor(logging.DEBUG):
                    self.logger.debug(
                        f"Removed SNP: {sid} (remaining: {len(optimized_ids)})"
                    )
                    if candidate_pair_dists.size > 0:
                        self.logger.debug(
                            f"New min distance: {candidate_pair_dists.min()}"
                        )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Elimination complete: removed {removed_count} SNPs, final set: {len(optimized_ids)} SNPs"
            )

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
        
        Args:
            vcf_data: VCF data containing SNP genotypes and cached scores
            min_distance: Minimum Hamming distance requirement between samples
            n_sets: Number of sets to find (0 or "unlimited" for all possible)
            
        Returns:
            List of SNP sets, each containing SNP IDs for one discriminatory set
        """
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

        disjoint_sets: List[List[str]] = []
        excluded_global: Set[str] = set()

        alt_index = 0

        while n_sets_int == 0 or n_sets_int > len(disjoint_sets):
            if self.shutdown_checker and self.shutdown_checker():
                if self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        "Graceful shutdown requested. Stopping SNP set search."
                    )
                break

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
