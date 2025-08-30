"""Hamming distance calculation optimized for large datasets."""

import logging
from typing import List, Protocol, Iterable, Optional, Callable

import numpy as np
from numpy.typing import NDArray
from tqdm import tqdm

from ..utils.memory_monitor import MemoryMonitor
from .vcf_parser import SNPData

__all__ = ["DistanceCalculator", "SNPDataMapping"]


class SNPDataMapping(Protocol):
    """Protocol for objects that provide SNP data access like dict[str, SNPData]."""

    def __getitem__(self, key: str) -> SNPData: ...


class DistanceCalculator:
    """Optimized Hamming distance calculations."""

    def __init__(
        self,
        memory_monitor: MemoryMonitor,
        logger: logging.Logger,
        show_progress: bool = True,
        shutdown_checker: Optional[Callable[[], bool]] = None,
    ):
        """Initialize distance calculator with memory monitoring and logging.

        Args:
            memory_monitor: MemoryMonitor instance for tracking memory usage
            logger: Logger instance for output
            show_progress: Whether to show progress bars for large datasets
            shutdown_checker: Optional function to check if shutdown was requested

        Example:
            >>> from src.utils.memory_monitor import MemoryMonitor
            >>> from src.utils.logging_setup import setup_logger
            >>> logger = setup_logger("distance_calc")
            >>> memory_monitor = MemoryMonitor(logger)
            >>> calc = DistanceCalculator(memory_monitor, logger, show_progress=True)
        """
        self.memory_monitor = memory_monitor
        self.logger = logger
        self.show_progress = show_progress
        self.shutdown_checker = shutdown_checker

    def calc_min_distance(self, snps: List[List[float]], *, log: bool = True) -> float:
        """Return minimal pairwise Hamming distance across samples for given SNPs.

        Optimized: compute only upper-triangle pairwise distances (j > i) and avoid
        self-pairs. Uses vectorized comparisons per anchor column.

        Args:
            snps: List of SNP genotype lists, each containing sample genotypes
            log: Whether to log distance statistics

        Returns:
            Minimum Hamming distance between any pair of samples

        Example:
            >>> calc = DistanceCalculator(memory_monitor, logger)
            >>> # Example with 3 SNPs and 4 samples
            >>> snp_genotypes = [
            ...     [0.0, 1.0, 0.0, 1.0],  # SNP 1: samples 0,1,2,3
            ...     [1.0, 1.0, 0.0, 0.0],  # SNP 2: samples 0,1,2,3
            ...     [0.0, 0.0, 1.0, 1.0],  # SNP 3: samples 0,1,2,3
            ... ]
            >>> min_dist = calc.calc_min_distance(snp_genotypes)
            >>> print(f"Minimum distance between samples: {min_dist}")
            Minimum distance between samples: 2.0
        """
        if log and self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Calculate pairwise distance based on {len(snps)} SNPs..."
            )

        # Monitor memory before creating large arrays
        self.memory_monitor.check_memory_and_warn("distance calculation start")

        snps_array = np.array([i for i in snps])  # shape: (num_snps, num_samples)
        num_samples = snps_array.shape[1]
        pairwise_distances: List[float] = []

        if num_samples <= 1:
            if log and self.logger.isEnabledFor(logging.INFO):
                self.logger.info("Distance between samples (min/med/avg/max): 0/0/0/0")
            return 0.0

        # Create progress bar for sample comparisons
        sample_range: Iterable[int] = range(num_samples - 1)
        if (
            self.show_progress and num_samples > 10 and log
        ):  # Only show for larger datasets
            sample_range = tqdm(
                sample_range,
                desc="Computing pairwise distances",
                unit="sample",
                leave=False,
            )

        for i in sample_range:
            # Check for graceful shutdown request
            if self.shutdown_checker and self.shutdown_checker():
                if log and self.logger.isEnabledFor(logging.INFO):
                    self.logger.info(
                        "Graceful shutdown requested. Stopping distance calculation."
                    )
                break

            col_i = snps_array[:, i][:, None]  # (num_snps, 1)
            rest = snps_array[:, i + 1 :]  # (num_snps, num_samples - i - 1)
            if rest.shape[1] == 0:
                continue
            valid_mask = (~np.isnan(col_i)) & (~np.isnan(rest))
            diffs = (col_i != rest) & valid_mask
            dists = diffs.sum(axis=0).astype(float)  # per pair distances vs column i
            pairwise_distances.extend(dists.tolist())

        # Handle early exit case
        if not pairwise_distances:
            if log and self.logger.isEnabledFor(logging.INFO):
                self.logger.info(
                    "Distance calculation incomplete due to shutdown request"
                )
            return float("inf")  # Return infinity to indicate incomplete calculation

        res: NDArray[np.float64] = np.array(pairwise_distances, dtype=float)

        if log and self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                "Distance between samples (min/med/avg/max): "
                f"{np.min(res)}/{np.median(res)}/{round(float(np.mean(res)), 1)}/{np.max(res)}"
            )

        # Monitor memory after distance calculation
        self.memory_monitor.check_memory_and_warn("distance calculation complete")

        return float(np.min(res))

    def calc_distance_for_snp_ids(
        self, snp_ids: List[str], snp_data: SNPDataMapping, *, log: bool = True
    ) -> float:
        """Helper computing minimal distance for a list of SNP IDs.

        Args:
            snp_ids: List of SNP identifiers to include in distance calculation
            snp_data: Mapping from SNP ID to SNPData
            log: Whether to log distance calculation details

        Returns:
            Minimum Hamming distance between any pair of samples

        Example:
            >>> calc = DistanceCalculator(memory_monitor, logger)
            >>> selected_snps = ["rs123", "rs456", "rs789"]
            >>> min_distance = calc.calc_distance_for_snp_ids(
            ...     selected_snps, vcf_data.snp_genotypes
            ... )
            >>> print(f"Minimum distance with selected SNPs: {min_distance}")
            Minimum distance with selected SNPs: 3.0
        """
        if not snp_ids:
            return 0.0
        genos = [snp_data[sid].genotypes for sid in snp_ids]
        return self.calc_min_distance(genos, log=log)
