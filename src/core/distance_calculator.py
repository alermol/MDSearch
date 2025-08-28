"""Hamming distance calculation optimized for large datasets."""

import logging
from typing import List, Protocol, Iterable

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
    ):
        self.memory_monitor = memory_monitor
        self.logger = logger
        self.show_progress = show_progress

    def calc_min_distance(self, snps: List[List[float]], *, log: bool = True) -> float:
        """Return minimal pairwise Hamming distance across samples for given SNPs.

        Optimized: compute only upper-triangle pairwise distances (j > i) and avoid
        self-pairs. Uses vectorized comparisons per anchor column.
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
        if self.show_progress and num_samples > 10 and log:  # Only show for larger datasets
            sample_range = tqdm(
                sample_range,
                desc="Computing pairwise distances",
                unit="sample",
                leave=False,
            )

        for i in sample_range:
            col_i = snps_array[:, i][:, None]  # (num_snps, 1)
            rest = snps_array[:, i + 1 :]  # (num_snps, num_samples - i - 1)
            if rest.shape[1] == 0:
                continue
            valid_mask = (~np.isnan(col_i)) & (~np.isnan(rest))
            diffs = (col_i != rest) & valid_mask
            dists = diffs.sum(axis=0).astype(float)  # per pair distances vs column i
            pairwise_distances.extend(dists.tolist())

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
        """Helper computing minimal distance for a list of SNP IDs."""
        if not snp_ids:
            return 0.0
        genos = [snp_data[sid].genotypes for sid in snp_ids]
        return self.calc_min_distance(genos, log=log)
