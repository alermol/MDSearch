"""Hamming distance calculation optimized for large datasets."""

import logging
from typing import List, Dict

import numpy as np
from numpy.typing import NDArray

from ..utils.memory_monitor import MemoryMonitor
from .vcf_parser import SNPData


class DistanceCalculator:
    """Optimized Hamming distance calculations."""
    
    def __init__(self, memory_monitor: MemoryMonitor, logger: logging.Logger):
        self.memory_monitor = memory_monitor
        self.logger = logger
    
    def calc_min_distance(self, snps: List[List[float]]) -> float:
        """Return minimal pairwise Hamming distance across samples for given SNPs.

        Optimized: compute only upper-triangle pairwise distances (j > i) and avoid
        self-pairs. Uses vectorized comparisons per anchor column.
        """
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Calculate pairwise distance based on {len(snps)} SNPs..."
            )

        # Monitor memory before creating large arrays
        self.memory_monitor.check_memory_and_warn("distance calculation start")

        snps_array = np.array([i for i in snps])  # shape: (num_snps, num_samples)
        num_samples = snps_array.shape[1]
        pairwise_distances: List[float] = []
        
        if num_samples <= 1:
            if self.logger.isEnabledFor(logging.INFO):
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
            
        res: NDArray[np.float64] = np.array(pairwise_distances, dtype=float)
        
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                "Distance between samples (min/med/avg/max): "
                f"{np.min(res)}/{np.median(res)}/{round(float(np.mean(res)), 1)}/{np.max(res)}"
            )

        # Monitor memory after distance calculation
        self.memory_monitor.check_memory_and_warn("distance calculation complete")

        return float(np.min(res))
    
    def calc_distance_for_snp_ids(self, snp_ids: List[str], snp_data: Dict[str, SNPData]) -> float:
        """Helper computing minimal distance for a list of SNP IDs."""
        if not snp_ids:
            return 0.0
        genos = [snp_data[sid].genotypes for sid in snp_ids]
        return self.calc_min_distance(genos)
