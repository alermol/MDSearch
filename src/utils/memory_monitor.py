"""Memory usage monitoring and warning system."""

import logging
import gc
from typing import Dict

import psutil

__all__ = ["MemoryMonitor"]


class MemoryMonitor:
    """Utility class for monitoring memory usage and providing warnings."""

    # Dynamic threshold percentages of total system memory
    WARNING_THRESHOLD_PERCENT = 50.0  # 50% of total system memory
    CRITICAL_THRESHOLD_PERCENT = 90.0  # 90% of total system memory

    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.process = psutil.Process()

        # Calculate dynamic thresholds based on total system memory
        total_memory_mb = psutil.virtual_memory().total / 1024 / 1024
        self.warning_threshold_mb = total_memory_mb * (
            self.WARNING_THRESHOLD_PERCENT / 100
        )
        self.critical_threshold_mb = total_memory_mb * (
            self.CRITICAL_THRESHOLD_PERCENT / 100
        )

        # Log the calculated thresholds for transparency
        self.logger.debug(
            f"Memory thresholds calculated: Warning={self.warning_threshold_mb:.1f}MB "
            f"({self.WARNING_THRESHOLD_PERCENT}%), Critical={self.critical_threshold_mb:.1f}MB "
            f"({self.CRITICAL_THRESHOLD_PERCENT}%) of {total_memory_mb:.1f}MB total"
        )

    def get_memory_usage_mb(self) -> float:
        """Get current memory usage in MB."""
        rss: int = self.process.memory_info().rss
        return float(rss / 1024 / 1024)

    def get_available_memory_mb(self) -> float:
        """Get available system memory in MB."""
        available: int = psutil.virtual_memory().available
        return float(available / 1024 / 1024)

    def check_memory_and_warn(self, operation: str = "operation") -> None:
        """Check current memory usage and warn if approaching limits."""
        current_mb = self.get_memory_usage_mb()
        available_mb = self.get_available_memory_mb()

        if current_mb > self.critical_threshold_mb:
            self.logger.warning(
                f"CRITICAL: High memory usage during {operation}: {current_mb:.1f}MB "
                f"(>{self.critical_threshold_mb:.1f}MB threshold). "
                f"Available: {available_mb:.1f}MB. Consider using smaller datasets "
                "or increase system memory."
            )
        elif current_mb > self.warning_threshold_mb:
            self.logger.warning(
                f"WARNING: Elevated memory usage during {operation}: {current_mb:.1f}MB "
                f"(>{self.warning_threshold_mb:.1f}MB threshold). "
                f"Available: {available_mb:.1f}MB. Monitor for potential issues."
            )
        else:
            self.logger.debug(f"Memory usage during {operation}: {current_mb:.1f}MB. ")

    def estimate_matrix_memory_mb(self, num_snps: int, num_samples: int) -> float:
        """Estimate memory usage for genotype matrix in MB."""
        # NumPy float64 = 8 bytes per element
        matrix_bytes = num_snps * num_samples * 8
        return matrix_bytes / 1024 / 1024

    def get_threshold_info(self) -> Dict[str, float]:
        """Get current memory threshold information."""
        return {
            "warning_threshold_mb": self.warning_threshold_mb,
            "critical_threshold_mb": self.critical_threshold_mb,
            "warning_percent": self.WARNING_THRESHOLD_PERCENT,
            "critical_percent": self.CRITICAL_THRESHOLD_PERCENT,
            "total_memory_mb": psutil.virtual_memory().total / 1024 / 1024,
        }

    def warn_for_large_dataset(self, num_snps: int, num_samples: int) -> None:
        """Warn user about potential memory issues with large datasets."""
        estimated_mb = self.estimate_matrix_memory_mb(num_snps, num_samples)
        available_mb = self.get_available_memory_mb()

        # Use dynamic thresholds for warnings
        if estimated_mb > available_mb * 0.8:  # Using 80% of available memory
            self.logger.warning(
                f"MEMORY WARNING: Dataset ({num_snps} SNPs × {num_samples} samples) "
                f"may require ~{estimated_mb:.1f}MB memory, but only {available_mb:.1f}MB available. "
                "Consider reducing dataset size or using a machine with more memory."
            )
        elif estimated_mb > self.warning_threshold_mb:  # Use dynamic warning threshold
            self.logger.warning(
                f"Large dataset detected ({num_snps} SNPs × {num_samples} samples). "
                f"Estimated memory usage ~{estimated_mb:.1f}MB exceeds warning threshold "
                f"({self.warning_threshold_mb:.1f}MB). Monitor memory usage carefully."
            )
        elif (
            estimated_mb > self.warning_threshold_mb * 0.5
        ):  # Half of warning threshold
            self.logger.info(
                f"Large dataset detected ({num_snps} SNPs × {num_samples} samples). "
                f"Estimated memory usage: ~{estimated_mb:.1f}MB"
            )

    def force_garbage_collection(self) -> None:
        """Force garbage collection to free up memory."""
        gc.collect()
        self.logger.debug("Forced garbage collection completed")
