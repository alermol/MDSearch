"""Utility modules for infrastructure and helpers."""

from .memory_monitor import MemoryMonitor
from .logging_setup import setup_logger
from .validation import validate_cli_arguments

__all__ = [
    "MemoryMonitor",
    "setup_logger",
    "validate_cli_arguments",
    "ensure_variant_index",
]


from typing import Optional
from pathlib import Path
import logging
import importlib


def ensure_variant_index(
    vpath: Path, fmt_letter: Optional[str], logger: Optional[logging.Logger]
) -> None:
    """Ensure an index exists for the given variant file.
    
    This is a convenience function that imports and calls the actual implementation
    from the indexing module to avoid circular imports.
    
    Args:
        vpath: Path to the variant file (VCF/BCF)
        fmt_letter: Format identifier (v, z, u, b) or None for auto-detection
        logger: Optional logger instance for progress messages
    """
    module = importlib.import_module(f"{__package__}.indexing")
    module.ensure_variant_index(vpath, fmt_letter, logger)
