"""Utility modules for infrastructure and helpers."""

from .memory_monitor import MemoryMonitor
from .logging_setup import setup_logger
from .validation import validate_cli_arguments, validate_overlap_constraints

__all__ = [
    "MemoryMonitor",
    "setup_logger",
    "validate_cli_arguments",
    "validate_overlap_constraints",
    "ensure_variant_index",
]


from typing import Optional
from pathlib import Path
import logging
import importlib


def ensure_variant_index(
    vpath: Path, fmt_letter: Optional[str], logger: Optional[logging.Logger]
) -> None:
    # Dynamic import to satisfy static analyzers while keeping runtime resolution
    module = importlib.import_module(f"{__package__}.indexing")
    module.ensure_variant_index(vpath, fmt_letter, logger)
