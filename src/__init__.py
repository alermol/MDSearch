"""MDSearch - Minimum Discriminatory SNPs set Search.

A Python package for identifying minimal sets of Single Nucleotide Polymorphisms (SNPs)
that can discriminate between samples in a VCF file based on a specified minimum Hamming distance.
"""

from .app import MDSearchApp, MDSearchConfig
from .core.vcf_parser import VCFParser, VCFData, SNPData
from .core.distance_calculator import DistanceCalculator
from .core.snp_selector import SNPSelector
from .utils.memory_monitor import MemoryMonitor

__version__ = "1.0.0"
__all__ = [
    "MDSearchApp",
    "MDSearchConfig",
    "VCFParser",
    "VCFData",
    "SNPData",
    "DistanceCalculator",
    "SNPSelector",
    "MemoryMonitor",
]
