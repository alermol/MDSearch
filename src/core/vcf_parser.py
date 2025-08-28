"""VCF file parsing and validation."""

import sys
import logging
from pathlib import Path
from typing import List, Dict, Set, Any, Tuple, Optional, Iterator
from dataclasses import dataclass, field
from collections import OrderedDict

from ..utils.memory_monitor import MemoryMonitor
from .genotype_utils import extract_gt, gt_to_value, calculate_maf

__all__ = [
    "VCFHeaders",
    "SNPData",
    "VCFData",
    "SNPMetadata",
    "LazyVCFData",
    "VCFParser",
]


@dataclass
class VCFHeaders:
    """VCF header information."""

    samples: List[str]
    fileformat_found: bool


@dataclass
class SNPData:
    """Data for a single SNP."""

    genotypes: List[float]
    sample_fields: List[str]
    maf: float


@dataclass
class VCFData:
    """Complete VCF data structure."""

    headers: VCFHeaders
    snp_genotypes: Dict[str, SNPData]
    snp_maf_cache: Dict[str, float]


@dataclass
class SNPMetadata:
    """Lightweight SNP metadata without genotype data for lazy loading."""

    snp_id: str
    maf: float
    file_position: int  # Byte position where this SNP line starts
    line_length: int  # Length of the line for seeking
    line_number: int  # Line number in file for error reporting


class LRUCache:
    """Simple LRU cache for SNP data."""

    def __init__(self, max_size: int = 1000):
        self.max_size = max_size
        self.cache: OrderedDict[str, SNPData] = OrderedDict()

    def get(self, key: str) -> Optional[SNPData]:
        """Get item from cache and move it to end (most recently used)."""
        if key in self.cache:
            # Move to end (most recently used)
            self.cache.move_to_end(key)
            return self.cache[key]
        return None

    def put(self, key: str, value: SNPData) -> None:
        """Add item to cache, evicting least recently used if needed."""
        if key in self.cache:
            # Update existing item and move to end
            self.cache[key] = value
            self.cache.move_to_end(key)
        else:
            # Add new item
            self.cache[key] = value
            # Evict least recently used if over capacity
            if len(self.cache) > self.max_size:
                self.cache.popitem(last=False)  # Remove first (least recently used)

    def clear(self) -> None:
        """Clear all cached data."""
        self.cache.clear()

    def size(self) -> int:
        """Return current cache size."""
        return len(self.cache)


@dataclass
class LazyVCFData:
    """VCF data with lazy loading support."""

    headers: VCFHeaders
    snp_metadata: Dict[str, SNPMetadata]
    snp_maf_cache: Dict[str, float]
    _file_path: Path
    _ploidy: int
    _convert_het: bool
    _cache: LRUCache = field(default_factory=lambda: LRUCache(max_size=1000))
    _memory_monitor: Optional[MemoryMonitor] = None
    _logger: Optional[logging.Logger] = None

    def get_snp_data(self, snp_id: str) -> SNPData:
        """Get SNP data with lazy loading and caching."""
        # Check cache first
        cached_data = self._cache.get(snp_id)
        if cached_data is not None:
            return cached_data

        # Load from file if not in cache
        if snp_id not in self.snp_metadata:
            raise KeyError(f"SNP '{snp_id}' not found in VCF data")

        metadata = self.snp_metadata[snp_id]
        snp_data = self._load_snp_from_file(metadata)

        # Cache the loaded data
        self._cache.put(snp_id, snp_data)

        if self._logger and self._logger.isEnabledFor(logging.DEBUG):
            self._logger.debug(
                f"Lazy loaded SNP '{snp_id}' (cache size: {self._cache.size()})"
            )

        return snp_data

    def _load_snp_from_file(self, metadata: SNPMetadata) -> SNPData:
        """Load a specific SNP from file using its metadata."""
        try:
            with open(self._file_path, "r") as vcf_file:
                # Seek to the SNP's position
                vcf_file.seek(metadata.file_position)
                line = vcf_file.read(metadata.line_length).strip()

                # Parse the line
                parts = line.split("\t")
                if len(parts) < 9:
                    raise ValueError(
                        f"Invalid VCF line format at position {metadata.file_position}"
                    )

                # Extract sample data
                format_field = parts[8]
                sample_fields = parts[9:]

                # Convert genotypes
                geno: List[float] = []
                for sample_field in sample_fields:
                    gt = extract_gt(format_field, sample_field)
                    geno.append(gt_to_value(gt, self._ploidy, self._convert_het))

                return SNPData(
                    genotypes=geno, sample_fields=sample_fields, maf=metadata.maf
                )

        except Exception as e:
            raise RuntimeError(
                f"Failed to load SNP '{metadata.snp_id}' from line {metadata.line_number}: {e}"
            ) from e

    @property
    def snp_genotypes(self) -> "LazyGenotypesDict":
        """Provide backward compatibility by creating a dict-like interface."""
        return LazyGenotypesDict(self)

    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache performance statistics."""
        return {
            "cache_size": self._cache.size(),
            "max_cache_size": self._cache.max_size,
            "total_snps": len(self.snp_metadata),
        }


class LazyGenotypesDict:
    """Dict-like interface for lazy loaded SNP genotypes."""

    def __init__(self, lazy_vcf_data: LazyVCFData):
        self._lazy_data = lazy_vcf_data

    def __getitem__(self, snp_id: str) -> SNPData:
        return self._lazy_data.get_snp_data(snp_id)

    def __contains__(self, snp_id: str) -> bool:
        return snp_id in self._lazy_data.snp_metadata

    def __iter__(self) -> Iterator[str]:
        return iter(self._lazy_data.snp_metadata.keys())

    def __len__(self) -> int:
        return len(self._lazy_data.snp_metadata)

    def items(self) -> Iterator[Tuple[str, SNPData]]:
        """Iterate over SNP ID and data pairs (loads data lazily)."""
        for snp_id in self._lazy_data.snp_metadata:
            yield snp_id, self._lazy_data.get_snp_data(snp_id)

    def keys(self) -> Iterator[str]:
        """Return SNP IDs."""
        return iter(self._lazy_data.snp_metadata.keys())

    def values(self) -> Iterator[SNPData]:
        """Return SNP data (loads all data - use carefully!)."""
        for snp_id in self._lazy_data.snp_metadata:
            yield self._lazy_data.get_snp_data(snp_id)

    def get(self, snp_id: str, default: Optional[SNPData] = None) -> Optional[SNPData]:
        """Get SNP data with optional default."""
        try:
            return self._lazy_data.get_snp_data(snp_id)
        except KeyError:
            return default


class VCFParser:
    """Handles VCF file parsing and validation."""

    def __init__(self, memory_monitor: MemoryMonitor, logger: logging.Logger):
        self.memory_monitor = memory_monitor
        self.logger = logger

    def parse_and_validate(
        self, vcf_path: Path, ploidy: int, convert_het: bool
    ) -> VCFData:
        """Parse VCF file and return structured data (loads all SNPs into memory)."""
        # First pass: validate headers and count samples
        headers = self._validate_headers(vcf_path)

        # Second pass: parse and validate data lines
        snp_genotypes, snp_maf_cache = self._validate_variants(
            vcf_path, headers, ploidy, convert_het
        )

        # Check memory usage after VCF parsing and warn for large datasets
        self.memory_monitor.check_memory_and_warn("VCF parsing")
        self.memory_monitor.warn_for_large_dataset(
            len(snp_genotypes), len(headers.samples)
        )

        return VCFData(
            headers=headers, snp_genotypes=snp_genotypes, snp_maf_cache=snp_maf_cache
        )

    def parse_and_validate_lazy(
        self, vcf_path: Path, ploidy: int, convert_het: bool, cache_size: int = 1000
    ) -> LazyVCFData:
        """Parse VCF file and return lazy-loaded data structure (minimal memory usage)."""
        # First pass: validate headers and count samples
        headers = self._validate_headers(vcf_path)

        # Second pass: build metadata without loading all genotype data
        snp_metadata, snp_maf_cache = self._build_snp_metadata(
            vcf_path, headers, ploidy, convert_het
        )

        # Create lazy VCF data structure
        lazy_data = LazyVCFData(
            headers=headers,
            snp_metadata=snp_metadata,
            snp_maf_cache=snp_maf_cache,
            _file_path=vcf_path,
            _ploidy=ploidy,
            _convert_het=convert_het,
            _cache=LRUCache(max_size=cache_size),
            _memory_monitor=self.memory_monitor,
            _logger=self.logger,
        )

        # Check memory usage after VCF metadata parsing
        self.memory_monitor.check_memory_and_warn("VCF metadata parsing")

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"Lazy VCF loading enabled: {len(snp_metadata)} variants, "
                f"{len(headers.samples)} samples, cache_size={cache_size}"
            )

        return lazy_data

    def _validate_headers(self, vcf_path: Path) -> VCFHeaders:
        """Validate VCF headers and extract sample information."""
        header_line = None
        samples = []
        fileformat_found = False
        line_number = 0

        with open(vcf_path) as vcf:
            for line in vcf:
                line_number += 1
                line = line.strip()

                if not line:
                    continue

                if line.startswith("##fileformat=VCF"):
                    fileformat_found = True
                elif line.startswith("#CHROM"):
                    header_line = line
                    parts = line.split("\t")
                    if len(parts) < 9:
                        sys.exit(
                            f"ERROR: Line {line_number}: Invalid header format. "
                            f"Expected at least 9 columns, found {len(parts)}."
                        )
                    samples = parts[9:]
                    break
                elif line.startswith("#"):
                    continue  # Skip other header lines
                else:
                    sys.exit(
                        f"ERROR: Line {line_number}: Found data line before required #CHROM header. "
                        "VCF must have proper header structure."
                    )

        # Validate headers
        if not fileformat_found:
            sys.exit(
                "ERROR: Missing required VCF header '##fileformat=VCFv4.x'. "
                "Ensure input is a valid VCF file."
            )

        if header_line is None:
            sys.exit(
                "ERROR: Missing required #CHROM header line. "
                "VCF must have column headers."
            )

        # Validate minimum sample count
        if len(samples) < 2:
            sys.exit(
                f"ERROR: Insufficient samples for distance calculation. "
                f"Found {len(samples)} samples, need at least 2."
            )

        return VCFHeaders(samples=samples, fileformat_found=fileformat_found)

    def _validate_variants(
        self, vcf_path: Path, headers: VCFHeaders, ploidy: int, convert_het: bool
    ) -> Tuple[Dict[str, SNPData], Dict[str, float]]:
        """Parse and validate variant lines."""
        seen_ids: Set[str] = set()
        data_line_count = 0
        snp_genotypes: Dict[str, SNPData] = {}
        snp_maf_cache: Dict[str, float] = {}

        with open(vcf_path) as vcf:
            line_number = 0
            for line in vcf:
                line_number += 1
                line = line.strip()

                if line.startswith("#") or not line:
                    continue

                data_line_count += 1
                parts = line.split("\t")

                # Validate line format
                expected_columns = 9 + len(headers.samples)
                if len(parts) != expected_columns:
                    sys.exit(
                        f"ERROR: Line {line_number}: Invalid number of columns. "
                        f"Expected {expected_columns}, found {len(parts)}."
                    )

                # Validate required columns exist
                if len(parts) < 5:
                    sys.exit(
                        f"ERROR: Line {line_number}: Missing required VCF columns "
                        "(CHROM, POS, ID, REF, ALT)."
                    )

                # Validate SNP ID
                snp_id = parts[2]
                if (not snp_id) or (snp_id == "."):
                    sys.exit(
                        f"ERROR: Line {line_number}: Missing or placeholder SNP ID. "
                        "Ensure IDs are present and non-'.'."
                    )

                if snp_id in seen_ids:
                    sys.exit(
                        f"ERROR: Line {line_number}: Duplicate SNP ID '{snp_id}'. "
                        "Ensure SNP IDs are unique."
                    )
                seen_ids.add(snp_id)

                # ALT-based multiallelic detection
                alt_field = parts[4]
                if "," in alt_field:
                    sys.exit(
                        f"ERROR: Line {line_number}: Multiallelic site detected (ALT='{alt_field}'). "
                        "Filter or split multiallelic sites before processing."
                    )

                # Validate FORMAT and sample fields
                format_field = parts[8] if len(parts) > 8 else "GT"
                sample_fields = parts[9:]

                geno: List[float] = []
                for sample_idx, sample_field in enumerate(sample_fields):
                    gt = extract_gt(format_field, sample_field)

                    # Detect multiallelic allele indices in GT
                    try:
                        alleles = [
                            int(x)
                            for x in gt.replace("|", "/").split("/")
                            if x not in ("", ".")
                        ]
                    except ValueError:
                        alleles = []

                    if any(a > 1 for a in alleles):
                        sys.exit(
                            f"ERROR: Line {line_number}, Sample {sample_idx + 1}: "
                            f"Multiallelic genotype detected (GT='{gt}'). "
                            "Filter or split multiallelic sites before processing."
                        )

                    geno.append(gt_to_value(gt, ploidy, convert_het))

                # Store validated data
                maf = calculate_maf(geno, ploidy)
                snp_data = SNPData(genotypes=geno, sample_fields=sample_fields, maf=maf)
                snp_genotypes[snp_id] = snp_data
                snp_maf_cache[snp_id] = maf

        # Final validation
        if data_line_count == 0:
            sys.exit(
                "ERROR: No data lines found in VCF. "
                "Ensure VCF contains variant records."
            )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"VCF validation complete: {data_line_count} variants, "
                f"{len(headers.samples)} samples processed successfully."
            )

        return snp_genotypes, snp_maf_cache

    def _build_snp_metadata(
        self, vcf_path: Path, headers: VCFHeaders, ploidy: int, convert_het: bool
    ) -> Tuple[Dict[str, SNPMetadata], Dict[str, float]]:
        """Build SNP metadata without loading all genotype data (for lazy loading)."""
        seen_ids: Set[str] = set()
        data_line_count = 0
        snp_metadata: Dict[str, SNPMetadata] = {}
        snp_maf_cache: Dict[str, float] = {}

        with open(
            vcf_path, "rb"
        ) as vcf_file:  # Binary mode for accurate byte positions
            line_number = 0

            for line_bytes in vcf_file:
                line_start_pos = vcf_file.tell() - len(line_bytes)
                line = line_bytes.decode("utf-8").strip()
                line_number += 1

                if line.startswith("#") or not line:
                    continue

                data_line_count += 1
                parts = line.split("\t")

                # Validate line format (same validation as _validate_variants)
                expected_columns = 9 + len(headers.samples)
                if len(parts) != expected_columns:
                    sys.exit(
                        f"ERROR: Line {line_number}: Invalid number of columns. "
                        f"Expected {expected_columns}, found {len(parts)}."
                    )

                # Validate required columns exist
                if len(parts) < 5:
                    sys.exit(
                        f"ERROR: Line {line_number}: Missing required VCF columns "
                        "(CHROM, POS, ID, REF, ALT)."
                    )

                # Validate SNP ID
                snp_id = parts[2]
                if (not snp_id) or (snp_id == "."):
                    sys.exit(
                        f"ERROR: Line {line_number}: Missing or placeholder SNP ID. "
                        "Ensure IDs are present and non-'.'."
                    )

                if snp_id in seen_ids:
                    sys.exit(
                        f"ERROR: Line {line_number}: Duplicate SNP ID '{snp_id}'. "
                        "Ensure SNP IDs are unique."
                    )
                seen_ids.add(snp_id)

                # ALT-based multiallelic detection
                alt_field = parts[4]
                if "," in alt_field:
                    sys.exit(
                        f"ERROR: Line {line_number}: Multiallelic site detected (ALT='{alt_field}'). "
                        "Filter or split multiallelic sites before processing."
                    )

                # Validate FORMAT and sample fields
                format_field = parts[8] if len(parts) > 8 else "GT"
                sample_fields = parts[9:]

                # Quick MAF calculation without storing all genotypes
                geno: List[float] = []
                for sample_idx, sample_field in enumerate(sample_fields):
                    gt = extract_gt(format_field, sample_field)

                    # Detect multiallelic allele indices in GT
                    try:
                        alleles = [
                            int(x)
                            for x in gt.replace("|", "/").split("/")
                            if x not in ("", ".")
                        ]
                    except ValueError:
                        alleles = []

                    if any(a > 1 for a in alleles):
                        sys.exit(
                            f"ERROR: Line {line_number}, Sample {sample_idx + 1}: "
                            f"Multiallelic genotype detected (GT='{gt}'). "
                            "Filter or split multiallelic sites before processing."
                        )

                    geno.append(gt_to_value(gt, ploidy, convert_het))

                # Calculate MAF for metadata
                maf = calculate_maf(geno, ploidy)

                # Store metadata (not full genotype data)
                metadata = SNPMetadata(
                    snp_id=snp_id,
                    maf=maf,
                    file_position=line_start_pos,
                    line_length=len(line_bytes),
                    line_number=line_number,
                )
                snp_metadata[snp_id] = metadata
                snp_maf_cache[snp_id] = maf

        # Final validation
        if data_line_count == 0:
            sys.exit(
                "ERROR: No data lines found in VCF. "
                "Ensure VCF contains variant records."
            )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(
                f"VCF metadata build complete: {data_line_count} variants, "
                f"{len(headers.samples)} samples (lazy loading enabled)"
            )

        return snp_metadata, snp_maf_cache
