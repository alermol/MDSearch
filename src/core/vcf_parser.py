"""VCF file parsing and validation."""

import sys
import logging
from pathlib import Path
from typing import List, Dict, Set, Any, Tuple, Optional, Iterator
from dataclasses import dataclass, field
from collections import OrderedDict

from ..utils.memory_monitor import MemoryMonitor
from .genotype_utils import extract_gt, gt_to_value, calculate_maf
import pysam

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
    """Lightweight SNP metadata without genotype data for lazy loading.

    For plain-text uncompressed VCFs, we store byte offsets for fast seeking.
    For compressed VCFs/BCFs, random access requires an index; we store
    chrom/pos to fetch on demand via pysam.
    """

    snp_id: str
    maf: float
    file_position: int  # Byte position where this SNP line starts (-1 if unknown)
    line_length: int  # Length of the line for seeking (0 if unknown)
    line_number: int  # Line number in file for error reporting
    chrom: Optional[str] = None
    pos: Optional[int] = None


class LRUCache:
    """Simple LRU cache for SNP data."""

    def __init__(self, max_size: int = 1000):
        """Initialize LRU cache with maximum size.

        Args:
            max_size: Maximum number of SNPs to keep in cache

        Example:
            >>> cache = LRUCache(max_size=500)
            >>> cache.max_size
            500
        """
        self.max_size = max_size
        self.cache: OrderedDict[str, SNPData] = OrderedDict()

    def get(self, key: str) -> Optional[SNPData]:
        """Get item from cache and move it to end (most recently used).

        Args:
            key: SNP ID to retrieve

        Returns:
            SNPData if found, None otherwise

        Example:
            >>> cache = LRUCache()
            >>> snp_data = SNPData(genotypes=[0.0, 1.0], sample_fields=["0/0", "0/1"], maf=0.5)
            >>> cache.put("rs123", snp_data)
            >>> retrieved = cache.get("rs123")
            >>> retrieved.maf
            0.5
        """
        if key in self.cache:
            # Move to end (most recently used)
            self.cache.move_to_end(key)
            return self.cache[key]
        return None

    def put(self, key: str, value: SNPData) -> None:
        """Add item to cache, evicting least recently used if needed.

        Args:
            key: SNP ID to store
            value: SNPData to cache

        Example:
            >>> cache = LRUCache(max_size=2)
            >>> snp1 = SNPData(genotypes=[0.0, 1.0], sample_fields=["0/0", "0/1"], maf=0.5)
            >>> snp2 = SNPData(genotypes=[1.0, 0.0], sample_fields=["1/1", "0/0"], maf=0.5)
            >>> snp3 = SNPData(genotypes=[0.0, 0.0], sample_fields=["0/0", "0/0"], maf=0.0)
            >>> cache.put("rs1", snp1)
            >>> cache.put("rs2", snp2)
            >>> cache.put("rs3", snp3)  # This will evict rs1
            >>> "rs1" in cache
            False
        """
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
        """Clear all cached data.

        Example:
            >>> cache = LRUCache()
            >>> cache.put("rs123", SNPData(genotypes=[0.0], sample_fields=["0/0"], maf=0.0))
            >>> cache.size()
            1
            >>> cache.clear()
            >>> cache.size()
            0
        """
        self.cache.clear()

    def size(self) -> int:
        """Return current cache size.

        Returns:
            Number of SNPs currently cached

        Example:
            >>> cache = LRUCache(max_size=1000)
            >>> cache.size()
            0
            >>> cache.put("rs123", SNPData(genotypes=[0.0], sample_fields=["0/0"], maf=0.0))
            >>> cache.size()
            1
        """
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
        """Get SNP data with lazy loading and caching.

        Args:
            snp_id: SNP identifier to retrieve

        Returns:
            SNPData containing genotypes and metadata

        Raises:
            KeyError: If SNP ID not found in metadata

        Example:
            >>> # Assuming lazy_vcf_data is initialized with VCF file
            >>> snp_data = lazy_vcf_data.get_snp_data("rs123")
            >>> print(f"SNP {snp_data.maf:.2f} MAF, {len(snp_data.genotypes)} samples")
            SNP 0.45 MAF, 100 samples
        """
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
            # If we have a valid byte offset, use fast text seek
            if metadata.file_position >= 0 and metadata.line_length > 0:
                with open(self._file_path, "r") as vcf_file:
                    vcf_file.seek(metadata.file_position)
                    line = vcf_file.read(metadata.line_length).strip()

                    parts = line.split("\t")
                    if len(parts) < 9:
                        raise ValueError(
                            f"Invalid VCF line format at position {metadata.file_position}"
                        )

                    format_field = parts[8]
                    sample_fields = parts[9:]

                    geno: List[float] = []
                    for sample_field in sample_fields:
                        gt = extract_gt(format_field, sample_field)
                        geno.append(gt_to_value(gt, self._ploidy, self._convert_het))

                    return SNPData(
                        genotypes=geno, sample_fields=sample_fields, maf=metadata.maf
                    )

            # Otherwise, try to fetch via pysam using chrom/pos (requires index)
            if metadata.chrom is None or metadata.pos is None:
                raise RuntimeError(
                    "Missing chrom/pos for random access in compressed input"
                )

            with pysam.VariantFile(str(self._file_path)) as ivcf:
                found = None
                for rec in ivcf.fetch(metadata.chrom, metadata.pos - 1, metadata.pos):
                    if rec.id == metadata.snp_id:
                        found = rec
                        break
                if found is None:
                    raise RuntimeError(
                        f"Record {metadata.snp_id} not found at {metadata.chrom}:{metadata.pos}"
                    )

                # Build genotype values from GT tuples
                sample_fields_list: List[str] = []
                geno_values: List[float] = []
                for sample in ivcf.header.samples:
                    data = found.samples[sample]
                    gt_tuple = data.get("GT")
                    if gt_tuple is None or all(g is None for g in gt_tuple):
                        gt_str = "."
                    else:
                        # Assume unphased for conversion
                        gt_parts = ["." if g is None else str(g) for g in gt_tuple]
                        sep = "/" if len(gt_parts) != 1 else "/"
                        gt_str = sep.join(gt_parts)
                    # Reconstruct sample field with GT-only
                    sample_fields_list.append(gt_str)
                    geno_values.append(
                        gt_to_value(gt_str, self._ploidy, self._convert_het)
                    )

                return SNPData(
                    genotypes=geno_values,
                    sample_fields=sample_fields_list,
                    maf=metadata.maf,
                )

        except Exception as e:
            raise RuntimeError(
                f"Failed to load SNP '{metadata.snp_id}' from line {metadata.line_number}: {e}"
            ) from e

    @property
    def snp_genotypes(self) -> "LazyGenotypesDict":
        """Provide backward compatibility by creating a dict-like interface.

        Returns:
            LazyGenotypesDict that provides dict-like access to SNP data

        Example:
            >>> # Access SNP data like a regular dictionary
            >>> snp_data = lazy_vcf_data.snp_genotypes["rs123"]
            >>> print(f"Found {len(snp_data.genotypes)} genotypes")
            Found 100 genotypes
        """
        return LazyGenotypesDict(self)

    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache performance statistics.

        Returns:
            Dictionary with cache size, max size, and total SNPs

        Example:
            >>> stats = lazy_vcf_data.get_cache_stats()
            >>> print(f"Cache: {stats['cache_size']}/{stats['max_cache_size']} SNPs")
            Cache: 45/1000 SNPs
        """
        return {
            "cache_size": self._cache.size(),
            "max_cache_size": self._cache.max_size,
            "total_snps": len(self.snp_metadata),
        }


class LazyGenotypesDict:
    """Dict-like interface for lazy loaded SNP genotypes."""

    def __init__(self, lazy_vcf_data: LazyVCFData):
        """Initialize with reference to lazy VCF data.

        Args:
            lazy_vcf_data: LazyVCFData instance to wrap
        """
        self._lazy_data = lazy_vcf_data

    def __getitem__(self, snp_id: str) -> SNPData:
        """Get SNP data by ID, triggering lazy loading if needed.

        Args:
            snp_id: SNP identifier

        Returns:
            SNPData for the requested SNP

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> snp_data = genotypes["rs123"]  # Triggers lazy loading
            >>> print(f"SNP has {len(snp_data.genotypes)} samples")
            SNP has 100 samples
        """
        return self._lazy_data.get_snp_data(snp_id)

    def __contains__(self, snp_id: str) -> bool:
        """Check if SNP ID exists in metadata.

        Args:
            snp_id: SNP identifier to check

        Returns:
            True if SNP exists, False otherwise

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> "rs123" in genotypes
            True
            >>> "nonexistent" in genotypes
            False
        """
        return snp_id in self._lazy_data.snp_metadata

    def __iter__(self) -> Iterator[str]:
        """Iterate over SNP IDs.

        Returns:
            Iterator of SNP identifiers

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> snp_ids = list(genotypes)[:5]  # Get first 5 SNP IDs
            >>> print(f"First 5 SNPs: {snp_ids}")
            First 5 SNPs: ['rs1', 'rs2', 'rs3', 'rs4', 'rs5']
        """
        return iter(self._lazy_data.snp_metadata.keys())

    def __len__(self) -> int:
        """Return total number of SNPs.

        Returns:
            Number of SNPs in the dataset

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> print(f"Total SNPs: {len(genotypes)}")
            Total SNPs: 1000
        """
        return len(self._lazy_data.snp_metadata)

    def items(self) -> Iterator[Tuple[str, SNPData]]:
        """Iterate over SNP ID and data pairs (loads data lazily).

        Returns:
            Iterator of (SNP ID, SNPData) tuples

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> for snp_id, snp_data in genotypes.items():
            ...     print(f"{snp_id}: MAF={snp_data.maf:.2f}")
            ...     break  # Just show first one
            rs1: MAF=0.45
        """
        for snp_id in self._lazy_data.snp_metadata:
            yield snp_id, self._lazy_data.get_snp_data(snp_id)

    def keys(self) -> Iterator[str]:
        """Return SNP IDs.

        Returns:
            Iterator of SNP identifiers

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> snp_ids = list(genotypes.keys())[:3]
            >>> print(f"SNP IDs: {snp_ids}")
            SNP IDs: ['rs1', 'rs2', 'rs3']
        """
        return iter(self._lazy_data.snp_metadata.keys())

    def values(self) -> Iterator[SNPData]:
        """Return SNP data (loads all data - use carefully!).

        Returns:
            Iterator of SNPData objects

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> # Be careful - this loads ALL SNP data into memory
            >>> all_snps = list(genotypes.values())
            >>> print(f"Loaded {len(all_snps)} SNPs")
            Loaded 1000 SNPs
        """
        for snp_id in self._lazy_data.snp_metadata:
            yield self._lazy_data.get_snp_data(snp_id)

    def get(self, snp_id: str, default: Optional[SNPData] = None) -> Optional[SNPData]:
        """Get SNP data with optional default.

        Args:
            snp_id: SNP identifier to retrieve
            default: Value to return if SNP not found

        Returns:
            SNPData if found, default value otherwise

        Example:
            >>> genotypes = lazy_vcf_data.snp_genotypes
            >>> snp_data = genotypes.get("rs123", None)
            >>> if snp_data:
            ...     print(f"Found SNP with MAF {snp_data.maf}")
            ... else:
            ...     print("SNP not found")
            Found SNP with MAF 0.45
        """
        try:
            return self._lazy_data.get_snp_data(snp_id)
        except KeyError:
            return default


class VCFParser:
    """Handles VCF file parsing and validation."""

    def __init__(self, memory_monitor: MemoryMonitor, logger: logging.Logger):
        """Initialize VCF parser with memory monitoring and logging.

        Args:
            memory_monitor: MemoryMonitor instance for tracking memory usage
            logger: Logger instance for output

        Example:
            >>> from src.utils.memory_monitor import MemoryMonitor
            >>> from src.utils.logging_setup import setup_logger
            >>> logger = setup_logger("vcf_parser")
            >>> memory_monitor = MemoryMonitor(logger)
            >>> parser = VCFParser(memory_monitor, logger)
        """
        self.memory_monitor = memory_monitor
        self.logger = logger

    def parse_and_validate(
        self, vcf_path: Path, ploidy: int, convert_het: bool
    ) -> VCFData:
        """Parse VCF file and return structured data (loads all SNPs into memory).

        Args:
            vcf_path: Path to VCF file to parse
            ploidy: Ploidy level (e.g., 2 for diploid)
            convert_het: Whether to convert heterozygous calls to missing

        Returns:
            VCFData containing all SNP information

        Example:
            >>> parser = VCFParser(memory_monitor, logger)
            >>> vcf_data = parser.parse_and_validate(
            ...     Path("sample.vcf"), ploidy=2, convert_het=False
            ... )
            >>> print(f"Parsed {len(vcf_data.snp_genotypes)} SNPs")
            >>> print(f"Found {len(vcf_data.headers.samples)} samples")
            Parsed 1000 SNPs
            Found 50 samples
        """
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
        """Parse VCF file and return lazy-loaded data structure (minimal memory usage).

        Args:
            vcf_path: Path to VCF file to parse
            ploidy: Ploidy level (e.g., 2 for diploid)
            convert_het: Whether to convert heterozygous calls to missing
            cache_size: Maximum number of SNPs to keep in memory cache

        Returns:
            LazyVCFData with lazy loading capabilities

        Example:
            >>> parser = VCFParser(memory_monitor, logger)
            >>> lazy_vcf_data = parser.parse_and_validate_lazy(
            ...     Path("large_sample.vcf"), ploidy=2, convert_het=False, cache_size=500
            ... )
            >>> print(f"Metadata for {len(lazy_vcf_data.snp_metadata)} SNPs loaded")
            >>> # Access individual SNPs as needed
            >>> snp_data = lazy_vcf_data.get_snp_data("rs123")
            Metadata for 10000 SNPs loaded
        """
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
        try:
            with pysam.VariantFile(str(vcf_path)) as vf:
                samples = list(vf.header.samples)
                # pysam provides header version via .header.version in newer releases
                fileformat_found = True if getattr(vf.header, "version", None) else True

                if len(samples) < 2:
                    sys.exit(
                        f"ERROR: Insufficient samples for distance calculation. "
                        f"Found {len(samples)} samples, need at least 2."
                    )

                return VCFHeaders(samples=samples, fileformat_found=fileformat_found)
        except Exception as e:
            sys.exit(f"ERROR: Failed to read VCF/BCF headers via pysam: {e}")

    def _validate_variants(
        self, vcf_path: Path, headers: VCFHeaders, ploidy: int, convert_het: bool
    ) -> Tuple[Dict[str, SNPData], Dict[str, float]]:
        """Parse and validate variant lines."""
        seen_ids: Set[str] = set()
        data_line_count = 0
        snp_genotypes: Dict[str, SNPData] = {}
        snp_maf_cache: Dict[str, float] = {}

        with pysam.VariantFile(str(vcf_path)) as vf:
            line_number = 0
            for rec in vf:
                line_number += 1
                # Skip header pseudo-lines; pysam yields only records
                data_line_count += 1

                snp_id = rec.id or ""
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
                if rec.alts and len(rec.alts) > 1:
                    sys.exit(
                        f"ERROR: Line {line_number}: Multiallelic site detected (ALT='{','.join(rec.alts)}'). "
                        "Filter or split multiallelic sites before processing."
                    )

                # Validate sample count
                if len(rec.samples) != len(headers.samples):
                    sys.exit(
                        f"ERROR: Line {line_number}: Invalid number of samples. "
                        f"Expected {len(headers.samples)}, found {len(rec.samples)}."
                    )

                geno: List[float] = []
                sample_fields_list: List[str] = []
                for sample_idx, sample in enumerate(headers.samples):
                    data = rec.samples[sample]
                    gt_tuple = data.get("GT")
                    if gt_tuple is None or all(g is None for g in gt_tuple):
                        gt = "."
                    else:
                        gt_parts = ["." if g is None else str(g) for g in gt_tuple]
                        # Assume unphased representation for numeric conversion
                        gt = "/".join(gt_parts)

                    # Detect multiallelic indices in GT
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

                    sample_fields_list.append(gt)
                    geno.append(gt_to_value(gt, ploidy, convert_het))

                maf = calculate_maf(geno, ploidy)
                snp_data = SNPData(
                    genotypes=geno, sample_fields=sample_fields_list, maf=maf
                )
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

        # Decide strategy based on file extension: uncompressed .vcf supports byte offsets
        is_plain_vcf = str(vcf_path).endswith(".vcf")

        if is_plain_vcf:
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

                    alt_field = parts[4]
                    if "," in alt_field:
                        sys.exit(
                            f"ERROR: Line {line_number}: Multiallelic site detected (ALT='{alt_field}'). "
                            "Filter or split multiallelic sites before processing."
                        )

                    format_field = parts[8] if len(parts) > 8 else "GT"
                    sample_fields = parts[9:]

                    geno: List[float] = []
                    for sample_idx, sample_field in enumerate(sample_fields):
                        gt = extract_gt(format_field, sample_field)
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

                    maf = calculate_maf(geno, ploidy)
                    metadata = SNPMetadata(
                        snp_id=snp_id,
                        maf=maf,
                        file_position=line_start_pos,
                        line_length=len(line_bytes),
                        line_number=line_number,
                        chrom=parts[0],
                        pos=int(parts[1]) if parts[1].isdigit() else None,
                    )
                    snp_metadata[snp_id] = metadata
                    snp_maf_cache[snp_id] = maf
        else:
            # Compressed VCF/BCF: iterate with pysam and store chrom/pos for fetch
            try:
                with pysam.VariantFile(str(vcf_path)) as vf:
                    line_number = 0
                    for rec in vf:
                        line_number += 1
                        data_line_count += 1

                        snp_id = rec.id or ""
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

                        if rec.alts and len(rec.alts) > 1:
                            sys.exit(
                                f"ERROR: Line {line_number}: Multiallelic site detected (ALT='{','.join(rec.alts)}'). "
                                "Filter or split multiallelic sites before processing."
                            )

                        geno_values: List[float] = []
                        for sample_idx, sample in enumerate(headers.samples):
                            data = rec.samples[sample]
                            gt_tuple = data.get("GT")
                            if gt_tuple is None or all(g is None for g in gt_tuple):
                                gt = "."
                            else:
                                gt_parts = [
                                    "." if g is None else str(g) for g in gt_tuple
                                ]
                                gt = "/".join(gt_parts)

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
                            geno_values.append(gt_to_value(gt, ploidy, convert_het))

                        maf = calculate_maf(geno_values, ploidy)
                        metadata = SNPMetadata(
                            snp_id=snp_id,
                            maf=maf,
                            file_position=-1,
                            line_length=0,
                            line_number=line_number,
                            chrom=rec.chrom,
                            pos=int(rec.pos),
                        )
                        snp_metadata[snp_id] = metadata
                        snp_maf_cache[snp_id] = maf
            except Exception as e:
                sys.exit(f"ERROR: Failed to read compressed VCF/BCF via pysam: {e}")

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
