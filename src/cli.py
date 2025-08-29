"""Command-line interface for MDSearch."""

import argparse
from pathlib import Path

from .app import MDSearchApp, MDSearchConfig
from .core.snp_selector import OverlapConstraints
from .utils.validation import validate_cli_arguments

__all__ = ["parser_resolve_path", "create_parser", "main"]


def parser_resolve_path(path: str) -> Path:
    """Resolve CLI-provided path string to an absolute Path.

    Args:
        path: Path string from command line

    Returns:
        Resolved absolute Path object

    Example:
        >>> parser_resolve_path("sample.vcf")
        PosixPath('/absolute/path/to/sample.vcf')
        >>> parser_resolve_path("./relative/path.vcf")
        PosixPath('/absolute/path/to/relative/path.vcf')
    """
    return Path(path).resolve()


def create_parser() -> argparse.ArgumentParser:
    """Create and configure argument parser.

    Returns:
        Configured ArgumentParser with all MDSearch options

    Example:
        >>> parser = create_parser()
        >>> args = parser.parse_args(["input.vcf", "output", "-pl", "2", "-md", "3"])
        >>> print(f"Input: {args.ivcf}, Ploidy: {args.pl}, Min distance: {args.md}")
        Input: input.vcf, Ploidy: 2, Min distance: 3
    """
    parser = argparse.ArgumentParser(
        description=(
            "Select minimal discriminatory SNP sets from a VCF given a minimal "
            "pairwise Hamming distance; supports multiple alternative sets and "
            "optional constraints."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Notes: Input VCF must be bi-allelic with SNP IDs present. "
            "Multiallelics are rejected via ALT or GT."
        ),
    )

    parser.add_argument(
        "ivcf",
        help="Input VCF file (bi-allelic, with SNP IDs)",
        type=parser_resolve_path,
        metavar="IVCF",
    )
    parser.add_argument(
        "ovcf_prefix",
        help="Prefix for output VCF(s)",
        type=parser_resolve_path,
        metavar="OVCF_PREFIX",
    )

    parser.add_argument(
        "-pl", help="VCF ploidy (default: 2)", default=2, type=int, metavar="PLOIDY"
    )
    parser.add_argument(
        "-ts",
        help="Total SNPs in output set (0 = keep minimal; default: 0)",
        default=0,
        type=int,
        metavar="TOTAL_SNP",
    )
    parser.add_argument(
        "-md",
        help="Minimal Hamming distance between samples (default: 1)",
        default=1,
        type=int,
        metavar="MIN_DIST",
    )
    parser.add_argument(
        "-ch",
        help="Convert heterozygous calls into NA (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-ns",
        help="Number of distinct SNP sets in output (default: 1)",
        default=1,
        type=int,
        metavar="N_SETS",
    )
    parser.add_argument(
        "-oMx",
        help=(
            "Maximum overlap count allowed with the base minimal set for "
            "alternative sets (-1 = unlimited; default: -1)"
        ),
        default=-1,
        type=int,
        metavar="OVERLAP_MAX_N",
    )
    parser.add_argument(
        "-oMf",
        help=(
            "Maximum overlap fraction allowed with the base minimal set for "
            "alternative sets (-1 = unlimited; default: -1.0)"
        ),
        default=-1.0,
        type=float,
        metavar="OVERLAP_MAX_FRAC",
    )
    parser.add_argument(
        "--quiet",
        help="Suppress progress output (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--log-level",
        help=(
            "Logging level (DEBUG, INFO, WARNING, ERROR); default depends on --quiet"
        ),
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default=None,
    )
    parser.add_argument(
        "--log-format",
        help="Logging format: text or json (default: text)",
        choices=["text", "json"],
        default="text",
    )
    parser.add_argument(
        "--summary-tsv",
        help=(
            "Write per-set summary TSV to this path (columns: set_index, output_vcf, "
            "num_snps, min_distance, snp_ids)."
        ),
        type=parser_resolve_path,
        default=None,
        metavar="SUMMARY_TSV",
    )

    # Memory optimization options
    parser.add_argument(
        "--lazy-loading",
        help="Use lazy loading for large VCF files to reduce memory usage",
        action="store_true",
    )
    parser.add_argument(
        "--cache-size",
        help="Number of SNPs to keep in memory cache when using lazy loading (default: 1000)",
        type=int,
        default=1000,
        metavar="CACHE_SIZE",
    )

    # IO formats (bcftools-style letters)
    parser.add_argument(
        "--input-format",
        help=(
            "Input format (bcftools-style): auto (default), v (VCF), z (VCF.gz), "
            "u (uncompressed BCF), b (compressed BCF)"
        ),
        choices=["auto", "v", "z", "u", "b"],
        default="auto",
    )
    parser.add_argument(
        "--output-format",
        help=(
            "Output format (bcftools-style): v (VCF), z (VCF.gz), "
            "u (uncompressed BCF), b (compressed BCF)"
        ),
        choices=["v", "z", "u", "b"],
        default="v",
    )

    return parser


def main() -> None:
    """CLI entry point.

    This function:
    1. Parses command line arguments
    2. Validates argument combinations
    3. Creates application configuration
    4. Runs the MDSearch pipeline

    Example:
        >>> # Command line usage:
        >>> # python -m src.cli input.vcf output -pl 2 -md 3 -ns 2
        >>> #
        >>> # Programmatic usage:
        >>> import sys
        >>> sys.argv = ["mdsearch", "input.vcf", "output", "-pl", "2", "-md", "3"]
        >>> main()
        >>> # Will process input.vcf and create output_1.vcf, output_2.vcf
    """
    parser = create_parser()
    args = parser.parse_args()

    # Validate arguments
    validate_cli_arguments(args)

    # Create configuration
    config = MDSearchConfig(
        input_vcf=args.ivcf,
        output_prefix=args.ovcf_prefix,
        ploidy=args.pl,
        max_snps=args.ts,
        min_distance=args.md,
        convert_het=args.ch,
        n_sets=args.ns,
        overlap_constraints=OverlapConstraints(
            max_number=args.oMx, max_fraction=args.oMf
        ),
        verbose=not args.quiet,
        log_level=args.log_level,
        log_format=args.log_format,
        summary_tsv=args.summary_tsv,
        lazy_loading=args.lazy_loading,
        cache_size=args.cache_size,
        input_format=args.input_format,
        output_format=args.output_format,
    )

    # Run application
    app = MDSearchApp(config)
    app.run()


if __name__ == "__main__":
    main()
