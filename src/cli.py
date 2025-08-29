"""Command-line interface for MDSearch."""

import argparse
from pathlib import Path
import subprocess
import platform

from .app import MDSearchApp, MDSearchConfig
from .utils.validation import validate_cli_arguments
from . import __version__

__all__ = ["parser_resolve_path", "create_parser", "main"]


def _get_git_commit() -> str:
    """Return short git commit hash if available, else 'unknown'."""
    try:
        res = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return res.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _build_version_string() -> str:
    """Compose version string with build and runtime info."""
    commit = _get_git_commit()
    py = platform.python_version()
    return f"MDSearch {__version__} (commit hash {commit})\nPython {py}"


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

    # Version flag (Misc)
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=_build_version_string(),
        help="Show program version, commit hash, and Python version, then exit",
    )

    # Positional (required) arguments
    parser.add_argument(
        "ivcf",
        help="Input VCF file (bi-allelic, with SNP IDs)",
        type=parser_resolve_path,
        metavar="IVCF",
    )
    parser.add_argument(
        "outdir",
        help="Path to output folder (will be created if absent); VCF files will be placed in 'mdss' subdirectory",
        type=parser_resolve_path,
        metavar="OUTPUT_FOLDER",
    )

    # Group: Selection parameters
    grp_select = parser.add_argument_group(
        "Selection", "Core parameters controlling selection"
    )
    grp_select.add_argument(
        "-pl",
        "--ploidy",
        dest="pl",
        help="VCF ploidy",
        default=2,
        type=int,
        metavar="PLOIDY",
    )
    grp_select.add_argument(
        "-md",
        "--min-distance",
        dest="md",
        help="Minimal Hamming distance between samples",
        default=1,
        type=int,
        metavar="MIN_DIST",
    )
    grp_select.add_argument(
        "-ns",
        "--num-sets",
        dest="ns",
        help="Number of distinct SNP sets in output",
        default=1,
        type=int,
        metavar="N_SETS",
    )

    # Group: Output & formatting
    grp_output = parser.add_argument_group("Output", "Output formatting and summary")
    grp_output.add_argument(
        "-ch",
        "--convert-het",
        dest="ch",
        help="Convert heterozygous calls into NA",
        action="store_true",
        default=False,
    )

    # Group: IO formats
    grp_io = parser.add_argument_group("IO formats", "Input and output VCF/BCF formats")
    grp_io.add_argument(
        "-I",
        "--input-format",
        help=(
            "Input format (bcftools-style): auto (default), v (VCF), z (VCF.gz), "
            "u (uncompressed BCF), b (compressed BCF)"
        ),
        choices=["auto", "v", "z", "u", "b"],
        default="auto",
    )
    grp_io.add_argument(
        "-O",
        "--output-format",
        help=(
            "Output format (bcftools-style): v (VCF), z (VCF.gz), "
            "u (uncompressed BCF), b (compressed BCF)"
        ),
        choices=["v", "z", "u", "b"],
        default="v",
    )

    # Group: Performance & memory
    # Lazy loading options removed - no longer needed

    # Group: Logging
    grp_log = parser.add_argument_group("Logging", "Logging verbosity and format")
    grp_log.add_argument(
        "-q",
        "--quiet",
        help="Suppress progress output",
        action="store_true",
        default=False,
    )
    grp_log.add_argument(
        "-L",
        "--log-level",
        help=(
            "Logging level (DEBUG, INFO, WARNING, ERROR); default depends on --quiet"
        ),
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default=None,
    )
    grp_log.add_argument(
        "-F",
        "--log-format",
        help="Logging format: text or json",
        choices=["text", "json"],
        default="text",
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
        output_prefix=args.outdir,
        ploidy=args.pl,
        min_distance=args.md,
        convert_het=args.ch,
        n_sets=args.ns,
        verbose=not args.quiet,
        log_level=args.log_level,
        log_format=args.log_format,
        # Lazy loading arguments removed - no longer needed
        input_format=args.input_format,
        output_format=args.output_format,
    )

    # Run application
    app = MDSearchApp(config)
    app.run()


if __name__ == "__main__":
    main()
