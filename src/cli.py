"""Command-line interface for MDSearch."""

import argparse
import signal
import sys
from pathlib import Path
import subprocess
import platform

from .app import MDSearchApp, MDSearchConfig
from .utils.validation import validate_cli_arguments
from . import __version__

__all__ = ["parser_resolve_path", "create_parser", "main", "is_shutdown_requested"]

_shutdown_requested = False


def _signal_handler(signum: int, frame: object) -> None:
    """Handle signals for graceful shutdown."""
    global _shutdown_requested
    if signum == signal.SIGINT:
        print("\nReceived interrupt signal (Ctrl+C). Shutting down gracefully...")
    elif signum == signal.SIGTERM:
        print("\nReceived termination signal. Shutting down gracefully...")
    else:
        print(f"\nReceived signal {signum}. Shutting down gracefully...")
    _shutdown_requested = True


def is_shutdown_requested() -> bool:
    """Check if graceful shutdown was requested."""
    return _shutdown_requested


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

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=_build_version_string(),
        help="Show program version, commit hash, and Python version, then exit",
    )

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
        help="Number of distinct SNP sets in output (0 or 'unlimited' for all possible sets)",
        default="0",
        type=str,
        metavar="N_SETS",
    )
    grp_select.add_argument(
        "-we",
        "--weight-entropy",
        dest="weight_entropy",
        help="Weight for entropy component in SNP scoring (0.0 to 1.0)",
        default=0.5,
        type=float,
        metavar="WEIGHT_ENTROPY",
    )
    grp_select.add_argument(
        "-wm",
        "--weight-maf",
        dest="weight_maf",
        help="Weight for MAF component in SNP scoring (0.0 to 1.0)",
        default=0.5,
        type=float,
        metavar="WEIGHT_MAF",
    )

    grp_output = parser.add_argument_group("Output", "Output formatting and summary")
    grp_output.add_argument(
        "-ch",
        "--convert-het",
        dest="ch",
        help="Convert heterozygous calls into NA",
        action="store_true",
        default=False,
    )

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
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)

    parser = create_parser()
    args = parser.parse_args()

    validate_cli_arguments(args)

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
        weight_entropy=args.weight_entropy,
        weight_maf=args.weight_maf,
        input_format=args.input_format,
        output_format=args.output_format,
    )

    print(f"Starting MDSearch {_build_version_string()}")
    print(f"Input VCF: {config.input_vcf}")
    print(f"Output directory: {config.output_prefix}")
    print(
        f"Configuration: ploidy={config.ploidy}, min_distance={config.min_distance}, n_sets={config.n_sets}"
    )
    if config.convert_het:
        print("Converting heterozygous calls to missing values")
    print(f"Weights: entropy={config.weight_entropy}, MAF={config.weight_maf}")
    print(f"Logging: level={config.log_level or 'INFO'}, format={config.log_format}")
    print("-" * 60)

    app = MDSearchApp(config, shutdown_checker=is_shutdown_requested)
    try:
        app.run()
    except KeyboardInterrupt:
        print("\nOperation interrupted by user. Exiting gracefully.")
        sys.exit(1)
    except Exception as e:
        if _shutdown_requested:
            print(f"\nGraceful shutdown completed. Error during shutdown: {e}")
            sys.exit(1)
        else:
            raise


if __name__ == "__main__":
    main()
