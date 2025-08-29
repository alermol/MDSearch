"""Input validation utilities."""

import sys
import argparse

__all__ = ["validate_cli_arguments"]


def validate_cli_arguments(args: argparse.Namespace) -> None:
    """Validate CLI argument combinations and constraints.

    Args:
        args: Parsed command line arguments
    """
    # Validate CLI numeric constraints early
    if args.ns is None or args.ns < 1:
        sys.exit("-ns (number of sets) must be >= 1")
