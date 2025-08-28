"""Input validation utilities."""

import sys
import argparse


def validate_cli_arguments(args: argparse.Namespace) -> None:
    """Validate CLI argument combinations and constraints."""
    # Validate mutually exclusive overlap constraints
    validate_overlap_constraints(args.oMx, args.oMf)
    
    # Validate CLI numeric constraints early
    if args.ts is not None and args.ts < 0:
        sys.exit("-ts (total SNPs) must be >= 0")
    if args.ns is None or args.ns < 1:
        sys.exit("-ns (number of sets) must be >= 1")


def validate_overlap_constraints(overlap_number: int, overlap_fraction: float) -> None:
    """Validate mutual exclusivity of overlap constraints."""
    if (overlap_number is not None and overlap_number >= 0) and (
        overlap_fraction is not None and overlap_fraction >= 0
    ):
        sys.exit(
            "Specify only one of -oMx (max overlap number) or -oMf "
            "(max overlap fraction), not both."
        )
