"""Input validation utilities."""

import sys
import argparse

__all__ = ["validate_cli_arguments", "validate_overlap_constraints"]


def validate_cli_arguments(args: argparse.Namespace) -> None:
    """Validate CLI argument combinations and constraints.

    Args:
        args: Parsed command line arguments

    Raises:
        SystemExit: If validation fails

    Example:
        >>> import argparse
        >>> parser = argparse.ArgumentParser()
        >>> parser.add_argument("-ts", type=int, default=0)
        >>> parser.add_argument("-ns", type=int, default=1)
        >>> parser.add_argument("-oMx", type=int, default=-1)
        >>> parser.add_argument("-oMf", type=float, default=-1.0)
        >>>
        >>> # Valid arguments
        >>> args = parser.parse_args(["-ts", "10", "-ns", "2"])
        >>> validate_cli_arguments(args)  # No error
        >>>
        >>> # Invalid arguments
        >>> args = parser.parse_args(["-ts", "-5"])  # Negative total SNPs
        >>> validate_cli_arguments(args)  # Will exit with error
    """
    # Validate mutually exclusive overlap constraints
    validate_overlap_constraints(args.oMx, args.oMf)

    # Validate CLI numeric constraints early
    if args.ts is not None and args.ts < 0:
        sys.exit("-ts (total SNPs) must be >= 0")
    if args.ns is None or args.ns < 1:
        sys.exit("-ns (number of sets) must be >= 1")


def validate_overlap_constraints(overlap_number: int, overlap_fraction: float) -> None:
    """Validate mutual exclusivity of overlap constraints.

    Args:
        overlap_number: Maximum overlap count (-1 for unlimited)
        overlap_fraction: Maximum overlap fraction (-1.0 for unlimited)

    Raises:
        SystemExit: If both constraints are specified

    Example:
        >>> # Valid: only number constraint
        >>> validate_overlap_constraints(5, -1.0)  # No error
        >>>
        >>> # Valid: only fraction constraint
        >>> validate_overlap_constraints(-1, 0.5)  # No error
        >>>
        >>> # Invalid: both constraints specified
        >>> validate_overlap_constraints(5, 0.5)  # Will exit with error
    """
    if (overlap_number is not None and overlap_number >= 0) and (
        overlap_fraction is not None and overlap_fraction >= 0
    ):
        sys.exit(
            "Specify only one of -oMx (max overlap number) or -oMf "
            "(max overlap fraction), not both."
        )
