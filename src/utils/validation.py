"""Input validation utilities."""

import sys
import argparse

__all__ = ["validate_cli_arguments"]


def validate_cli_arguments(args: argparse.Namespace) -> None:
    """Validate CLI argument combinations and constraints.

    Args:
        args: Parsed command line arguments
    """
    if args.ns == "unlimited" or args.ns == "0":
        args.ns = 0
    elif isinstance(args.ns, str) and args.ns.isdigit() and int(args.ns) >= 1:
        args.ns = int(args.ns)
    elif isinstance(args.ns, int) and args.ns >= 1:
        pass
    else:
        sys.exit("-ns (number of sets) must be >= 1, 0, or 'unlimited'")

    if not (0.0 <= args.weight_entropy <= 1.0):
        sys.exit("-we (weight_entropy) must be between 0.0 and 1.0")

    if not (0.0 <= args.weight_maf <= 1.0):
        sys.exit("-wm (weight_maf) must be between 0.0 and 1.0")

    if args.weight_entropy + args.weight_maf > 1.0:
        sys.exit("-we (weight_entropy) and -wm (weight_maf) must sum to <= 1.0")
