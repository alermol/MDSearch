"""Index verification and creation for variant files (VCF/BCF)."""

from pathlib import Path
from typing import Optional
import logging
import subprocess

import pysam

__all__ = ["ensure_variant_index", "infer_format_letter"]


def infer_format_letter(path: Path) -> str:
    """Infer bcftools-style format letter from filename.

    Returns one of: v (VCF), z (VCF.gz), u/b (BCF). Defaults to 'v'.
    """
    name = str(path)
    if name.endswith(".vcf.gz") or name.endswith(".vcfz"):
        return "z"
    if name.endswith(".bcf"):
        # Compression type isn't encoded in extension; use 'b' as default
        return "b"
    if name.endswith(".vcf"):
        return "v"
    # Fallback: try to open and inspect
    try:
        with pysam.VariantFile(name) as vf:
            text_header = str(vf.header)
            # Heuristic: bgzipped VCFs usually need tabix; plain VCFs are not compressed
            if name.endswith(".vcf.gz"):
                return "z"
            if "BCF" in text_header:
                return "b"
    except Exception:
        pass
    return "v"


def ensure_variant_index(
    vpath: Path, fmt_letter: Optional[str], logger: Optional[logging.Logger]
) -> None:
    """Ensure an index exists for the given variant file.

    - For 'z' (VCF.gz): ensure .tbi or .csi exists; create with tabix if missing
    - For 'b'/'u' (BCF): ensure .csi exists; create with pysam.index if missing
    - For 'v' (VCF): no index is created; log and return
    - For 'auto': infer from extension
    """
    fmt = fmt_letter or "auto"
    if fmt == "auto":
        fmt = infer_format_letter(vpath)

    if fmt == "v":
        if logger and logger.isEnabledFor(logging.INFO):
            logger.info("Input is plain VCF; no index required")
        return

    if fmt in ("b", "u", "z"):
        # CSI naming convention appends .csi to the full filename
        csi = Path(str(vpath) + ".csi")
        if csi.exists():
            return
        if logger and logger.isEnabledFor(logging.INFO):
            logger.info("Index not found; creating CSI index via bcftools...")
        # Create CSI index using bcftools
        _run_bcftools_index(vpath, logger)
        return

    # Unknown: do nothing
    if logger and logger.isEnabledFor(logging.WARNING):
        logger.warning(f"Unknown input format '{fmt}'; skipping index verification")


def _run_bcftools_index(vpath: Path, logger: Optional[logging.Logger]) -> None:
    """Run bcftools index --csi on the given file, raising on failure."""
    cmd = ["bcftools", "index", "--csi", "-f", str(vpath)]
    try:
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if logger and logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"bcftools index stdout: {res.stdout.strip()}")
    except FileNotFoundError:
        msg = "bcftools not found in PATH; please install bcftools to enable CSI indexing"
        if logger:
            logger.error(msg)
        raise RuntimeError(msg)
    except subprocess.CalledProcessError as e:
        if logger:
            logger.error(f"bcftools index failed: {e.stderr.strip()}")
        raise RuntimeError(f"bcftools index failed with code {e.returncode}")


