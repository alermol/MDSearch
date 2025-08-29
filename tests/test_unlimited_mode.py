"""Tests for unlimited mode functionality in MDSearch."""

import sys
import subprocess
from pathlib import Path

from .helpers import (
    write_vcf,
    assert_discriminative,
    save_out_prefix_vcfs,
    get_mdss_vcf_path,
)


def run_mdsearch(ivcf: Path, out_prefix: Path, **kwargs):
    """Run MDSearch with given parameters."""
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
    ]
    # Map kwargs to CLI
    cli_map = {
        "ploidy": "-pl",
        "min_dist": "-md",
        "convert_het": "-ch",
        "n_sets": "-ns",
    }
    for k, v in kwargs.items():
        flag = cli_map[k]
        if isinstance(v, bool):
            if v:
                args.append(flag)
        else:
            args.extend([flag, str(v)])
    subprocess.run(args, check=True)


def test_unlimited_mode_finds_all_sets(tmp_path: Path):
    """Test that unlimited mode finds all possible disjoint SNP sets."""
    samples = ["S1", "S2", "S3", "S4"]
    # Construct multiple alternative minimal sets that are disjoint
    variants = [
        # Set1: A,B - discriminates samples with pattern [0,0], [0,1], [1,0], [1,1]
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "1/1"],
        },
        # Set2: C,D - discriminates samples with pattern [0,0], [0,1], [1,0], [1,1]
        {
            "chrom": "1",
            "pos": 300,
            "id": "C",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "1/1", "0/0"],
        },
        {
            "chrom": "1",
            "pos": 400,
            "id": "D",
            "ref": "A",
            "alt": "T",
            "genotypes": ["1/1", "0/0", "1/1", "0/0"],
        },
        # Set3: E,F - discriminates samples with pattern [0,0], [0,1], [1,0], [1,1]
        {
            "chrom": "1",
            "pos": 500,
            "id": "E",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "0/0"],
        },
        {
            "chrom": "1",
            "pos": 600,
            "id": "F",
            "ref": "A",
            "alt": "T",
            "genotypes": ["1/1", "0/0", "1/1", "1/1"],
        },
        # Additional SNPs to ensure enough are available
        {
            "chrom": "1",
            "pos": 700,
            "id": "G",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "0/0"],
        },
        {
            "chrom": "1",
            "pos": 800,
            "id": "H",
            "ref": "A",
            "alt": "T",
            "genotypes": ["1/1", "1/1", "0/0", "1/1"],
        },
        # More SNPs to ensure we have enough
        {
            "chrom": "1",
            "pos": 900,
            "id": "I",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "1/1", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 1000,
            "id": "J",
            "ref": "A",
            "alt": "T",
            "genotypes": ["1/1", "0/0", "0/0", "0/0"],
        },
    ]
    original_vcf = tmp_path / "orig_unlimited.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_unlimited"

    # Test unlimited mode with -ns 0
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        min_dist=1,
        n_sets=0,
    )

    # Should find 3 disjoint sets
    produced1 = get_mdss_vcf_path(out_prefix, 1)
    produced2 = get_mdss_vcf_path(out_prefix, 2)
    produced3 = get_mdss_vcf_path(out_prefix, 3)

    assert produced1.exists(), "First set should exist"
    assert produced2.exists(), "Second set should exist"
    assert produced3.exists(), "Third set should exist"

    # Verify sets are discriminative
    assert_discriminative(produced1, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(produced2, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(produced3, ploidy=2, min_dist=1, convert_het=False)

    # Save artifacts if enabled
    save_out_prefix_vcfs(out_prefix, subdir="unlimited")


def test_unlimited_mode_with_string(tmp_path: Path):
    """Test that unlimited mode works with 'unlimited' string."""
    samples = ["S1", "S2", "S3", "S4"]
    # Simple case with 2 possible sets
    variants = [
        # Set1: A,B
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "1/1"],
        },
        # Set2: C,D
        {
            "chrom": "1",
            "pos": 300,
            "id": "C",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "1/1", "0/0"],
        },
        {
            "chrom": "1",
            "pos": 400,
            "id": "D",
            "ref": "A",
            "alt": "T",
            "genotypes": ["1/1", "0/0", "1/1", "0/0"],
        },
    ]
    original_vcf = tmp_path / "orig_unlimited_str.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_unlimited_str"

    # Test unlimited mode with -ns unlimited
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        min_dist=1,
        n_sets="unlimited",
    )

    # Should find 2 disjoint sets
    produced1 = get_mdss_vcf_path(out_prefix, 1)
    produced2 = get_mdss_vcf_path(out_prefix, 2)

    assert produced1.exists(), "First set should exist"
    assert produced2.exists(), "Second set should exist"

    # Verify sets are discriminative
    assert_discriminative(produced1, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(produced2, ploidy=2, min_dist=1, convert_het=False)

    # Save artifacts if enabled
    save_out_prefix_vcfs(out_prefix, subdir="unlimited_str")


def test_unlimited_mode_limited_by_available_snps(tmp_path: Path):
    """Test that unlimited mode stops when no more disjoint sets can be found."""
    samples = ["S1", "S2", "S3", "S4"]
    # Only enough SNPs for 1 set
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "1/1"],
        },
    ]
    original_vcf = tmp_path / "orig_limited.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_limited"

    # Test unlimited mode
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        min_dist=1,
        n_sets=0,
    )

    # Should find only 1 set (limited by available SNPs)
    produced1 = get_mdss_vcf_path(out_prefix, 1)
    produced2 = get_mdss_vcf_path(out_prefix, 2)

    assert produced1.exists(), "First set should exist"
    assert (
        not produced2.exists()
    ), "Second set should not exist (no more SNPs available)"

    # Verify set is discriminative
    assert_discriminative(produced1, ploidy=2, min_dist=1, convert_het=False)

    # Save artifacts if enabled
    save_out_prefix_vcfs(out_prefix, subdir="limited")
