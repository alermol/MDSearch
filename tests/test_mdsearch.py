import sys
import subprocess
from pathlib import Path

import pytest

from .helpers import write_vcf, assert_discriminative, get_snp_ids, save_out_prefix_vcfs


def run_mdsearch(ivcf: Path, out_prefix: Path, **kwargs):
    args = [sys.executable, str(Path(__file__).resolve().parents[1] / "mdsearch.py"), str(ivcf), str(out_prefix)]
    # Map kwargs to CLI
    cli_map = {
        "seed": "-s",
        "tries": "-t",
        "cpus": "-c",
        "ploidy": "-pl",
        "total_snps": "-ts",
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


def test_ploidy_handling(tmp_path: Path):
    # Samples
    samples = ["S1", "S2", "S3", "S4"]

    # Generator: Create original VCF with one unique discriminative set of 2 SNPs (A,B)
    # and add extra non-discriminative SNPs with het and missing.
    variants = [
        # Discriminative set (unique minimal): A and B
        # Codewords across samples (diploid homozygous genotypes as 0/0 or 1/1)
        # A: S1=0, S2=0, S3=1, S4=1
        {"chrom": "1", "pos": 100, "id": "A", "ref": "A", "alt": "T", "genotypes": ["0/0", "0/0", "1/1", "1/1"]},
        # B: S1=0, S2=1, S3=0, S4=1 -> Unique patterns with A
        {"chrom": "1", "pos": 200, "id": "B", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "0/0", "1/1"]},
        # Non-discriminative C: identical to A (redundant)
        {"chrom": "1", "pos": 300, "id": "C", "ref": "A", "alt": "T", "genotypes": ["0/0", "0/0", "1/1", "1/1"]},
        # Non-discriminative D: has heterozygotes and missing
        {"chrom": "1", "pos": 400, "id": "D", "ref": "A", "alt": "T", "genotypes": ["0/1", "./.", "0/1", "./."]},
    ]

    original_vcf = tmp_path / "orig_ploidy.vcf"
    expected_vcf = tmp_path / "expected_ploidy.vcf"
    write_vcf(original_vcf, samples, variants)

    # Expected: only A and B
    write_vcf(
        expected_vcf,
        samples,
        [v for v in variants if v["id"] in {"A", "B"}],
    )

    out_prefix = tmp_path / "out_ploidy"

    # Run for diploid
    run_mdsearch(original_vcf, out_prefix, seed=1, tries=30, cpus=2, ploidy=2, total_snps=0, min_dist=1, n_sets=1)
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="ploidy")

    # Validate: Discriminative and equals expected set of IDs
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=False)
    assert set(get_snp_ids(produced)) == {"A", "B"}

    # Also test haploid behavior: build a haploid original where the same logical set works with 0 and 1 alleles
    variants_hap = [
        {"chrom": "1", "pos": 100, "id": "A", "ref": "A", "alt": "T", "genotypes": ["0", "0", "1", "1"]},
        {"chrom": "1", "pos": 200, "id": "B", "ref": "A", "alt": "T", "genotypes": ["0", "1", "0", "1"]},
        {"chrom": "1", "pos": 400, "id": "D", "ref": "A", "alt": "T", "genotypes": [".", ".", "1", "."]},
    ]
    original_vcf_hap = tmp_path / "orig_ploidy_hap.vcf"
    write_vcf(original_vcf_hap, samples, variants_hap)
    out_prefix_hap = tmp_path / "out_ploidy_hap"
    run_mdsearch(original_vcf_hap, out_prefix_hap, seed=1, tries=30, cpus=2, ploidy=1, total_snps=0, min_dist=1, n_sets=1)
    produced_hap = out_prefix_hap.with_name(out_prefix_hap.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix_hap, subdir="ploidy_hap")
    assert_discriminative(produced_hap, ploidy=1, min_dist=1, convert_het=False)
    assert set(get_snp_ids(produced_hap)) == {"A", "B"}


def test_total_snp_count_expansion(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Minimal discriminative set: A,B
    # Additional polymorphic SNPs E,F,G with varying MAF for PIC scoring
    variants = [
        {"chrom": "1", "pos": 100, "id": "A", "ref": "A", "alt": "T", "genotypes": ["0/0", "0/0", "1/1", "1/1"]},
        {"chrom": "1", "pos": 200, "id": "B", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "0/0", "1/1"]},
        # Candidate expansion SNPs
        # E has MAF ~ 0.5 (best PIC)
        {"chrom": "1", "pos": 500, "id": "E", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "1/1", "0/0"]},
        # F has MAF ~ 0.25
        {"chrom": "1", "pos": 600, "id": "F", "ref": "A", "alt": "T", "genotypes": ["0/0", "0/0", "0/0", "1/1"]},
        # G has MAF ~ 0.75 (same PIC as 0.25, comes later, but ordering follows file order)
        {"chrom": "1", "pos": 700, "id": "G", "ref": "A", "alt": "T", "genotypes": ["1/1", "1/1", "1/1", "0/0"]},
        # Add a noisy non-discriminative with het and missing
        {"chrom": "1", "pos": 800, "id": "H", "ref": "A", "alt": "T", "genotypes": ["0/1", "./.", "0/1", "./."]},
    ]
    original_vcf = tmp_path / "orig_ts.vcf"
    write_vcf(original_vcf, samples, variants)

    # Expected: minimal set A,B (only discriminative variants). We'll assert produced set is a superset and has size==total.
    expected_vcf = tmp_path / "expected_ts.vcf"
    write_vcf(expected_vcf, samples, [v for v in variants if v["id"] in {"A", "B"}])

    out_prefix = tmp_path / "out_ts"
    total = 4  # expand from 2 to 4
    run_mdsearch(original_vcf, out_prefix, seed=1, tries=30, cpus=2, ploidy=2, total_snps=total, min_dist=1, n_sets=1)
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="ts")

    # Validate size and discriminative property
    ids = get_snp_ids(produced)
    assert len(ids) == total
    assert set(["A", "B"]).issubset(ids)
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=False)


def test_min_hamming_distance_requirement(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Minimal set of 3 SNPs reaching min pairwise distance >=2 using codewords 000, 011, 101, 110
    variants = [
        {"chrom": "1", "pos": 100, "id": "X1", "ref": "A", "alt": "T", "genotypes": ["0/0", "0/0", "1/1", "1/1"]},  # bit1
        {"chrom": "1", "pos": 200, "id": "X2", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "0/0", "1/1"]},  # bit2
        {"chrom": "1", "pos": 300, "id": "X3", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "1/1", "0/0"]},  # bit3
        # Add noise with het and missing
        {"chrom": "1", "pos": 400, "id": "N1", "ref": "A", "alt": "T", "genotypes": ["0/1", "./.", "0/1", "./."]},
    ]
    original_vcf = tmp_path / "orig_md.vcf"
    write_vcf(original_vcf, samples, variants)

    expected_vcf = tmp_path / "expected_md.vcf"
    write_vcf(expected_vcf, samples, [v for v in variants if v["id"] in {"X1", "X2", "X3"}])

    out_prefix = tmp_path / "out_md"
    run_mdsearch(original_vcf, out_prefix, seed=1, tries=30, cpus=2, ploidy=2, total_snps=0, min_dist=2, n_sets=1)
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="md")

    assert_discriminative(produced, ploidy=2, min_dist=2, convert_het=False)
    assert set(get_snp_ids(produced)) == {"X1", "X2", "X3"}


def test_convert_het_handling(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Design where one SNP Y is heterozygous for several samples and thus should not contribute when -ch is used.
    # Z1 and Z2 form the minimal set when ignoring hets.
    variants = [
        # Heterozygous-heavy SNP (ignored when -ch)
        {"chrom": "1", "pos": 150, "id": "Y", "ref": "A", "alt": "T", "genotypes": ["0/1", "0/1", "1/1", "./."]},
        # Discriminative core under -ch (all homozygous patterns)
        {"chrom": "1", "pos": 210, "id": "Z1", "ref": "A", "alt": "T", "genotypes": ["0/0", "0/0", "1/1", "1/1"]},
        {"chrom": "1", "pos": 220, "id": "Z2", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "0/0", "1/1"]},
        # Extra noise
        {"chrom": "1", "pos": 230, "id": "N2", "ref": "A", "alt": "T", "genotypes": ["0/1", "./.", "0/1", "./."]},
    ]
    original_vcf = tmp_path / "orig_ch.vcf"
    write_vcf(original_vcf, samples, variants)

    expected_vcf = tmp_path / "expected_ch.vcf"
    write_vcf(expected_vcf, samples, [v for v in variants if v["id"] in {"Z1", "Z2"}])

    out_prefix = tmp_path / "out_ch"
    run_mdsearch(original_vcf, out_prefix, seed=1, tries=30, cpus=2, ploidy=2, total_snps=0, min_dist=1, convert_het=True, n_sets=1)
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="ch")

    # Distance computed ignoring het should be >=1 and set should be Z1,Z2
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=True)
    assert set(get_snp_ids(produced)) == {"Z1", "Z2"}


@pytest.mark.xfail(reason="Current algorithm may output only one minimal set; future improvement expected.")
def test_multiple_sets_generation(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Construct two alternative minimal sets:
    # Set1: A,B  and Set2: C,D
    # Patterns mirror each other to ensure both are valid minimal sets.
    variants = [
        # Set1
        {"chrom": "1", "pos": 100, "id": "A", "ref": "A", "alt": "T", "genotypes": ["0/0", "0/0", "1/1", "1/1"]},
        {"chrom": "1", "pos": 200, "id": "B", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "0/0", "1/1"]},
        # Set2
        {"chrom": "1", "pos": 300, "id": "C", "ref": "A", "alt": "T", "genotypes": ["0/0", "1/1", "1/1", "0/0"]},
        {"chrom": "1", "pos": 400, "id": "D", "ref": "A", "alt": "T", "genotypes": ["1/1", "0/0", "1/1", "0/0"]},
        # Noise
        {"chrom": "1", "pos": 500, "id": "N1", "ref": "A", "alt": "T", "genotypes": ["0/1", "./.", "0/1", "./."]},
    ]
    original_vcf = tmp_path / "orig_ns.vcf"
    write_vcf(original_vcf, samples, variants)

    # Expected sets
    expected1 = tmp_path / "expected_ns_1.vcf"
    expected2 = tmp_path / "expected_ns_2.vcf"
    write_vcf(expected1, samples, [v for v in variants if v["id"] in {"A", "B"}])
    write_vcf(expected2, samples, [v for v in variants if v["id"] in {"C", "D"}])

    out_prefix = tmp_path / "out_ns"
    # Increase tries to improve chances of discovering both sets
    run_mdsearch(original_vcf, out_prefix, seed=42, tries=200, cpus=2, ploidy=2, total_snps=0, min_dist=1, n_sets=2)

    produced1 = out_prefix.with_name(out_prefix.name + "_1.vcf")
    produced2 = out_prefix.with_name(out_prefix.name + "_2.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="ns")

    # Both files should exist and be discriminative
    assert produced1.exists()
    assert produced2.exists()
    assert_discriminative(produced1, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(produced2, ploidy=2, min_dist=1, convert_het=False)

    s1 = set(get_snp_ids(produced1))
    s2 = set(get_snp_ids(produced2))
    assert s1 in ({"A", "B"}, {"C", "D"})
    assert s2 in ({"A", "B"}, {"C", "D"})
    assert s1 != s2, "Expected two distinct discriminative sets"


