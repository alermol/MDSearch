import sys
import subprocess
from pathlib import Path


from .helpers import write_vcf, assert_discriminative, get_snp_ids, save_out_prefix_vcfs


def run_mdsearch(ivcf: Path, out_prefix: Path, **kwargs):
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
    ]
    # Map kwargs to CLI
    cli_map = {
        "ploidy": "-pl",
        "total_snps": "-ts",
        "min_dist": "-md",
        "convert_het": "-ch",
        "n_sets": "-ns",
        "overlap_max_number": "-oMx",
        "overlap_max_fraction": "-oMf",
        "summary_tsv": "--summary-tsv",
    }
    for k, v in kwargs.items():
        flag = cli_map[k]
        if isinstance(v, bool):
            if v:
                args.append(flag)
        else:
            args.extend([flag, str(v)])
    subprocess.run(args, check=True)


def test_ploidy_handling_diploid_and_haploid(tmp_path: Path):
    # Samples
    samples = ["S1", "S2", "S3", "S4"]

    # Generator: Create original VCF with one unique discriminative set of 2 SNPs (A,B)
    # and add extra non-discriminative SNPs with het and missing.
    variants = [
        # Discriminative set (unique minimal): A and B
        # Codewords across samples (diploid homozygous genotypes as 0/0 or 1/1)
        # A: S1=0, S2=0, S3=1, S4=1
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },
        # B: S1=0, S2=1, S3=0, S4=1 -> Unique patterns with A
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "1/1"],
        },
        # Non-discriminative C: identical to A (redundant)
        {
            "chrom": "1",
            "pos": 300,
            "id": "C",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },
        # Non-discriminative D: has heterozygotes and missing
        {
            "chrom": "1",
            "pos": 400,
            "id": "D",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
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
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="ploidy")

    # Validate: Discriminative and equals expected set of IDs
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=False)
    assert set(get_snp_ids(produced)) == {"A", "B"}

    # Also test haploid behavior: build a haploid original where the same logical set works with 0 and 1 alleles
    variants_hap = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0", "0", "1", "1"],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0", "1", "0", "1"],
        },
        {
            "chrom": "1",
            "pos": 400,
            "id": "D",
            "ref": "A",
            "alt": "T",
            "genotypes": [".", ".", "1", "."],
        },
    ]
    original_vcf_hap = tmp_path / "orig_ploidy_hap.vcf"
    write_vcf(original_vcf_hap, samples, variants_hap)
    out_prefix_hap = tmp_path / "out_ploidy_hap"
    run_mdsearch(
        original_vcf_hap,
        out_prefix_hap,
        ploidy=1,
        total_snps=0,
        min_dist=1,
        n_sets=1,
    )
    produced_hap = out_prefix_hap.with_name(out_prefix_hap.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix_hap, subdir="ploidy_hap")
    assert_discriminative(produced_hap, ploidy=1, min_dist=1, convert_het=False)
    assert set(get_snp_ids(produced_hap)) == {"A", "B"}


def test_total_snp_count_expansion_adds_by_pic(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Minimal discriminative set: A,B
    # Additional polymorphic SNPs E,F,G with varying MAF for PIC scoring
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
        # Candidate expansion SNPs
        # E has MAF ~ 0.5 (best PIC)
        {
            "chrom": "1",
            "pos": 500,
            "id": "E",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "1/1", "0/0"],
        },
        # F has MAF ~ 0.25
        {
            "chrom": "1",
            "pos": 600,
            "id": "F",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "0/0", "1/1"],
        },
        # G has MAF ~ 0.75 (same PIC as 0.25, comes later, but ordering follows file order)
        {
            "chrom": "1",
            "pos": 700,
            "id": "G",
            "ref": "A",
            "alt": "T",
            "genotypes": ["1/1", "1/1", "1/1", "0/0"],
        },
        # Add a noisy non-discriminative with het and missing
        {
            "chrom": "1",
            "pos": 800,
            "id": "H",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
    ]
    original_vcf = tmp_path / "orig_ts.vcf"
    write_vcf(original_vcf, samples, variants)

    # Expected: minimal set A,B (only discriminative variants). We'll assert produced set is a superset and has size==total.
    expected_vcf = tmp_path / "expected_ts.vcf"
    write_vcf(expected_vcf, samples, [v for v in variants if v["id"] in {"A", "B"}])

    out_prefix = tmp_path / "out_ts"
    total = 4  # expand from 2 to 4
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=total,
        min_dist=1,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="ts")

    # Validate size and discriminative property
    ids = get_snp_ids(produced)
    assert len(ids) == total
    assert set(["A", "B"]).issubset(ids)
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=False)


def test_min_hamming_distance_requirement_all_pairs(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Minimal set of 3 SNPs reaching min pairwise distance >=2 using codewords 000, 011, 101, 110
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "X1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },  # bit1
        {
            "chrom": "1",
            "pos": 200,
            "id": "X2",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "1/1"],
        },  # bit2
        {
            "chrom": "1",
            "pos": 300,
            "id": "X3",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "1/1", "0/0"],
        },  # bit3
        # Add noise with het and missing
        {
            "chrom": "1",
            "pos": 400,
            "id": "N1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
    ]
    original_vcf = tmp_path / "orig_md.vcf"
    write_vcf(original_vcf, samples, variants)

    expected_vcf = tmp_path / "expected_md.vcf"
    write_vcf(
        expected_vcf, samples, [v for v in variants if v["id"] in {"X1", "X2", "X3"}]
    )

    out_prefix = tmp_path / "out_md"
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=0,
        min_dist=2,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="md")

    assert_discriminative(produced, ploidy=2, min_dist=2, convert_het=False)
    assert set(get_snp_ids(produced)) == {"X1", "X2", "X3"}


def test_convert_het_handling_ignores_heterozygotes(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Design where one SNP Y is heterozygous for several samples and thus should not contribute when -ch is used.
    # Z1 and Z2 form the minimal set when ignoring hets.
    variants = [
        # Heterozygous-heavy SNP (ignored when -ch)
        {
            "chrom": "1",
            "pos": 150,
            "id": "Y",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "0/1", "1/1", "./."],
        },
        # Discriminative core under -ch (all homozygous patterns)
        {
            "chrom": "1",
            "pos": 210,
            "id": "Z1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 220,
            "id": "Z2",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "1/1"],
        },
        # Extra noise
        {
            "chrom": "1",
            "pos": 230,
            "id": "N2",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
    ]
    original_vcf = tmp_path / "orig_ch.vcf"
    write_vcf(original_vcf, samples, variants)

    expected_vcf = tmp_path / "expected_ch.vcf"
    write_vcf(expected_vcf, samples, [v for v in variants if v["id"] in {"Z1", "Z2"}])

    out_prefix = tmp_path / "out_ch"
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        convert_het=True,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="ch")

    # Distance computed ignoring het should be >=1 and set should be Z1,Z2
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=True)
    assert set(get_snp_ids(produced)) == {"Z1", "Z2"}


def test_multiple_sets_generation_no_overlap(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Construct two alternative minimal sets:
    # Set1: A,B  and Set2: C,D
    # Patterns mirror each other to ensure both are valid minimal sets.
    variants = [
        # Set1
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
        # Set2
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
        # Noise
        {
            "chrom": "1",
            "pos": 500,
            "id": "N1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
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
    # Require no overlap between sets (disjoint sets)
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=2,
        overlap_max_number=0,
    )

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
    # Expect disjoint discriminative pairs
    assert s1.isdisjoint(s2)
    assert s1 != s2, "Expected two distinct discriminative sets"


def test_multiple_sets_with_max_overlap_number_cap(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Design three informative SNPs where multiple minimal pairs exist
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
        {
            "chrom": "1",
            "pos": 300,
            "id": "C",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "1/1", "0/0"],
        },
        # Noise with het/missing
        {
            "chrom": "1",
            "pos": 400,
            "id": "N1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
    ]
    original_vcf = tmp_path / "orig_ovl.vcf"
    write_vcf(original_vcf, samples, variants)

    # Without overlap constraint, we can get two distinct pairs
    out_prefix1 = tmp_path / "out_ovl1"
    run_mdsearch(
        original_vcf,
        out_prefix1,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=2,
    )
    p1a = out_prefix1.with_name(out_prefix1.name + "_1.vcf")
    p1b = out_prefix1.with_name(out_prefix1.name + "_2.vcf")
    save_out_prefix_vcfs(out_prefix1, subdir="ovl_none")
    s1 = set(get_snp_ids(p1a))
    s2 = set(get_snp_ids(p1b))
    assert s1 != s2
    assert_discriminative(p1a, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(p1b, ploidy=2, min_dist=1, convert_het=False)

    # With overlap_max_number=1, ensure sets share at most one SNP
    out_prefix2 = tmp_path / "out_ovl2"
    run_mdsearch(
        original_vcf,
        out_prefix2,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=2,
        overlap_max_number=1,
    )
    p2a = out_prefix2.with_name(out_prefix2.name + "_1.vcf")
    p2b = out_prefix2.with_name(out_prefix2.name + "_2.vcf")
    save_out_prefix_vcfs(out_prefix2, subdir="ovl_n1")
    s1b = set(get_snp_ids(p2a))
    s2b = set(get_snp_ids(p2b))
    base_overlap = len(s1b.intersection(s2b))
    assert base_overlap <= 1
    assert s1b != s2b
    assert_discriminative(p2a, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(p2b, ploidy=2, min_dist=1, convert_het=False)


def test_multiple_sets_with_max_overlap_fraction_cap(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
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
            "id": "N1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
    ]
    original_vcf = tmp_path / "orig_ovlf.vcf"
    write_vcf(original_vcf, samples, variants)

    # No cap (default): ensure two distinct sets exist
    out_prefix_f1 = tmp_path / "out_ovlf1"
    run_mdsearch(
        original_vcf,
        out_prefix_f1,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=2,
    )
    pfa = out_prefix_f1.with_name(out_prefix_f1.name + "_1.vcf")
    pfb = out_prefix_f1.with_name(out_prefix_f1.name + "_2.vcf")
    save_out_prefix_vcfs(out_prefix_f1, subdir="ovl_fnc")
    assert pfa.exists() and pfb.exists()
    sfa = set(get_snp_ids(pfa))
    sfb = set(get_snp_ids(pfb))
    assert sfa != sfb
    assert_discriminative(pfa, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(pfb, ploidy=2, min_dist=1, convert_het=False)

    # With overlap_max_fraction=0.5, intersection must be <= floor(0.5 * |base|) = 1 (since base size is 2)
    out_prefix_fm = tmp_path / "out_ovlfm"
    run_mdsearch(
        original_vcf,
        out_prefix_fm,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=2,
        overlap_max_fraction=0.5,
    )
    pfm_a = out_prefix_fm.with_name(out_prefix_fm.name + "_1.vcf")
    pfm_b = out_prefix_fm.with_name(out_prefix_fm.name + "_2.vcf")
    save_out_prefix_vcfs(out_prefix_fm, subdir="ovl_f05")
    assert pfm_a.exists() and pfm_b.exists()
    sfm_a = set(get_snp_ids(pfm_a))
    sfm_b = set(get_snp_ids(pfm_b))
    assert len(sfm_a.intersection(sfm_b)) <= 1
    assert_discriminative(pfm_a, ploidy=2, min_dist=1, convert_het=False)
    assert_discriminative(pfm_b, ploidy=2, min_dist=1, convert_het=False)

    # With overlap_max_fraction=0.0, disjoint requirement leaves only one set in this toy design
    out_prefix_f2 = tmp_path / "out_ovlf2"
    run_mdsearch(
        original_vcf,
        out_prefix_f2,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=2,
        overlap_max_fraction=0.0,
    )
    pf2a = out_prefix_f2.with_name(out_prefix_f2.name + "_1.vcf")
    pf2b = out_prefix_f2.with_name(out_prefix_f2.name + "_2.vcf")
    save_out_prefix_vcfs(out_prefix_f2, subdir="ovl_f00")
    assert pf2a.exists()
    assert not pf2b.exists()
    assert_discriminative(pf2a, ploidy=2, min_dist=1, convert_het=False)


def test_phased_vcf_basic_pipe_separator_handled(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0|0", "0|0", "1|1", "1|1"],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0|0", "1|1", "0|0", "1|1"],
        },
        # Noise phased het/missing
        {
            "chrom": "1",
            "pos": 300,
            "id": "N1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0|1", ".|.", "0|1", ".|."],
        },
    ]
    original_vcf = tmp_path / "orig_phased.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_phased"
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="phased_basic")

    assert produced.exists()
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=False)
    assert set(get_snp_ids(produced)) == {"A", "B"}


def test_phased_vcf_with_ch_converts_hets(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Core discriminators are Z1 and Z2 (both homozygous). P contains phased hets.
    variants = [
        {
            "chrom": "1",
            "pos": 210,
            "id": "Z1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0|0", "0|0", "1|1", "1|1"],
        },
        {
            "chrom": "1",
            "pos": 220,
            "id": "Z2",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0|0", "1|1", "0|0", "1|1"],
        },
        {
            "chrom": "1",
            "pos": 150,
            "id": "P",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0|1", "0|1", "1|1", "0|0"],
        },
        # Extra phased noise
        {
            "chrom": "1",
            "pos": 230,
            "id": "N2",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0|1", ".|.", "0|1", ".|."],
        },
    ]
    original_vcf = tmp_path / "orig_phased_ch.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_phased_ch"
    # Force inclusion of P via total_snps expansion
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=3,
        min_dist=1,
        convert_het=True,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="phased_ch")

    assert produced.exists()
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=True)
    ids = set(get_snp_ids(produced))
    assert {"Z1", "Z2"}.issubset(ids)
    assert "P" in ids  # P added during expansion
    # Heterozygous phased calls must be converted to '.|.' in output for selected P
    txt = Path(produced).read_text()
    assert ".|." in txt


def test_unphased_vcf_with_ch_converts_hets(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Core discriminators are Z1 and Z2; Q contains unphased hets
    variants = [
        {
            "chrom": "1",
            "pos": 210,
            "id": "Z1",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "0/0", "1/1", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 220,
            "id": "Z2",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1", "0/0", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 150,
            "id": "Q",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "0/1", "1/1", "0/0"],
        },
    ]
    original_vcf = tmp_path / "orig_unphased_ch.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_unphased_ch"
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=3,
        min_dist=1,
        convert_het=True,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    save_out_prefix_vcfs(out_prefix, subdir="unphased_ch")

    assert produced.exists()
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=True)
    txt = Path(produced).read_text()
    assert "./." in txt


def test_multiallelic_site_causes_error_exit(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Include a genotype with allele index 2 to trigger multiallelic error
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "M1",
            "ref": "A",
            "alt": "T,C",
            "genotypes": ["0/0", "0/0", "1/2", "2/2"],
        },
    ]
    original_vcf = tmp_path / "orig_multiallelic.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_multiallelic"
    import pytest

    with pytest.raises(subprocess.CalledProcessError):
        run_mdsearch(
            original_vcf,
            out_prefix,
            ploidy=2,
            total_snps=0,
            min_dist=1,
            n_sets=1,
        )


def test_empty_vcf_causes_error_exit(tmp_path: Path):
    samples = ["S1", "S2"]
    # Write VCF with headers only, no variant lines
    empty_vcf = tmp_path / "orig_empty.vcf"
    write_vcf(empty_vcf, samples, [])

    out_prefix = tmp_path / "out_empty"
    import pytest

    with pytest.raises(subprocess.CalledProcessError):
        run_mdsearch(
            empty_vcf,
            out_prefix,
            ploidy=2,
            total_snps=0,
            min_dist=1,
            n_sets=1,
        )


def test_single_sample_vcf_causes_error_exit(tmp_path: Path):
    # Single-sample header should fail during header validation (need >= 2 samples)
    samples = ["S1"]
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0"],
        }
    ]
    single_vcf = tmp_path / "orig_single.vcf"
    write_vcf(single_vcf, samples, variants)

    out_prefix = tmp_path / "out_single"
    import pytest

    with pytest.raises(subprocess.CalledProcessError):
        run_mdsearch(
            single_vcf,
            out_prefix,
            ploidy=2,
            total_snps=0,
            min_dist=1,
            n_sets=1,
        )


def test_total_snps_lower_than_minimal_kept_at_min(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
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
        {
            "chrom": "1",
            "pos": 300,
            "id": "N",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
    ]
    original_vcf = tmp_path / "orig_tslt.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_tslt"
    # Request fewer SNPs than minimal size (1 < 2). Expect still 2 SNPs.
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=2,
        total_snps=1,
        min_dist=1,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    ids = get_snp_ids(produced)
    assert len(ids) == 2
    assert set(ids) == {"A", "B"}


def test_overlap_flags_mutually_exclusive_error(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
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
    original_vcf = tmp_path / "orig_ovl_err.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_ovl_err"
    import pytest

    with pytest.raises(subprocess.CalledProcessError):
        run_mdsearch(
            original_vcf,
            out_prefix,
            ploidy=2,
            total_snps=0,
            min_dist=1,
            n_sets=2,
            overlap_max_number=1,
            overlap_max_fraction=0.5,
        )


def test_triploid_ploidy_handling_with_mixed_genotypes(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Triploid: use 3 alleles per sample, with only 0 and 1 allele indices
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0/0", "0/0/0", "1/1/1", "1/1/1"],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0/0", "1/1/1", "0/0/0", "1/1/1"],
        },
        # Heterozygous triploid
        {
            "chrom": "1",
            "pos": 300,
            "id": "H",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1/1", "././.", "0/1/1", "././."],
        },
    ]
    original_vcf = tmp_path / "orig_pl3.vcf"
    write_vcf(original_vcf, samples, variants)

    out_prefix = tmp_path / "out_pl3"
    run_mdsearch(
        original_vcf,
        out_prefix,
        ploidy=3,
        total_snps=0,
        min_dist=1,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    assert produced.exists()
    assert_discriminative(produced, ploidy=3, min_dist=1, convert_het=False)


def _write_multifield_vcf(path: Path, samples: list[str], variants: list[dict]):
    """Write a VCF with FORMAT=GT:DP:GQ and sample fields accordingly."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=mdsearch-tests\n")
        f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples)
            + "\n"
        )
        for v in variants:
            fmt = "GT:DP:GQ"
            sample_cols = [
                f"{gt}:{dp}:{gq}" for gt, dp, gq in zip(v["gts"], v["dps"], v["gqs"])
            ]
            line = [
                v.get("chrom", "1"),
                str(v.get("pos", 1)),
                v["id"],
                v.get("ref", "A"),
                v.get("alt", "T"),
                ".",
                "PASS",
                ".",
                fmt,
            ] + sample_cols
            f.write("\t".join(line) + "\n")


def test_gt_only_parsing_with_multifield_format(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "gts": ["0/0", "0/0", "1/1", "1/1"],
            "dps": [5, 7, 9, 11],
            "gqs": [50, 60, 70, 80],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "gts": ["0/0", "1/1", "0/0", "1/1"],
            "dps": [6, 8, 10, 12],
            "gqs": [55, 65, 75, 85],
        },
        # Noise with het and missing GTs but with DP/GQ present
        {
            "chrom": "1",
            "pos": 300,
            "id": "N",
            "ref": "A",
            "alt": "T",
            "gts": ["0/1", "./.", "0/1", "./."],
            "dps": [3, 0, 4, 0],
            "gqs": [20, 0, 25, 0],
        },
    ]
    orig = tmp_path / "orig_multifield.vcf"
    _write_multifield_vcf(orig, samples, variants)
    out_prefix = tmp_path / "out_multifield"
    run_mdsearch(orig, out_prefix, ploidy=2, total_snps=0, min_dist=1, n_sets=1)
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    assert produced.exists()
    assert_discriminative(produced, ploidy=2, min_dist=1, convert_het=False)
    assert set(get_snp_ids(produced)) == {"A", "B"}


def test_writer_preserves_format_and_replaces_only_gt_on_ch(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
    # Core discriminators Z1,Z2 homozygous; P has hets and should be included via expansion
    variants = [
        {
            "chrom": "1",
            "pos": 210,
            "id": "Z1",
            "ref": "A",
            "alt": "T",
            "gts": ["0/0", "0/0", "1/1", "1/1"],
            "dps": [10, 11, 12, 13],
            "gqs": [90, 91, 92, 93],
        },
        {
            "chrom": "1",
            "pos": 220,
            "id": "Z2",
            "ref": "A",
            "alt": "T",
            "gts": ["0/0", "1/1", "0/0", "1/1"],
            "dps": [14, 15, 16, 17],
            "gqs": [94, 95, 96, 97],
        },
        {
            "chrom": "1",
            "pos": 150,
            "id": "P",
            "ref": "A",
            "alt": "T",
            "gts": ["0/1", "0/1", "1/1", "0/0"],
            "dps": [21, 22, 23, 24],
            "gqs": [70, 71, 72, 73],
        },
    ]
    orig = tmp_path / "orig_multifield_ch.vcf"
    _write_multifield_vcf(orig, samples, variants)
    out_prefix = tmp_path / "out_multifield_ch"
    run_mdsearch(
        orig,
        out_prefix,
        ploidy=2,
        total_snps=3,
        min_dist=1,
        convert_het=True,
        n_sets=1,
    )
    produced = out_prefix.with_name(out_prefix.name + "_1.vcf")
    assert produced.exists()
    # Validate DP/GQ preserved while GT for hets becomes ./.
    text = Path(produced).read_text().splitlines()
    records = [line for line in text if not line.startswith("#")]
    # Find record P
    p = [line for line in records if line.split("\t")[2] == "P"][0]
    cols = p.split("\t")
    assert cols[8] == "GT:DP:GQ"
    # Sample 1 and 2 were hets in input → should be converted to ./.
    s1 = cols[9].split(":")
    s2 = cols[10].split(":")
    assert s1[0] == "./." and int(s1[1]) == 21 and int(s1[2]) == 70
    assert s2[0] == "./." and int(s2[1]) == 22 and int(s2[2]) == 71
    # For a homozygous sample (S3), GT stays and DP/GQ preserved
    s3 = cols[11].split(":")
    assert s3[0] == "1/1" and int(s3[1]) == 23 and int(s3[2]) == 72


def test_summary_tsv_emitted_and_correct(tmp_path: Path):
    samples = ["S1", "S2", "S3", "S4"]
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
        {
            "chrom": "1",
            "pos": 300,
            "id": "N",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/1", "./.", "0/1", "./."],
        },
    ]
    orig = tmp_path / "orig_summary.vcf"
    write_vcf(orig, samples, variants)
    out_prefix = tmp_path / "out_summary"
    summary_path = tmp_path / "summary.tsv"
    run_mdsearch(
        orig,
        out_prefix,
        ploidy=2,
        total_snps=0,
        min_dist=1,
        n_sets=1,
        summary_tsv=summary_path,
    )
    # Save artifacts if enabled
    save_out_prefix_vcfs(out_prefix, subdir="summary")
    # Validate TSV exists
    assert summary_path.exists()
    # Read and validate content
    lines = summary_path.read_text().strip().splitlines()
    assert lines[0].split("\t") == [
        "set_index",
        "output_vcf",
        "num_snps",
        "min_distance",
        "snp_ids",
    ]
    # Only one set expected; minimal set is A,B with min_distance >=1
    cols = lines[1].split("\t")
    assert cols[0] == "1"
    assert cols[1].endswith("out_summary_1.vcf")
    assert cols[2] == "2"
    assert int(cols[3]) >= 1
    ids = set(cols[4].split(","))
    assert ids == {"A", "B"}
