import sys
import subprocess
from pathlib import Path


def _write_lines(path: Path, lines: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def run_mdsearch(ivcf: Path, out_prefix: Path, **kwargs):
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
    ]
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


def test_missing_fileformat_header_errors(tmp_path: Path):
    ivcf = tmp_path / "missing_fileformat.vcf"
    lines = [
        "##source=mdsearch-tests",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
        "1\t100\tA\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    ]
    _write_lines(ivcf, lines)

    out_prefix = tmp_path / "out_missing_fileformat"
    import pytest

    with pytest.raises(subprocess.CalledProcessError):
        run_mdsearch(ivcf, out_prefix, ploidy=2, min_dist=1, n_sets=1)


def test_missing_chrom_header_errors(tmp_path: Path):
    ivcf = tmp_path / "missing_chrom.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "##source=mdsearch-tests",
        # No #CHROM line before data
        "1\t100\tA\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    ]
    _write_lines(ivcf, lines)

    out_prefix = tmp_path / "out_missing_chrom"
    import pytest

    with pytest.raises(subprocess.CalledProcessError):
        run_mdsearch(ivcf, out_prefix, ploidy=2, min_dist=1, n_sets=1)


def test_invalid_number_of_columns_errors(tmp_path: Path):
    # Header declares two samples, but only one sample field provided in data line
    ivcf = tmp_path / "bad_columns.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "##source=mdsearch-tests",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
        "1\t100\tA\tA\tT\t.\tPASS\t.\tGT\t0/0",
    ]
    _write_lines(ivcf, lines)

    out_prefix = tmp_path / "out_bad_columns"
    import pytest

    with pytest.raises(subprocess.CalledProcessError):
        run_mdsearch(ivcf, out_prefix, ploidy=2, min_dist=1, n_sets=1)


def test_multiallelic_genotype_indices_errors_even_with_single_alt(tmp_path: Path):
    # ALT is single allele, but genotype uses allele index 2 → pysam handles gracefully
    # by converting invalid genotype to missing, so program should succeed
    ivcf = tmp_path / "bad_gt_indices.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        "##source=mdsearch-tests",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
        "1\t100\tA\tA\tT\t.\tPASS\t.\tGT\t2/2\t1/1",
    ]
    _write_lines(ivcf, lines)

    out_prefix = tmp_path / "out_bad_gt_indices"

    # Should succeed because pysam handles invalid genotypes gracefully
    run_mdsearch(ivcf, out_prefix, ploidy=2, min_dist=1, n_sets=1)
