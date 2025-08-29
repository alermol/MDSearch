import subprocess
from pathlib import Path
import pytest
import sys
import json

from .helpers import write_vcf


def _run_raw(args: list[str]) -> None:
    subprocess.run(args, check=True)


def _make_minimal_vcf(tmp_path: Path) -> Path:
    samples = ["S1", "S2"]
    variants = [
        {
            "chrom": "1",
            "pos": 100,
            "id": "A",
            "ref": "A",
            "alt": "T",
            "genotypes": ["0/0", "1/1"],
        },
        {
            "chrom": "1",
            "pos": 200,
            "id": "B",
            "ref": "A",
            "alt": "T",
            "genotypes": ["1/1", "0/0"],
        },
    ]
    ivcf = tmp_path / "cli_min.vcf"
    write_vcf(ivcf, samples, variants)
    return ivcf


def _run(ivcf: Path, out_prefix: Path, **kwargs):

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


def test_negative_total_snps_removed(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_tsneg"
    # Argument removed: running without it should succeed
    _run(ivcf, out_prefix, ploidy=2, min_dist=1, n_sets=1)


def test_zero_sets_errors(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_ns0"

    with pytest.raises(subprocess.CalledProcessError):
        _run(ivcf, out_prefix, ploidy=2, min_dist=1, n_sets=0)


def test_quiet_flag_runs(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_quiet"
    _run(ivcf, out_prefix, ploidy=2, min_dist=1, n_sets=1)
    # Now with quiet flag
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
        "--quiet",
    ]
    _run_raw(args)


def test_log_level_and_format_flags_run(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_logfmt"
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
        "--log-level",
        "DEBUG",
        "--log-format",
        "json",
    ]
    _run_raw(args)


def test_invalid_log_level_choice_errors(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_badlog"
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
        "--log-level",
        "FOO",
    ]
    with pytest.raises(subprocess.CalledProcessError):
        _run_raw(args)


def test_summary_tsv_always_created(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_summary"
    _run(
        ivcf,
        out_prefix,
        ploidy=2,
        min_dist=1,
        n_sets=1,
    )
    # Summary should always be created in output directory
    summary = out_prefix / "summary.tsv"
    assert summary.exists()


# Test removed: lazy loading functionality was removed from the tool


def test_json_logging_emits_structured_lines(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_jsonlog"
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
        "--log-level",
        "INFO",
        "--log-format",
        "json",
    ]
    # Capture stderr where logs are typically sent
    proc = subprocess.run(
        args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    stderr = proc.stderr.strip().splitlines()
    # Find at least one JSON object line
    json_like = [line for line in stderr if line.startswith("{") and line.endswith("}")]
    assert json_like, f"expected JSON log lines, got: {stderr[:5]}"
    # Basic keys presence check

    obj = json.loads(json_like[0])
    assert "level" in obj and obj["level"] in {"INFO", "DEBUG", "WARNING", "ERROR"}
    assert "message" in obj


def test_text_logging_emits_bracketed_levels(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)
    out_prefix = tmp_path / "out_cli_textlog"
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
        "--log-level",
        "INFO",
        "--log-format",
        "text",
    ]
    proc = subprocess.run(
        args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    stderr = proc.stderr
    assert "[INFO]" in stderr


# New tests for format letters and fixtures
def _run_cli(ivcf: Path, out_prefix: Path, extra: list[str] | None = None) -> None:
    args = [
        sys.executable,
        str(Path(__file__).resolve().parents[1] / "mdsearch.py"),
        str(ivcf),
        str(out_prefix),
        "-pl",
        "2",
        "-md",
        "1",
        "-ns",
        "1",
    ]
    if extra:
        args.extend(extra)
    subprocess.run(args, check=True)


def test_cli_output_format_letters(tmp_path: Path):
    ivcf = _make_minimal_vcf(tmp_path)

    for letter, ext in [("v", ".vcf"), ("z", ".vcf.gz"), ("u", ".bcf"), ("b", ".bcf")]:
        out_prefix = tmp_path / f"out_{letter}"
        _run_cli(ivcf, out_prefix, ["--output-format", letter])
        produced = out_prefix / "mdss" / f"minimal_set_1{ext}"
        assert (
            produced.exists()
        ), f"expected output with extension {ext} for letter {letter}"


# def test_cli_input_format_letters_with_fixtures(tmp_path: Path):
#     fixtures = Path(__file__).resolve().parent / "fixtures"
#     cases = [
#         (fixtures / "test1_vcf_uc.vcf", "v"),
#         (fixtures / "test1_vcf_c.vcf.gz", "z"),
#         (fixtures / "test1_vcf_uc.bcf", "u"),
#         (fixtures / "test1_vcf_c.bcf", "b"),
#     ]
#     for ivcf, letter in cases:
#         out_prefix = tmp_path / f"fmt_{letter}"
#         _run_cli(ivcf, out_prefix, ["--input-format", letter, "--output-format", "v"])
#         produced = Path(f"{out_prefix}_1.vcf")
#         assert produced.exists(), f"expected text VCF output for input format {letter}"
