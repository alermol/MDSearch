# **MDSearch**

**MDSearch** (**M**inimum **D**iscriminatory SNPs set **Search**) is a Python script designed to identify minimal sets of Single Nucleotide Polymorphisms (SNPs) that can discriminate between samples in a VCF (Variant Call Format) file based on a specified minimum Hamming distance. The tool is particularly useful for applications such as genetic barcoding, sample identification, and quality control in genomic studies.

## Description
**Key features**
* Correctly handles **missing genotypes** and **polyploid SNPs**
* Customizable **minimal Hamming distance** between samples
* Selects **core discriminatory SNP sets** and can **expand them** with additional polymorphic SNPs
* Generates **multiple distinct SNP sets** for validation
* **Overlap constraints** for alternative set generation (`-oMx`, `-oMf`)
* Flexible handling of **heterozygous calls** (convertible to NA)
* **Structured logging** with text/JSON output formats
* **Optional TSV summaries** with per-set statistics


## Table of Contents
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Notes](#notes)
- [License](#license)


## How it works

⚠️ **Important Notes**    
SNPs in the input VCF must be properly filtered before being processed with MDSearch. Typical steps include:
- Removing or splitting multiallelic sites.
- Removing SNPs with at least one missing genotype (optional, as MDSearch correctly handles missing data).
- Removing SNPs with MAF (Minor Allele Frequency) ≤ 10%.
- Removing SNPs with excessive heterozygosity.
- Performing LD pruning.

MDSearch employs a two-phase optimization approach:

1. **Forward Selection Phase**  
   Sequentially adds SNPs with highest Minor Allele Frequency (MAF) until all samples are distinguishable at the specified minimal Hamming distance. SNPs with MAF closest to 0.5 are prioritized for maximum discriminative power.

2. **Backward Elimination Phase**  
   Iteratively removes SNPs from the preliminary set while maintaining the required discriminative ability. Uses deterministic elimination to find minimal sets.

The final step optionally expands the minimal set with additional polymorphic SNPs (using Polymorphism Information Content, PIC) to reach a user-specified size.


## Requirements
- Python 3.8+
- NumPy
- Standard Python libraries (pathlib, logging, time, json, sys, itertools, math)


## Installation
```bash
git clone https://github.com/alermol/MDSearch.git
cd MDSearch
conda create --name mdsearch python numpy pytest
conda activate mdsearch
```

## Usage
### Basic Command
```bash
python mdsearch.py input.vcf output_prefix
```

### Advanced Options
```bash
usage: mdsearch.py [-h] [-pl PLOIDY] [-ts TOTAL_SNP] [-md MIN_DIST] [-ch] [-ns N_SETS]
                   [-oMx OVERLAP_MAX_N] [-oMf OVERLAP_MAX_FRAC] [--quiet]
                   [--log-level {DEBUG,INFO,WARNING,ERROR}] [--log-format {text,json}]
                   [--summary-tsv SUMMARY_TSV]
                   IVCF OVCF_PREFIX

positional arguments:
  IVCF                  Input VCF file (bi-allelic, with SNP IDs)
  OVCF_PREFIX           Prefix for output VCF(s)

options:
  -h, --help            Show help message and exit
  -pl PLOIDY            VCF ploidy (default: 2)
  -ts TOTAL_SNP         Total SNPs in output set (0 = keep minimal; default: 0)
  -md MIN_DIST          Minimal Hamming distance between samples (default: 1)
  -ch                   Convert heterozygous calls into NA (default: False)
  -ns N_SETS            Number of distinct SNP sets in output (default: 1)
  -oMx OVERLAP_MAX_N    Maximum overlap count for alternative sets (-1 = unlimited)
  -oMf OVERLAP_MAX_FRAC Maximum overlap fraction for alternative sets (-1 = unlimited)
  --quiet               Suppress progress output (default: False)
  --log-level {DEBUG,INFO,WARNING,ERROR}
                        Logging level (default depends on --quiet)
  --log-format {text,json}
                        Logging format: text or json (default: text)
  --summary-tsv SUMMARY_TSV
                        Write per-set summary TSV to specified path
```

### Examples

#### Basic Usage
```bash
python mdsearch.py data.vcf results -pl 2 -ts 50 -md 2 -ch -ns 3
```

#### Advanced Usage with Logging and Summaries
```bash
python mdsearch.py data.vcf results -md 2 -ns 3 -oMx 1 --log-level INFO --summary-tsv results_summary.tsv
```

#### JSON Logging for Automated Pipelines
```bash
python mdsearch.py data.vcf results --log-format json --quiet --summary-tsv results.tsv
```

The first command will:
- Assume a ploidy of 2
- Expand the final SNP set to 50 SNPs
- Ensure a minimum Hamming distance of 2 between samples
- Convert heterozygous calls to NA
- Generate 3 distinct SNP sets

## Output
Generates one or more VCF files containing only selected discriminative SNPs:
- `results_1.vcf` (Minimal discriminative set)
- `results_2.vcf` (Alternative minimal set)
- `results_3.vcf` (Third alternative set)

Each file maintains original VCF format while including only selected SNPs.

### Optional TSV Summary
When `--summary-tsv` is specified, generates a summary file with columns:
- `set_index`: Set number (1, 2, 3...)
- `output_vcf`: Path to corresponding VCF file
- `num_snps`: Number of SNPs in the set
- `min_distance`: Achieved minimum Hamming distance
- `snp_ids`: Comma-separated list of selected SNP IDs

Only `results_1.vcf` will be generated with `-ns 1` option.

## Notes
- Input VCF must contain SNP IDs in the third column (ID field)
- Use `-ch` when heterozygous calls shouldn't contribute to discrimination
- Typical minimal distances:
  - `-md 1` for basic discrimination
  - `-md 2-3` for more robust discrimination
- **Overlap constraints for multiple sets:**
  - `-oMx N` limits alternative sets to at most N SNPs in common with the base set
  - `-oMf F` limits alternative sets to at most fraction F overlap with the base set
  - These options are mutually exclusive
- **Logging options:**
  - `--quiet` suppresses progress output
  - `--log-level` controls verbosity (DEBUG, INFO, WARNING, ERROR)
  - `--log-format json` outputs structured logs for automated pipelines
- The script `generate_snp_passport.py` in the `scripts` folder will generate a human-readable passport in the following format:
  > sample_name: Chr:Coordinate(SNP_genotype); Chr:Coordinate(SNP_genotype); Chr:Coordinate(SNP_genotype)

## Testing
Run the test suite (make sure the `mdsearch` environment is active):
```bash
# Run all tests
pytest -v

# Run tests without saving VCF artifacts  
MDSEARCH_SAVE_VCFS=0 pytest -v
```

The comprehensive test suite (17 tests) covers:
- Different genotype ploidy handling (including haploid)
- Expansion to a target total SNP count (`-ts`)
- Enforcing minimal Hamming distance (`-md`)
- Correct handling of heterozygous genotypes in distance calculation (`-ch`)
- Multiple distinct discriminative SNP sets (`-ns`)
- Overlap constraints (`-oMx`, `-oMf`)
- Phased and unphased VCF handling
- Multi-field FORMAT parsing (GT:DP:GQ)
- Multiallelic site detection and rejection
- TSV summary generation

## License
This project is open-source and available under the MIT License.
