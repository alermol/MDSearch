# **MDSearch**

**MDSearch today is a major advancement over MDSearch described in [this article](https://doi.org/10.3390/agronomy14081802), featuring extensive improvements and new capabilities—so much so that it can be considered a distinct tool. If you use MDSearch in your work, please cite the original article linked above.**


**MDSearch** (**M**inimum **D**iscriminatory SNPs set **Search**) is a Python tool designed to identify minimal sets of Single Nucleotide Polymorphisms (SNPs) that can discriminate between samples in a VCF (Variant Call Format) file based on a specified minimum Hamming distance. The tool is particularly useful for applications such as genetic barcoding, sample identification, and quality control in genomic studies.

## Description
**Key features**
* Correctly handles **missing genotypes** and **polyploid SNPs**
* Customizable **minimal Hamming distance** between samples
* Selects **core discriminatory SNP sets**
* Generates **multiple distinct SNP sets** for validation
* Flexible handling of **heterozygous calls** (convertible to NA)
* **Structured logging** with text/JSON output formats
* **Comprehensive output** with human-readable SNP profiles, run information, and statistics
* **Advanced SNP scoring** with customizable entropy and MAF weights
* **Unlimited mode** to find all possible disjoint SNP sets
* **Best set detection** based on Shannon entropy for optimal chromosome distribution

## Table of Contents
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Notes](#notes)
- [Testing](#testing)
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
   Sequentially adds SNPs with highest combined entropy and MAF scores until all samples are distinguishable at the specified minimal Hamming distance. The scoring system balances information content (entropy) with allele frequency distribution (MAF) for optimal SNP selection.

2. **Backward Elimination Phase**  
   Iteratively removes SNPs from the preliminary set while maintaining the required discriminative ability. Uses deterministic elimination to find minimal sets.

The algorithm finds minimal discriminating sets only; expansion to a user-specified size is not performed.

## Requirements

### Dependencies
- **Python 3.8+** - Required for type hints and modern features
- **NumPy** - For numerical operations and array handling
- **pysam** - For VCF/BCF file reading and writing
- **psutil** - For memory monitoring and system information
- **bcftools** - For VCF/BCF indexing (required for compressed formats)
- **pytest** - For running the test suite

### Standard Python Libraries
- pathlib, logging, time, json, sys, itertools, math, dataclasses, typing

## Installation

### Conda Environment (Recommended)
```bash
git clone https://github.com/alermol/MDSearch.git
cd MDSearch
conda env create -f environment.yml
conda activate mdsearch
```

## Usage
### Basic Command
```bash
python mdsearch.py input.vcf output_folder
```

### CLI Reference (grouped)

Positional arguments:
- **IVCF**: Input VCF file (bi-allelic, with SNP IDs)
- **OUTPUT_FOLDER**: Path to output folder (will be created if absent)

Selection:
- `-pl PLOIDY` (default: 2) — VCF ploidy
- `-md MIN_DIST` (default: 1) — Minimal Hamming distance between samples
- `-ns N_SETS` (default: 1) — Number of distinct SNP sets in output (0 or 'unlimited' for all possible sets)
- `-we WEIGHT_ENTROPY` (default: 0.5) — Weight for entropy component in SNP scoring (0.0 to 1.0)
- `-wm WEIGHT_MAF` (default: 0.5) — Weight for MAF component in SNP scoring (0.0 to 1.0)

Output:
- `-ch` — Convert heterozygous calls into NA

IO formats:
- `--input-format {auto,v,z,u,b}` (default: auto) — Input format (VCF/VCF.gz/BCF uncompressed/BCF compressed)
- `--output-format {v,z,u,b}` (default: v) — Output format (VCF/VCF.gz/BCF uncompressed/BCF compressed)

Logging:
- `--quiet` — Suppress progress output
- `--log-level {DEBUG,INFO,WARNING,ERROR}` — Logging level
- `--log-format {text,json}` (default: text) — Logging format

Misc:
- `--version` — Show program version, commit hash, and Python version and exit
- `-h, --help` — Show help message and exit

### Examples

#### Basic Usage
```bash
python mdsearch.py data.vcf results -pl 2 -md 2 -ch -ns 3
```

#### Custom SNP Scoring Weights
```bash
python mdsearch.py data.vcf results -pl 2 -md 2 -we 0.8 -wm 0.2
```

#### Unlimited Mode (Find All Possible Sets)
```bash
python mdsearch.py data.vcf results -pl 2 -md 2 -ch -ns 0
# or
python mdsearch.py data.vcf results -pl 2 -md 2 -ch -ns unlimited
```

#### Advanced Usage with Logging
```bash
python mdsearch.py data.vcf results -md 2 -ns 3 --log-level INFO
```

#### JSON Logging for Automated Pipelines
```bash
python mdsearch.py data.vcf results --log-format json --quiet
```

#### BCF Format Processing
```bash
python mdsearch.py input.bcf results --input-format b --output-format b
```

The first command will:
- Assume a ploidy of 2
- Generate minimal discriminating set(s)
- Ensure a minimum Hamming distance of 2 between samples
- Convert heterozygous calls to NA
- Generate 3 distinct SNP sets

## Output
Creates a comprehensive output folder structure:

```
results/
├── mdss/                           # SNP set VCF files
│   ├── minimal_set_1.vcf           # First minimal discriminative set
│   ├── minimal_set_2.vcf           # Alternative minimal set
│   ├── minimal_set_3.vcf           # Third alternative set    
├── snp_profiles/                   # Human-readable SNP profiles
│   ├── set_1_profile.txt           # Profile for set 1
│   ├── set_2_profile.txt           # Profile for set 2
│   └── set_3_profile.txt           # Profile for set 3
├── summary.tsv                     # Comprehensive statistics for each SNP set
├── run_info.txt                    # Detailed execution information
├── output_structure.txt            # Output organization guide
├── best_set.vcf                    # Optimal set with highest entropy
└── best_set_profile.txt            # Human-readable format profile for the optimal SNP set
```


### Summary TSV Contents
The summary file includes columns:
- `set_index`: Set number (1, 2, 3...)
- `output_vcf`: Path to corresponding VCF file
- `num_snps`: Number of SNPs in the set
- `min_distance`: Achieved minimum Hamming distance
- `shannon_entropy`: Chromosome distribution balance score

## Notes
- Input VCF must contain SNP IDs in the third column (ID field)
- Use `-ch` when heterozygous calls shouldn't contribute to discrimination
- Typical minimal distances:
  - `-md 1` for basic discrimination
  - `-md 2-3` for more robust discrimination
- **Logging options:**
  - `--quiet` suppresses progress output
  - `--log-level` controls verbosity (DEBUG, INFO, WARNING, ERROR)
  - `--log-format json` outputs structured logs for automated pipelines
- **SNP scoring:**
  - `-we` controls entropy weight (information content)
  - `-wm` controls MAF weight (allele frequency distribution)
  - Combined scoring provides more informative SNP sets
- **File format support:**
  - Input formats: auto-detection, VCF (v), compressed VCF (z), uncompressed BCF (u), compressed BCF (b)
  - Output formats: VCF (v), compressed VCF (z), uncompressed BCF (u), compressed BCF (b)
- **Best set selection:**
  - Automatically identifies the SNP set with highest Shannon entropy
  - Provides optimal chromosome distribution balance

## Testing
Run the test suite (make sure the `mdsearch` environment is active):
```bash
# Run all tests
pytest -v

# Run tests without saving VCF artifacts  
MDSEARCH_SAVE_VCFS=0 pytest -v
```

The comprehensive test suite covers:
- Different genotype ploidy handling (including haploid)
- Enforcing minimal Hamming distance (`-md`)
- Correct handling of heterozygous genotypes in distance calculation (`-ch`)
- Multiple distinct discriminative SNP sets (`-ns`)
- Phased and unphased VCF handling
- Multi-field FORMAT parsing (GT:DP:GQ)
- TSV summary generation

## License
This project is open-source and available under the MIT License.
