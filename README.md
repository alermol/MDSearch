# **MDSearch**

**MDSearch** (**M**inimum **D**iscriminatory SNPs set **Search**) is a Python script designed to identify minimal sets of Single Nucleotide Polymorphisms (SNPs) that can discriminate between samples in a VCF (Variant Call Format) file based on a specified minimum Hamming distance. The tool is particularly useful for applications such as genetic barcoding, sample identification, and quality control in genomic studies.

## Description
**Key features**
* Correctly handles **missing genotypes** and **polyploid SNPs**
* Customizable **minimal Hamming distance** between samples
* Selects **core discriminatory SNP sets** and can **expand them** with additional polymorphic SNPs
* Generates **multiple distinct SNP sets** for validation
* Efficient **parallel processing** with multi-threading support
* Flexible handling of **heterozygous calls** (convertible to NA)


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
   Iteratively removes SNPs from the preliminary set while maintaining the required discriminative ability. Uses random elimination with multiple trials to find minimal sets.

The final step optionally expands the minimal set with additional polymorphic SNPs (using Polymorphism Information Content, PIC) to reach a user-specified size.


## Requirements
- Python 3.x
- NumPy
- Standard Python libraries (multiprocessing, re, random, pathlib)


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
usage: mdsearch.py [-h] [-pl PLOIDY] [-ts TOTAL SNP]
                   [-md MIN DIST] [-ch] [-ns N SETS] ivcf ovcf_prefix

positional arguments:
  ivcf           input vcf file
  ovcf_prefix    prefix of output vcf file

options:
  -h, --help     Show help message
  -pl PLOIDY     Sample ploidy (default: 2)
  -ts TOTAL SNP  Total SNPs in output (minimal set)
  -md MIN DIST   Minimal Hamming distance (default: 1)
  -ch            Convert heterozygotes to NA
  -ns N SETS     Number of output sets (default: 1)
```

### Example
```bash
python mdsearch.py data.vcf results -pl 2 -ts 50 -md 2 -ch -ns 3
```

This command will:
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

Only `results_1.vcf` will be generated with `-ns 1` option.

## Notes
- Input VCF must contain SNP IDs in the third column (ID field)
- Use `-ch` when heterozygous calls shouldn't contribute to discrimination
- Typical minimal distances:
  - `-md 1` for basic discrimination
  - `-md 2-3` for more robust discrimination
- The script `generate_snp_passport.py` in the `scripts` folder will generate a human-readable passport in the following format:
  > sample_name: Chr:Coordinate(SNP_genotype); Chr:Coordinate(SNP_genotype); Chr:Coordinate(SNP_genotype)

## Testing
Run the test suite (make sure the `mdsearch` environment is active):
```bash
pytest -q
```

The tests generate synthetic VCFs exercising:
- Different genotype ploidy handling (`--pl`)
- Expansion to a target total SNP count (`-ts`)
- Enforcing minimal Hamming distance (`-md`)
- Correct handling of heterozygous genotypes in distance calculation (`-ch`)
- Producing multiple distinct discriminative SNP sets (`-ns`)

## License
This project is open-source and available under the MIT License.
