# MDSearch
**MDSearch** is a tool for a **M**inimum **D**iscriminatory SNPs set **Search**


**Main features** of MDSearch:
* MDSearch correctly handles **missing genotypes**
* MDSearch allows you to customize **minimal distance** between samples in the selected SNPs set
* MDSearch correctly handles **polyploid SNPs**
* MDSearch allows you **to select the core SNPs set** (minimum discriminatory set), as well as **complement it** with additional polymorphic SNPs to the required amount
* MDSearch allows you to select **the required number** of the best (minimal) SNP sets


## Table of contents
- [How it works](#how-it-works)
- [Requrements](#requirements)
- [Installation](#installation)
- [Usage](#usage)


## How it works
MDSearch finds a minimum discriminatory set of SNPs of given samples in VCF file, prioritizing them based on MAF (minor allele frequency). More suitable (polymorphic) SNPs have MAF closer to 0.5.

MDSearch performs in two main steps:
1. ***Forward selection***    
During this step, SNPs with the highest MAF are subsequently added to the discriminatory set until all samples become distinguishable with required minimal distance between samples.
2. ***Backward selection (minimization)***    
During this step, SNPs from the primary set are removed one-by-one, and the discriminatory ability of a new set is assessed. If the discriminatory ability of the reduced set (with required minimal distance between samples) disappears, the SNP is returned to the set, and the other SNP is removed.

Final MDS, if required, complemets with additional polymorphic SNPs (based on PIC) to the required size.


## Requirements
File requirements.txt contains a list of all required packages for MDSearch.


## Installation
1. Clone git repository
```bash
git clone https://github.com/alermol/MDSearch.git
```

2. Create conda enviroment and activate it
```bash
cd MDSearch
conda create --name mdsearch --file requirements.txt
conda activate mdsearch
```

3. MDSearch is ready to use
```bash
python3 mdsearch.py -h
```

## Usage
```bash
usage: mdsearch.py [-h] [-s SEED] [-e STEPS] [-t TRIES] [-c CPU] [-pl PLOIDY] 
                   [-ts TOTAL SNP] [-md MIN DIST] [-ch] [-ns N SETS] ivcf ovcf_prefix

positional arguments:
  ivcf           input vcf file
  ovcf_prefix    prefix of output vcf file

options:
  -h, --help     show this help message and exit
  -s SEED        random seed (default: 810491)
  -e STEPS       number of backward one-by-one elimination steps (default: 10000)
  -t TRIES       number of tries to find minimal SNP set (default: 1000)
  -c CPU         number of CPUs (default: 4)
  -pl PLOIDY     VCF ploidy (default: 2)
  -ts TOTAL SNP  Total number of SNPs in output set (Default: minimal discriminative set)
  -md MIN DIST   Minimal hamming distance between samples (Default: 1)
  -ch            Convert heterozygous calls into NA
  -ns N SETS     Number of distinct SNP sets in output (Default: 1)
```
