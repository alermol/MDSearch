# MDSearch
<!-- The repository contains script for the **M**inimum **D**iscriminatory SNPs set **Search** for barcoding. -->

**M**inimum **D**iscriminatory SNPs set **Search**

## Table of contents
- [How it works](#how-it-works)
- [Requrements](#requirements)
- [Installation](#installation)
- [Usage](#usage)


## How it works
MDSearch finds a minimum discriminatory set of SNPs, prioritizing them based on MAF (minor allele frequency). More suitable SNPs have MAF closer to 0.5.


MDSearch performs in two main steps:
1. Forward selection    
During this step, SNPs with the highest MAF are subsequently added to the discriminatory set until all samples become distinguishable.
2. Backward selection (minimization)    
During this step, SNPs from the primary set are removed one-by-one, and the discriminatory ability of a new set is assessed. If the discriminatory ability of the reduced set disappears, the SNP is returned to the set, and the other SNP is removed.


## Requirements
Tested using python 3.9.2.   
Additional requirements: numpy


## Installation
1. Create conda enviroment with numpy and activate it
```bash
conda create -n mdsearch numpy
conda activate mdsearch
```
2. Clone git repository
```bash
git clone https://github.com/alermol/MDSearch.git
```

3. (Optional) Run test
```bash
cd MDSearch/test
make test
```
4. MDSearch is ready to use
```bash
python3 mdsearch.py -h
```

## Usage
```
mdsearch.py [-h] [-s SEED] [-e STEPS] [-t TRIES] [-c CPU] ivcf ovcf

positional arguments:
  ivcf        input vcf file
  ovcf        output vcf file

optional arguments:
  -h, --help  show this help message and exit
  -s SEED     random seed (default: 810491)
  -e STEPS    number of backward one-by-one elimination steps (default: 10000)
  -t TRIES    number of tries to find minimal SNP set (default: 1000)
  -c CPU      number of CPUs (default: 4)
  -p          VCF is phased
```
