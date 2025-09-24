# MDSearch: Minimum Discriminatory SNP Set Search

**MDSearch** is a Python tool designed to identify minimal sets of Single Nucleotide Polymorphisms (SNPs) that can discriminate between samples in a VCF (Variant Call Format) file based on a specified minimum Hamming distance. The tool is particularly useful for applications such as genetic barcoding, sample identification, and quality control in genomic studies.

## Key Features

- **Robust genotype handling**: Correctly processes missing genotypes and polyploid SNPs
- **Customizable discrimination**: Set minimum Hamming distance between samples
- **Multiple SNP sets**: Generate several distinct (complitely disjoint) discriminatory SNP sets
- **Advanced scoring**: Customizable entropy and Minor Allele Frequency (MAF) weights
- **Comprehensive output**: Human-readable SNP profiles, statistics, and VCF files
- **Memory monitoring**: Dynamic memory usage tracking with system-aware thresholds

## Core Techniques

The codebase implements several sophisticated algorithms:

- **Two-phase optimization**: Forward selection followed by backward elimination for minimal SNP sets
- **Vectorized Hamming distance calculation**: Optimized pairwise distance computation using NumPy arrays
- **Entropy-based SNP scoring**: Combines Shannon entropy and MAF for optimal SNP selection
- **Incremental distance updates**: Efficient backward elimination without full set reconstruction
- **Memory-aware processing**: Dynamic memory monitoring with graceful degradation
- **Deterministic SNP selection**: Stable, reproducible results with tie-breaking by SNP ID

## Key Technologies & Libraries

### Core Bioinformatics Libraries
- **[pysam](https://pysam.readthedocs.io/)**: High-performance VCF/BCF file reading and writing
- **[NumPy](https://numpy.org/)**: Vectorized numerical operations for genotype matrices
- **[bcftools](https://samtools.github.io/bcftools/)**: VCF/BCF indexing and format support

### System & Performance
- **[psutil](https://psutil.readthedocs.io/)**: Cross-platform system and process monitoring

### Quality
- **[pytest](https://docs.pytest.org/)**: Comprehensive testing framework with 29+ test cases

## Project Structure

```
MDSearch/
├── src/                              # Core application code
│   ├── core/                         # Core algorithms and data structures
│   │   ├── vcf_parser.py             # VCF file parsing and validation
│   │   ├── distance_calculator.py    # Optimized Hamming distance computation
│   │   ├── genotype_utils.py         # Genotype conversion and utility functions
│   │   └── snp_selector.py           # SNP selection and optimization algorithms
│   ├── io/                           # Input/output operations
│   │   ├── vcf_writer.py             # VCF file output with format support
│   │   ├── summary_writer.py         # Statistical summary generation
│   │   ├── run_info_writer.py        # Execution information and metadata
│   │   ├── snp_profile_writer.py     # Human-readable SNP profiles
│   │   └── structure_info_writer.py  # Output directory documentation
│   └── utils/                        # Utility functions
│       ├── memory_monitor.py         # Memory usage monitoring and warnings
│       ├── logging_setup.py          # Structured logging configuration
│       ├── indexing.py               # VCF indexing utilities
│       └── validation.py             # Input validation functions
├── tests/                            # Comprehensive test suite
├── environment.yml                   # Conda environment specification
└── mdsearch.py                       # Main command-line interface
```

**Key directories:**
- `src/core/`: Contains the core algorithms for SNP selection, distance calculation, and genotype processing
- `src/io/`: Handles all output generation including VCF files, summaries, and human-readable profiles
- `tests/`: Comprehensive test suite covering genotype handling, distance calculations, and CLI validation

## Installation & Usage

### Quick Start
```bash
git clone https://github.com/alermol/MDSearch.git
cd MDSearch
conda env create -f environment.yml
conda activate mdsearch
python mdsearch.py input.vcf results -md 2 -ns 3
```

### Basic Command
```bash
python mdsearch.py input.vcf output_folder -pl 2 -md 2 -ns 3
```

This will:
- Assume that input VCF file `input.vcf` is diploid
- Ensure that minimum Hamming distance of 2 between samples in all final discriminatory VCF sets
- Generate 3 distinct discriminatory SNP sets
- Create comprehensive output in the `output_folder`

## Command Line Interface

### Basic Syntax
```bash
python mdsearch.py IVCF OUTPUT_FOLDER [options]
```

### Positional Arguments
- **IVCF**: Input VCF file (bi-allelic, with SNP IDs)
- **OUTPUT_FOLDER**: Path to output folder (will be created if absent)

### Selection Parameters
- `-pl, --ploidy PLOIDY` (default: 2) — VCF ploidy
- `-md, --min-distance MIN_DIST` (default: 1) — Minimal Hamming distance between samples
- `-ns, --num-sets N_SETS` (default: 0) — Number of distinct SNP sets in output (0 or 'unlimited' for all possible sets)
- `-we, --weight-entropy WEIGHT_ENTROPY` (default: 0.5) — Weight for entropy component in SNP scoring (0.0 to 1.0)
- `-wm, --weight-maf WEIGHT_MAF` (default: 0.5) — Weight for MAF component in SNP scoring (0.0 to 1.0)

### Output Options
- `-ch, --convert-het` — Convert heterozygous calls into NA

### Input/Output Formats
- `-I, --input-format {auto,v,z,u,b}` (default: auto) — Input format (VCF/VCF.gz/BCF uncompressed/BCF compressed)
- `-O, --output-format {v,z,u,b}` (default: v) — Output format (VCF/VCF.gz/BCF uncompressed/BCF compressed)

### Logging Options
- `-q, --quiet` — Suppress progress output
- `-L, --log-level {DEBUG,INFO,WARNING,ERROR}` — Logging level
- `-F, --log-format {text,json}` (default: text) — Logging format

### Information
- `-V, --version` — Show program version, commit hash, and Python version and exit
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

#### Unlimited Mode (Find All Possible Disjoint SNP Sets)
```bash
python mdsearch.py data.vcf results -pl 2 -md 2 -ch -ns 0
# or
python mdsearch.py data.vcf results -pl 2 -md 2 -ch -ns unlimited
```

#### Advanced Usage with Logging
```bash
python mdsearch.py data.vcf results -md 2 -ns 3 -L INFO
```

#### JSON Logging for Automated Pipelines
```bash
python mdsearch.py data.vcf results -F json -q
```

#### BCF Format Processing
```bash
python mdsearch.py input.bcf results -I b -O b
```

## Output Structure

The tool creates a well-organized output directory with comprehensive analysis results:

```
results/
├── mdss/                           # SNP set VCF files
│   ├── minimal_set_1.vcf           # First minimal discriminative set
│   ├── minimal_set_2.vcf           # Alternative minimal set
│   └── minimal_set_3.vcf           # Third alternative set
├── snp_profiles/                   # Human-readable SNP profiles
│   ├── set_1_profile.txt           # Profile for set 1
│   ├── set_2_profile.txt           # Profile for set 2
│   └── set_3_profile.txt           # Profile for set 3
├── summary.tsv                     # Comprehensive statistics for each SNP set
├── run_info.txt                    # Detailed execution information
├── output_structure.txt            # Output organization guide
├── best_set.vcf                    # Optimal set with highest entropy
└── best_set_profile.txt            # Profile for the optimal SNP set
```

### Detailed File Descriptions

#### `summary.tsv` - Statistical Summary
A tab-separated values file containing comprehensive statistics for each generated SNP set:

| Column | Description |
|--------|-------------|
| `set_index` | 1-based index of the SNP set (1, 2, 3, ...) |
| `output_vcf` | Path to the corresponding VCF file (e.g., "minimal_set_1.vcf") |
| `num_snps` | Number of SNPs in the set |
| `min_distance` | Achieved minimum Hamming distance between all sample pairs |
| `shannon_entropy` | Shannon entropy score measuring chromosome distribution balance |

The Shannon entropy indicates how well SNPs are distributed across chromosomes - higher values indicate more balanced distribution, which is often desirable for genetic studies.

#### `snp_profiles/` - Human-Readable Profiles
Each profile file (`set_N_profile.txt`) provides a compact, human-readable representation of SNP sets:

**Header Section:**
- SNP set identification and chromosome distribution summary
- Format: `#SNP Set N Profile` and `#Chromosome Distribution: Chr1(2); Chr2(3); Chr3(1)`

**Sample Profiles:**
Each sample is represented as a line showing all SNP genotypes in the format:
```
Sample_name: Chr1:Position1(Genotype1) Chr2:Position2(Genotype2) ...
```

**Genotype Format:**
- **Diploid**: `AA`, `AT`, `TT` (homozygous reference, heterozygous, homozygous alternate)
- **Triploid**: `AAA`, `AAT`, `TTT` (no separators between alleles)
- **Missing**: `--` (for diploid) or `---` (for triploid)
- **Polyploid**: Concatenated alleles without separators

**Example Profile Content:**
```
#SNP Set 1 Profile
#Chromosome Distribution: Chr1(2); Chr2(1); Chr3(2)
--------------------------------------------------
Sample_001: Chr1:12345(AA) Chr1:67890(TT) Chr2:11111(AT) Chr3:22222(AA) Chr3:33333(TT)
Sample_002: Chr1:12345(AT) Chr1:67890(AT) Chr2:11111(AA) Chr3:22222(AT) Chr3:33333(AT)
```

#### `best_set.vcf` and `best_set_profile.txt`
The optimal SNP set automatically selected based on highest Shannon entropy:
- `best_set.vcf`: Standard VCF format of the optimal SNP set
- `best_set_profile.txt`: Human-readable profile of the optimal set (same format as other profiles)

#### `run_info.txt` - Execution Metadata
Comprehensive execution information including:
- Program version and build details
- Input file specifications and parameters
- Output statistics and timing information
- System memory usage and thresholds
- Complete configuration used for the run

#### `output_structure.txt` - Directory Guide
Detailed explanation of all output files and their purposes, serving as a reference for understanding the results.

## Advanced Features

### Custom SNP Scoring
```bash
python mdsearch.py input.vcf results -we 0.8 -wm 0.2
```
Control the balance between entropy (`-we`) and MAF (`-wm`) components in SNP selection.

### Unlimited Mode
```bash
python mdsearch.py input.vcf results -ns unlimited
```
Find all possible disjoint discriminatory SNP sets.

### Structured Logging
```bash
python mdsearch.py input.vcf results -F json -q
```
Generate JSON-formatted logs for automated pipeline integration.

## Testing

Run the comprehensive test suite:
```bash
pytest -v
```

The test suite covers genotype handling, distance calculations, CLI validation, and VCF format compliance.

## License

This project is open-source and available under the MIT License.

---

**Note**: MDSearch today represents a major advancement over the version described in the [original publication](https://doi.org/10.3390/agronomy14081802), featuring extensive improvements and new capabilities. If you use MDSearch in your work, please cite the original article.
