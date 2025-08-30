"""SNP Profile Writer module for generating human-readable SNP set formats."""

import logging
from pathlib import Path
from typing import Dict, List, Optional

from ..core.vcf_parser import VCFData


class SNPProfileWriter:
    """Handles generation of human-readable SNP profiles from SNP sets."""

    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize SNP Profile Writer.

        Args:
            logger: Optional logger instance for output messages
        """
        self.logger = logger or logging.getLogger(__name__)

    def write_snp_profiles(
        self,
        output_prefix: Path,
        snp_sets: List[List[str]],
        vcf_data: VCFData,
        ploidy: int = 2,
        chr_recode_file: Optional[Path] = None,
        sample_recode_file: Optional[Path] = None,
    ) -> Path:
        """Generate human-readable SNP profiles for all SNP sets.

        Args:
            output_prefix: Path to output directory
            snp_sets: List of SNP sets, each containing SNP IDs
            vcf_data: Parsed VCF data containing SNP and sample information
            chr_recode_file: Optional chromosome recoding file (TSV)
            sample_recode_file: Optional sample recoding file (TSV)

        Returns:
            Path to the snp_profiles directory

        Example:
            >>> writer = SNPProfileWriter()
            >>> profiles_dir = writer.write_snp_profiles(
            ...     Path("output"), [["rs1", "rs2"], ["rs3", "rs4"]], vcf_data
            ... )
        """
        # Create snp_profiles subdirectory
        profiles_dir = output_prefix / "snp_profiles"
        profiles_dir.mkdir(parents=True, exist_ok=True)

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Generating SNP profiles in: {profiles_dir}")

        # Apply recoding if files provided
        sample_names = vcf_data.headers.samples.copy()
        if sample_recode_file is not None:
            sample_names = self._recode_samples(sample_names, sample_recode_file)

        # Generate profile for each SNP set
        for set_index, snp_set in enumerate(snp_sets, 1):
            profile_file = profiles_dir / f"snp_set_{set_index}_profile.txt"
            self._write_single_profile(
                profile_file,
                snp_set,
                vcf_data,
                sample_names,
                set_index,
                ploidy,
                chr_recode_file,
            )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Generated {len(snp_sets)} SNP profile files")

        return profiles_dir

    def write_best_set_profile(
        self,
        output_prefix: Path,
        best_set_vcf_path: Path,
        vcf_data: VCFData,
        chr_recode_file: Optional[Path] = None,
        sample_recode_file: Optional[Path] = None,
    ) -> Path:
        """Create human-readable profile for the best SNP set from best_set.vcf.

        Args:
            output_prefix: Output directory prefix
            best_set_vcf_path: Path to best_set.vcf file
            vcf_data: Parsed VCF data
            chr_recode_file: Optional chromosome recoding file
            sample_recode_file: Optional sample recoding file

        Returns:
            Path to the best set profile file
        """
        if not best_set_vcf_path.exists():
            if self.logger.isEnabledFor(logging.WARNING):
                self.logger.warning(f"Best set VCF not found: {best_set_vcf_path}")
            return output_prefix

        # Parse the best set VCF to extract SNP IDs
        snp_ids = []
        try:
            with open(best_set_vcf_path) as vcf_file:
                for line in vcf_file:
                    if not line.startswith("#"):
                        parts = line.strip().split("\t")
                        if len(parts) >= 3:
                            snp_ids.append(parts[2])  # SNP ID is in column 3
        except Exception as e:
            if self.logger.isEnabledFor(logging.ERROR):
                self.logger.error(f"Failed to parse best set VCF: {e}")
            return output_prefix

        if not snp_ids:
            if self.logger.isEnabledFor(logging.WARNING):
                self.logger.warning("No SNP IDs found in best set VCF")
            return output_prefix

        # Apply recoding if files provided
        sample_names = vcf_data.headers.samples.copy()
        if sample_recode_file is not None:
            sample_names = self._recode_samples(sample_names, sample_recode_file)

        # Generate profile for the best set directly in output folder
        profile_file = output_prefix / "best_set_profile.txt"
        self._write_single_profile(
            profile_file,
            snp_ids,
            vcf_data,
            sample_names,
            0,  # Use 0 to indicate best set
            vcf_data.headers.ploidy if hasattr(vcf_data.headers, "ploidy") else 2,
            chr_recode_file,
        )

        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(f"Generated best set profile: {profile_file}")

        return profile_file

    def _write_single_profile(
        self,
        profile_file: Path,
        snp_set: List[str],
        vcf_data: VCFData,
        sample_names: List[str],
        set_index: int,
        ploidy: int,
        chr_recode_file: Optional[Path] = None,
    ) -> None:
        """Write a single SNP set profile to file.

        Args:
            profile_file: Output file path for this profile
            snp_set: List of SNP IDs in this set
            vcf_data: Parsed VCF data
            sample_names: List of sample names (potentially recoded)
            set_index: Index number of this SNP set
            ploidy: Ploidy level for genotype interpretation
            chr_recode_file: Optional chromosome recoding file
        """
        # Load chromosome recoding if provided
        chr_recode_scheme = {}
        if chr_recode_file is not None:
            chr_recode_scheme = self._load_recode_scheme(chr_recode_file)

        with open(profile_file, "w") as output:
            # Write compact header with chromosome distribution only
            chr_counts: Dict[str, int] = {}

            for snp_id in snp_set:
                if snp_id in vcf_data.snp_genotypes:
                    snp_data = vcf_data.snp_genotypes[snp_id]
                    chrom = snp_data.chromosome
                    if chr_recode_scheme:
                        chrom = chr_recode_scheme.get(chrom, chrom)
                    chr_counts[chrom] = chr_counts.get(chrom, 0) + 1

            # Write chromosome distribution
            chr_dist = "; ".join(
                [f"{chrom}({count})" for chrom, count in sorted(chr_counts.items())]
            )

            # Add title based on whether this is the best set or a regular SNP set
            if set_index == 0:
                output.write("#Best SNP Set Profile\n")
                output.write(f"#Chromosome Distribution: {chr_dist}\n")
            else:
                output.write(f"#SNP Set {set_index} Profile\n")
                output.write(f"#Chromosome Distribution: {chr_dist}\n")

            output.write("-" * 50 + "\n")

            # Write sample profiles

            for sample_idx, sample_name in enumerate(sample_names):
                sample_genotypes = []

                for snp_id in snp_set:
                    if snp_id in vcf_data.snp_genotypes:
                        snp_data = vcf_data.snp_genotypes[snp_id]
                        chrom = snp_data.chromosome
                        if chr_recode_scheme:
                            chrom = chr_recode_scheme.get(chrom, chrom)

                        # Get the original GT string from sample_fields
                        gt_string = snp_data.sample_fields[sample_idx]

                        if gt_string == "." or gt_string == "":
                            # Missing genotype - repeat "-" by ploidy
                            genotype_str = "-" * ploidy
                        else:
                            # Parse the GT string to get allele indices
                            gt_parts = gt_string.replace("|", "/").split("/")
                            allele_indices = []

                            for part in gt_parts:
                                if part == "." or part == "":
                                    allele_indices.append(None)
                                else:
                                    try:
                                        allele_indices.append(int(part))
                                    except ValueError:
                                        allele_indices.append(None)

                            # Convert allele indices to nucleotides
                            nucleotides = []
                            for allele_idx in allele_indices:
                                if allele_idx is None:
                                    nucleotides.append("-")
                                elif allele_idx == 0:
                                    nucleotides.append(snp_data.reference_allele)
                                elif allele_idx == 1:
                                    nucleotides.append(snp_data.alternate_allele)
                                else:
                                    # Handle other alleles if they exist
                                    nucleotides.append(str(allele_idx))

                            # Join nucleotides based on ploidy
                            if ploidy == 1:
                                genotype_str = nucleotides[0] if nucleotides else "-"
                            else:
                                # For polyploid, concatenate without separator
                                genotype_str = "".join(nucleotides[:ploidy])

                        sample_genotypes.append(
                            f"{snp_data.chromosome}:{snp_data.position}({genotype_str})"
                        )

                output.write(f"{sample_name}: {'; '.join(sample_genotypes)}\n")

    def _load_recode_scheme(self, recode_file: Path) -> Dict[str, str]:
        """Load recoding scheme from TSV file.

        Args:
            recode_file: TSV file with old_name -> new_name mappings

        Returns:
            Dictionary mapping old names to new names
        """
        recode_scheme = {}
        try:
            with open(recode_file) as file:
                for line in file:
                    line_parts = line.strip().split("\t")
                    if len(line_parts) >= 2:
                        recode_scheme[line_parts[0]] = line_parts[1]
        except Exception as e:
            if self.logger.isEnabledFor(logging.WARNING):
                self.logger.warning(f"Failed to load recode file {recode_file}: {e}")

        return recode_scheme

    def _recode_samples(self, sample_names: List[str], recode_file: Path) -> List[str]:
        """Rename sample names using a two-column TSV mapping file.

        Args:
            sample_names: Original sample names
            recode_file: TSV file with old_name -> new_name mappings

        Returns:
            List with recoded sample names
        """
        recode_scheme = self._load_recode_scheme(recode_file)
        recoded_names = []

        for sample in sample_names:
            recoded_names.append(recode_scheme.get(sample, sample))

        return recoded_names
