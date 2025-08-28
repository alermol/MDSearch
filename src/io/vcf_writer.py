"""VCF file output operations."""

from pathlib import Path
from typing import List, TextIO
from dataclasses import dataclass

from ..core.genotype_utils import extract_gt, is_het

__all__ = ["WriteConfig", "VCFWriter"]


@dataclass
class WriteConfig:
    """Configuration for VCF writing."""
    ploidy: int
    convert_het: bool


class VCFWriter:
    """Handles VCF file output operations."""
    
    def write_snp_sets(self, 
                      input_vcf: Path,
                      output_prefix: Path, 
                      snp_sets: List[List[str]],
                      config: WriteConfig) -> None:
        """Write each SNP set to separate VCF files."""
        for si, s in enumerate(snp_sets, start=1):
            output_file = f"{output_prefix}_{si}.vcf"
            
            with open(input_vcf) as invcf, open(output_file, "w") as outvcf:
                for vcf_line in invcf:
                    if vcf_line.startswith("#"):
                        outvcf.write(vcf_line)
                    else:
                        line = vcf_line.strip().split("\t")
                        if (line[2] in s) and config.convert_het:
                            self._write_line_with_het_conversion(line, outvcf, config)
                        elif (line[2] in s) and (not config.convert_het):
                            outvcf.write(vcf_line)
                        else:
                            continue
                            
    def _write_line_with_het_conversion(self, 
                                       line: List[str], 
                                       outvcf: TextIO, 
                                       config: WriteConfig) -> None:
        """Write VCF line with heterozygous calls converted to missing."""
        format_field = line[8] if len(line) > 8 else "GT"
        keys = format_field.split(":") if format_field else []
        gt_index = keys.index("GT") if "GT" in keys else -1
        
        # Determine separator from first sample's GT
        first_gt = None
        if gt_index >= 0 and len(line) > 9:
            first_parts = line[9].split(":")
            if gt_index < len(first_parts):
                first_gt = first_parts[gt_index]
        sep = "/" if (first_gt and "/" in first_gt) else "|"
        missing_gt = (
            sep.join(["."] * config.ploidy)
            if (config.ploidy and config.ploidy > 1)
            else "."
        )
        
        new_samples = []
        for sample_field in line[9:]:
            if gt_index == -1:
                new_samples.append(sample_field)
                continue
                
            parts = sample_field.split(":")
            gt_val = (
                parts[gt_index] if gt_index < len(parts) else "."
            )
            
            if is_het(gt_val):
                if gt_index < len(parts):
                    parts[gt_index] = missing_gt
                    new_samples.append(":".join(parts))
                else:
                    new_samples.append(sample_field)
            else:
                new_samples.append(sample_field)
                
        line = line[:9] + new_samples
        outvcf.write("\t".join(line) + "\n")
