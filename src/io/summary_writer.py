"""TSV summary file output operations."""

from pathlib import Path
from typing import List, Union
from dataclasses import dataclass

from ..core.distance_calculator import DistanceCalculator
from ..core.vcf_parser import VCFData, LazyVCFData

__all__ = ["SetStatistics", "SummaryWriter"]


@dataclass
class SetStatistics:
    """Statistics for a SNP set."""
    set_index: int
    output_vcf: str
    num_snps: int
    min_distance: float
    snp_ids: str


class SummaryWriter:
    """Handles TSV summary file output."""
    
    def __init__(self, distance_calc: DistanceCalculator):
        self.distance_calc = distance_calc
    
    def write_summary(self, 
                     snp_sets: List[List[str]], 
                     output_prefix: Path,
                     output_path: Path,
                     vcf_data: Union[VCFData, LazyVCFData]) -> None:
        """Write per-set summary statistics to TSV."""
        header = [
            "set_index",
            "output_vcf",
            "num_snps", 
            "min_distance",
            "snp_ids",
        ]
        
        lines = ["\t".join(header)]
        
        for si, s in enumerate(snp_sets, start=1):
            out_vcf = f"{output_prefix}_{si}.vcf"
            min_d = self.distance_calc.calc_distance_for_snp_ids(s, vcf_data.snp_genotypes)
            snp_ids_str = ",".join(sorted(s))
            row = [str(si), out_vcf, str(len(s)), str(int(min_d)), snp_ids_str]
            lines.append("\t".join(row))
            
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
