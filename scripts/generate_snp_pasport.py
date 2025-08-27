from itertools import repeat
from pathlib import Path
import re
from typing import Dict, List


def parser_resolve_path(path: str) -> Path:
    """Resolve CLI-provided path string to an absolute Path."""
    return Path(path).resolve()


def vcf_to_json_transpose(vcf: Path) -> Dict[str, List[List[str]]]:
    """Transpose VCF samples into a per-sample list of variant records."""
    with open(vcf) as file:
        for line in file:
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                break
    res = {i: [] for i in samples}
    with open(vcf) as file:
        for line in file:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split("\t")
                genotypes = tuple(
                    i[0] + [i[-1]] for i in zip(repeat(line[:2] + line[3:5]), line[9:])
                )
                geno_samples = {i: j for i, j in zip(samples, genotypes)}
                for k, v in geno_samples.items():
                    res[k].append(v)
    return res


def recode_samples(
    genotypes: Dict[str, List[List[str]]], recode_file: Path
) -> Dict[str, List[List[str]]]:
    """Rename sample keys using a two-column TSV mapping file."""
    recode_scheme = {}
    with open(recode_file) as file:
        for line in file:
            line = line.strip().split("\t")
            recode_scheme[line[0]] = line[1]
    result = {recode_scheme[k]: v for k, v in genotypes.items()}
    return result


def recode_chr(
    genotypes: Dict[str, List[List[str]]], recode_file: Path
) -> Dict[str, List[List[str]]]:
    """Recode chromosome names in each record using a two-column TSV mapping."""
    recode_scheme = {}
    with open(recode_file) as file:
        for line in file:
            line = line.strip().split("\t")
            recode_scheme[line[0]] = line[1]
    result = {
        k: [[recode_scheme[i[0]]] + i[1:] for i in v] for k, v in genotypes.items()
    }
    return result


def recode_genotypes(
    genotypes: Dict[str, List[List[str]]],
) -> Dict[str, List[List[str]]]:
    """Expand numeric genotype to allele strings; replace missing with dashes."""
    result = {}
    for k, v in genotypes.items():
        for snp in v:
            # remove separators '/' and '|' from genotype string
            genotype = re.sub(r"[\/|]", "", snp[-1])
            if "." in genotype:
                snp[-1] = "-" * len(genotype)
            else:
                ref, alt = snp[-3], snp[-2]
                genotype = genotype.replace("0", ref).replace("1", alt)
                snp[-1] = genotype
        result[k] = v
    return result


def main(
    invcf: Path,
    outfile: Path,
    chr_recode_file: Path | None,
    sample_recode_file: Path | None,
) -> None:
    """Generate a simple SNP passport text file from a VCF."""
    print("Start working...")
    genotypes = vcf_to_json_transpose(invcf)
    if sample_recode_file is not None:
        genotypes = recode_samples(genotypes, sample_recode_file)
    if chr_recode_file is not None:
        genotypes = recode_chr(genotypes, chr_recode_file)
    genotypes = recode_genotypes(genotypes)

    with open(outfile, "w") as output:
        for sample, genotype in genotypes.items():
            output.write(f"{sample}: ")
            output.write(
                "; ".join([f"{i[0]}:{i[1]}({i[-1]})" for i in genotype]) + "\n"
            )
    print(f"Done. Output file: {outfile}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("ivcf", help="input vcf file", type=parser_resolve_path)
    parser.add_argument("ofile", help="output file", type=parser_resolve_path)
    parser.add_argument(
        "-sr",
        help="sample recode file",
        default=None,
        type=parser_resolve_path,
        metavar="SRF",
    )
    parser.add_argument(
        "-cr",
        help="chromosome recode file",
        default=None,
        type=parser_resolve_path,
        metavar="CRF",
    )

    args = parser.parse_args()

    main(args.ivcf, args.ofile, args.cr, args.sr)
