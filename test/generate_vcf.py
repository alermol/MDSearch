import random
from pathlib import Path
from itertools import compress
from string import ascii_lowercase


def parser_resolve_path(path):
    return Path(path).resolve()


def generate_sample_names(n: int):
    number = [
        [int(j) for j in list(bin(i)[2:].zfill(len(ascii_lowercase)))]
        for i in random.sample(range(50000, 999999), k=n)
    ]
    return "\t".join("".join(list(compress(ascii_lowercase, i))) for i in number)


def main(ovcf, nsamples, nsnp, ploidy, phased):
    with open(ovcf, "w") as ov:
        ov.write("##fileformat=VCFv4.3\n")
        ov.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + generate_sample_names(n=nsamples)
            + "\n"
        )
        for i in range(nsnp):
            sep = "|" if phased else "/"
            hom_gen = [(f"{str(i)}{sep}" * ploidy).strip(sep) for i in range(ploidy)]
            ov.write(
                f"1\t{i + 1}\tSNP{i}\tA\tG\t60\tPASS\t.\tGT\t"
                + "\t".join(
                    random.choices(hom_gen + [(f".{sep}" * ploidy).strip(sep)], k=nsnp)
                )
                + "\n"
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("ovcf", help="output vcf file", type=parser_resolve_path)
    parser.add_argument("--ploidy", help="required ploidy", type=int, default=2)
    parser.add_argument("--seed", help="random seed", default=625530192, type=int)
    parser.add_argument("--nsamples", help="number of samples", default=200, type=int)
    parser.add_argument("--nsnp", help="number of snp", default=1000, type=int)
    parser.add_argument("--phased", help="generate phased VCF", action="store_true")

    args = parser.parse_args()

    random.seed(args.seed)
    main(
        ovcf=args.ovcf,
        nsamples=args.nsamples,
        nsnp=args.nsnp,
        ploidy=args.ploidy,
        phased=args.phased,
    )
