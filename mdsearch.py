import random
import numpy as np
from multiprocessing import Pool
from pathlib import Path

import sys


def parser_resolve_path(path):
    return Path(path).resolve()


class MDSearch:
    def __init__(
        self,
        in_vcf,
        out_vcf,
        seed=None,
        elimination_steps=None,
        tries=None,
        ncups=None,
        ploidy=None,
    ):
        random.seed(seed)
        self.in_vcf = in_vcf
        self.out_vcf = out_vcf
        self.elimination_steps = elimination_steps
        self.tries = tries
        self.ncups = ncups
        self.ploidy = ploidy

        # calculate target number of genotypes and create list containing genotype for each SNP
        self.snp_genotypes = {}
        with open(self.in_vcf) as vcf:
            for l in vcf.readlines():
                l = l.strip()
                if l.startswith("#CHROM"):
                    self.target_gen_n = len(l.split("\t")[9:])
                elif l.startswith("#"):
                    continue
                else:
                    snp_id = l.split("\t")[2]
                    geno = l.split("\t")[9:]
                    self.snp_genotypes[snp_id] = geno
        self.main()

    @staticmethod
    def _calculate_maf(geno: list, ploidy: int):
        geno = "".join(geno)
        total_alleles = len(geno)
        maf = 1 / ploidy
        for i in range(ploidy):
            if geno.count(str(i)) == 0:
                continue
            else:
                if (geno.count(str(i)) / total_alleles) < maf:
                    maf = geno.count(str(i)) / total_alleles
                else:
                    continue
        return maf

    @staticmethod
    def _calc_N_distinct(snps: list):
        g = np.unique(np.array([i for i in snps]), axis=1)
        return g.shape[1]

    def select_first_snp(self):
        # select SNP with max MAF
        snp_id, snp_maf = None, None
        for sid, sg in self.snp_genotypes.items():
            if (snp_id is None) and (snp_maf is None):
                snp_id = sid
                snp_maf = self._calculate_maf(sg, self.ploidy)
            else:
                maf = self._calculate_maf(sg, self.ploidy)
                if maf > snp_maf:
                    snp_id = sid
                    snp_maf = maf
                else:
                    continue
        return snp_id

    @staticmethod
    def worker(elimination_steps, snp_set, snp_genotypes, target_gen_n):
        tested_snp_set = snp_set.copy()
        for i in range(elimination_steps):
            snp_to_remove = random.choice(tested_snp_set)
            tested_snp_set.remove(snp_to_remove)  # remove random snp
            tested_geno = [
                snp_genotypes[i] for i in tested_snp_set
            ]  # extract genotypes of new snp set
            if np.unique(np.array(tested_geno), axis=1).shape[1] < target_gen_n:
                tested_snp_set.append(snp_to_remove)
            else:
                continue
        return tested_snp_set

    def optimal_snp_set_search(self):
        current_snp_set = []
        current_snps_geno = []

        # choose first snp
        current_snp = self.select_first_snp()
        current_snps_geno.append(self.snp_genotypes[current_snp])
        current_snp_set.append(current_snp)

        print("Primary SNP selection...")

        # identify primary set of SNPs
        try:
            while self._calc_N_distinct(current_snps_geno) < self.target_gen_n:
                parent_nodes_info = []
                for s, g in self.snp_genotypes.items():
                    if s in current_snp_set:
                        continue
                    else:
                        maf = self._calculate_maf(g, self.ploidy)
                        snp_info = (s, maf)
                        parent_nodes_info.append(snp_info)
                current_snp = sorted(
                    parent_nodes_info, key=lambda x: x[1], reverse=True
                )[0][0]
                current_snps_geno.append(self.snp_genotypes[current_snp])
                current_snp_set.append(current_snp)
        except IndexError:
            sys.exit("Not enough polymorphic SNP to discriminate samples. Exit.")

        print(f"After 1st step {len(current_snp_set)} primary SNP selected")

        print("Backward one-by-one elimination...")

        # one-by-one elimination
        with Pool(processes=self.ncups) as pool:
            res = pool.starmap(
                self.worker,
                [
                    (
                        self.elimination_steps,
                        current_snp_set,
                        self.snp_genotypes,
                        self.target_gen_n,
                    )
                    for _ in range(self.tries)
                ],
            )
        best_snp_set = min(res, key=len)
        print(f"{len(best_snp_set)} SNPs selected. Done")
        return best_snp_set

    def write_vcf(self, snp_list):
        with open(self.in_vcf) as invcf, open(self.out_vcf, "w") as outvcf:
            for l in invcf.readlines():
                if l.startswith("#"):
                    outvcf.write(l)
                else:
                    line = l.strip().split("\t")
                    if line[2] in snp_list:
                        outvcf.write(l)

    def main(self):
        selected_snps = self.optimal_snp_set_search()
        print("Writing selected SNPs in VCF...")
        self.write_vcf(selected_snps)
        print("Done")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("ivcf", help="input vcf file", type=parser_resolve_path)
    parser.add_argument("ovcf", help="output vcf file", type=parser_resolve_path)
    parser.add_argument(
        "-s",
        help="random seed (default: 810491)",
        default="810491",
        type=int,
        metavar="SEED",
    )
    parser.add_argument(
        "-e",
        help="number of backward one-by-one elimination steps (default: 10000)",
        type=int,
        default=10000,
        metavar="STEPS",
    )
    parser.add_argument(
        "-t",
        help="number of tries to find minimal SNP set (default: 1000)",
        default=1000,
        type=int,
        metavar="TRIES",
    )
    parser.add_argument(
        "-c", help="number of CPUs (default: 4)", default=4, type=int, metavar="CPU"
    )
    parser.add_argument(
        "-pl", help="VCF ploidy (default: 2)", default=2, type=int, metavar="PLOIDY"
    )

    args = parser.parse_args()

    MDSearch(
        in_vcf=args.ivcf,
        out_vcf=args.ovcf,
        elimination_steps=args.e,
        seed=args.s,
        ncups=args.c,
        tries=args.t,
        ploidy=args.pl,
    )
