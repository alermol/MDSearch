import random
import numpy as np
from multiprocessing import Pool
from pathlib import Path
import re

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
        max_snps=None,
        min_dist=None,
        convert_het=None,
        n_sets=None
    ):
        random.seed(seed)
        self.in_vcf = in_vcf
        self.out_vcf = out_vcf
        self.elimination_steps = elimination_steps
        self.tries = tries
        self.ncups = ncups
        self.ploidy = ploidy
        self.max_snps = max_snps
        self.min_dist = min_dist
        self.convert_het = convert_het
        self.n_sets = n_sets

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
                    geno = []
                    for i in l.split("\t")[9:]:
                        if i.count('1') == 0:
                            geno.append(0)
                        elif i.count('1') == self.ploidy:
                            geno.append(1)
                        else:
                            if self.convert_het:
                                geno.append(np.nan)
                            else:
                                geno.append(round(1 / i.count('1'), 4))
                    self.snp_genotypes[snp_id] = [geno, l.split("\t")[9:]]
        self.main()

    @staticmethod
    def _calculate_maf(geno: list, ploidy: int):
        geno = "".join([str(i) for i in geno])
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
    def _calc_min_dist(snps: list):
        print(
            f'Calculate pairwise distance based on {len(snps)} SNPs...', end=' ')
        snps_array = np.array([i for i in snps])
        n_samples = snps_array.shape[1]
        distances = np.zeros((n_samples, n_samples), dtype=int)
        for i in range(n_samples):
            for j in range(n_samples):
                col_i = snps_array[:, i]
                col_j = snps_array[:, j]
                valid_mask = (~np.isnan(col_i) & ~np.isnan(col_j))
                distance = np.sum(
                    np.array(col_i[valid_mask]) != np.array(col_j[valid_mask]))
                distances[i, j] = distance
        res = min(distances[np.triu_indices(n_samples, k=1)])
        print(f'Minimal distance between samples: {res}')
        return res

    def select_first_snp(self):
        # select SNP with max MAF
        snp_id, snp_maf = None, None
        for sid, sg in self.snp_genotypes.items():
            if (snp_id is None) and (snp_maf is None):
                snp_id = sid
                snp_maf = self._calculate_maf(sg[1], self.ploidy)
            else:
                maf = self._calculate_maf(sg[1], self.ploidy)
                if maf > snp_maf:
                    snp_id = sid
                    snp_maf = maf
                else:
                    continue
        return snp_id

    def is_het(self, genotype: str):
        return len(set(re.findall(r'[0-9]+', genotype))) > 1

    @staticmethod
    def worker(elimination_steps, snp_set, snp_genotypes, min_dist, seed):
        def _calc_min_dist(snps: list):
            snps_array = np.array([i for i in snps])
            n_samples = snps_array.shape[1]
            distances = np.zeros((n_samples, n_samples), dtype=int)
            for i in range(n_samples):
                for j in range(n_samples):
                    col_i = snps_array[:, i]
                    col_j = snps_array[:, j]
                    valid_mask = (~np.isnan(col_i) & ~np.isnan(col_j))
                    distance = np.sum(
                        np.array(col_i[valid_mask]) != np.array(col_j[valid_mask]))
                    distances[i, j] = distance
            res = min(distances[np.triu_indices(n_samples, k=1)])
            return res

        random.seed(seed)
        tested_snp_set = snp_set.copy()
        for i in range(elimination_steps):
            snp_to_remove = random.choice(tested_snp_set)
            tested_snp_set.remove(snp_to_remove)  # remove random snp
            tested_geno = [
                snp_genotypes[i][0] for i in tested_snp_set
            ]  # extract genotypes of new snp set
            if _calc_min_dist(tested_geno) < min_dist:
                tested_snp_set.append(snp_to_remove)
            else:
                continue
        return tuple(tested_snp_set)

    def optimal_snp_set_search(self):
        current_snp_set = []
        current_snps_geno = []

        # choose first snp
        current_snp = self.select_first_snp()
        current_snps_geno.append(self.snp_genotypes[current_snp][0])
        current_snp_set.append(current_snp)

        print("Primary SNP selection...")

        # identify primary set of SNPs
        try:
            while self._calc_min_dist(current_snps_geno) < self.min_dist:
                print(
                    f'Current SNP set contains {len(current_snp_set)} SNPs...')
                parent_nodes_info = []
                for s, g in self.snp_genotypes.items():
                    if s in current_snp_set:
                        continue
                    else:
                        maf = self._calculate_maf(g[1], self.ploidy)
                        snp_info = (s, maf)
                        parent_nodes_info.append(snp_info)
                current_snp = sorted(
                    parent_nodes_info, key=lambda x: x[1], reverse=True
                )[0][0]
                current_snps_geno.append(self.snp_genotypes[current_snp][0])
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
                        self.min_dist,
                        random.uniform(10000, 10000000)
                    )
                    for _ in range(self.tries)
                ]
            )
        best_snp_sets = [list(i) for i in sorted(set([tuple(sorted(list(i))) for i in res]), key=len)[:self.n_sets]]
        print(f"{len(best_snp_sets)} discriminating SNP sets selected.")

        best_snp_sets_final = []
        for si, s in enumerate(best_snp_sets, start=1):
            orig_snp_number = len(s)
            if self.max_snps > len(s):
                snp_maf = {sid: self._calculate_maf(
                    g[1], self.ploidy) for sid, g in self.snp_genotypes.items() if sid not in s}
                snp_pic = sorted([(sid, 1 - ((maf ** 2) + ((1 - maf) ** 2)),) for sid, maf in snp_maf.items()],
                                 key=lambda x: x[1])[::-1]
                n_snps_to_add = self.max_snps - len(s)
                s += [i[0] for i in snp_pic[:n_snps_to_add]]
                print(
                    f"{n_snps_to_add} SNPs added to set {si} (orignial set contains {orig_snp_number} SNPs, total number of SNPs: {len(s)}).")
                best_snp_sets_final.append(s)
            else:
                best_snp_sets_final.append(s)
                print(
                    f"Addition of SNPs to discriminating set {si} is not required (total number of SNPs: {len(s)}).")
        return best_snp_sets_final

    def write_vcf(self, snp_sets_list):
        for si, s in enumerate(snp_sets_list, start=1):
            with open(self.in_vcf) as invcf, open(f'{self.out_vcf}_{si}.vcf', "w") as outvcf:
                for l in invcf.readlines():
                    if l.startswith("#"):
                        outvcf.write(l)
                    else:
                        line = l.strip().split("\t")
                        if (line[2] in s) & self.convert_het:
                            sep = '/' if '/' in line[9:][0] else '|'
                            line = line[:9] + [(f'.{sep}' * self.ploidy).rstrip(
                                sep) if self.is_het(i) else i for i in line[9:]]
                            line = '\t'.join(line) + '\n'
                            outvcf.write(line)
                        elif (line[2] in s) & (not self.convert_het):
                            outvcf.write(l)
                        else:
                            continue

    def main(self):
        selected_snps = self.optimal_snp_set_search()
        print("Writing selected SNPs in VCF...")
        self.write_vcf(selected_snps)
        print("Done")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("ivcf", help="input vcf file",
                        type=parser_resolve_path)
    parser.add_argument("ovcf_prefix", help="prefix of output vcf file",
                        type=parser_resolve_path)
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
    parser.add_argument(
        "-ts", help="Total number of SNPs in output set (Default: minimal discriminative set)", default=0, type=int, metavar="TOTAL SNP"
    )
    parser.add_argument(
        "-md", help="Minimal hamming distance between samples (Default: 1)", default=1, type=int, metavar="MIN DIST"
    )
    parser.add_argument(
        "-ch", help="Convert heterozygous calls into NA", action='store_true', default=False
    )
    parser.add_argument(
        "-ns", help="Number of distinct SNP sets in output (Default: 1)", default=1, type=int, metavar="N SETS"
    )

    args = parser.parse_args()

    MDSearch(
        in_vcf=args.ivcf,
        out_vcf=args.ovcf_prefix,
        elimination_steps=args.e,
        seed=args.s,
        ncups=args.c,
        tries=args.t,
        ploidy=args.pl,
        max_snps=args.ts,
        min_dist=args.md,
        convert_het=args.ch,
        n_sets=args.ns
    )
