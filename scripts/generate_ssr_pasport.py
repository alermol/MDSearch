from itertools import repeat, product
from pathlib import Path
import re
import numpy as np


def parser_resolve_path(path):
    return Path(path).resolve()


def main(table, outfile):
    print('Start working...')
    alleles = []
    with open(table) as file:
        for line in file:
            if not line.startswith('Name'):
                alleles.append(line.strip().split('\t')[1:])
    alleles = np.array(alleles)

    samples, markers = [], []
    with open(table) as file:
        for line in file:
            if not line.startswith('Name'):
                samples.append(line.strip().split('\t')[0])
            else:
                markers = line.strip().split('\t')[1:]
    markers = {str(m): [f'{i[0]}/{i[1]}' for i in zip(list(alleles[:,i[0]]), list(alleles[:,i[1]]))] 
               for m, i in zip(np.unique([i.split('_')[0] for i in markers]),
                               [(i, i + 1) for i in range(0, len(markers), 2)])}

    with open(outfile, 'w') as of:
        for i, s in enumerate(samples):
            alleles = '; '.join([f'{m}({g[i]})' for m, g in markers.items()])
            formula = (f'{s}: {alleles}').replace('**', '-')
            of.write(formula + '\n')
    print(f'Done. Output file: {outfile}')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("ifile", help="input file", type=parser_resolve_path)
    parser.add_argument("ofile", help="output file", type=parser_resolve_path)

    args = parser.parse_args()
    
    main(args.ifile, args.ofile)