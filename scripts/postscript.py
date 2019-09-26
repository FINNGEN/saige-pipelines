#!/usr/bin/env python3

import sys
import gzip

def run():
    '''
    Reads a GWAS summary stat file and writes a new one in PheWeb format to stdout filtering out variants
    1) if they are on the given blacklist or
    2) if chr is not numeric or X or
    3) if p-value == 0
    '''

    blacklist = { line.strip():"" for line in open(sys.argv[2], 'r').readlines() }

    filename = sys.argv[1]
    if filename.endswith('.gz') or filename.endswith('.bgz'):
        fh = gzip.open(filename, 'rt')
    else:
        fh = open(filename, 'rt')

    header = fh.readline().rstrip().split('\t')
    print('#chrom\tpos\tref\talt\tpval\tbeta\tsebeta\tmaf\tmaf_cases\tmaf_controls')
    for line in fh:
        split = line.rstrip().split('\t')
        try:
            chr = int(split[0].replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25'))
            if split[2] not in blacklist and float(split[11]) > 0:
                print('\t'.join([split[0].replace('chr', ''), split[1], split[3], split[4], split[11], split[8], split[9], split[6], split[18], split[19]]))
        except ValueError as e:
            pass
    fh.close()

if __name__ == '__main__':
    run()
