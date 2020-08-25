#!/usr/bin/env python3

import sys
import gzip

def run():
    '''
    Reads a GWAS summary stat file and writes a new one in PheWeb format to stdout filtering out variants
    1) if they are on the given blacklist or
    2) if chr is not numeric or X/Y/MT or
    3) if p-value == 0 or
    4) if imputation info is less than a given value (if given)
    '''

    blacklist = { line.strip():"" for line in open(sys.argv[2], 'r').readlines() }
    try:
        info_threshold = float(sys.argv[3])
    except IndexError:
        info_threshold = -1

    filename = sys.argv[1]
    if filename.endswith('.gz') or filename.endswith('.bgz'):
        fh = gzip.open(filename, 'rt')
    else:
        fh = open(filename, 'rt')

    header = [h.lower() for h in fh.readline().rstrip().split('\t')]
    hd = {h:i for i, h in enumerate(header)}
    if 'af.cases' in hd:
        print('#chrom\tpos\tref\talt\tpval\tbeta\tsebeta\tmaf\tmaf_cases\tmaf_controls\tn_hom_cases\tn_het_cases\tn_hom_controls\tn_het_controls')
    else:
        print('#chrom\tpos\tref\talt\tpval\tbeta\tsebeta\tmaf\tn_hom_cases\tn_het_cases\tn_hom_controls\tn_het_controls')
    for line in fh:
        split = line.rstrip().split('\t')
        try:
            chr = int(split[hd['chr']].replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25'))
            if split[hd['snpid']] not in blacklist and float(split[hd['p.value']]) > 0 and float(split[hd['imputationinfo']]) >= info_threshold:
                if 'af.cases' in hd:
                    print('\t'.join([
                        split[hd['chr']].replace('chr', ''),
                        split[hd['pos']],
                        split[hd['allele1']],
                        split[hd['allele2']],
                        split[hd['p.value']],
                        split[hd['beta']],
                        split[hd['se']],
                        split[hd['af_allele2']],
                        split[hd['af.cases']],
                        split[hd['af.controls']],
                        split[hd['homn_allele2_cases']],
                        split[hd['hetn_allele2_cases']],
                        split[hd['homn_allele2_ctrls']],
                        split[hd['hetn_allele2_ctrls']]]))
                else:
                    print('\t'.join([
                        split[hd['chr']].replace('chr', ''),
                        split[hd['pos']],
                        split[hd['allele1']],
                        split[hd['allele2']],
                        split[hd['p.value']],
                        split[hd['beta']],
                        split[hd['se']],
                        split[hd['af_allele2']],
                        split[hd['homn_allele2_cases']],
                        split[hd['hetn_allele2_cases']],
                        split[hd['homn_allele2_ctrls']],
                        split[hd['hetn_allele2_ctrls']]]))
        except ValueError as e:
            # invalid chromosome
            pass
    fh.close()

if __name__ == '__main__':
    run()
