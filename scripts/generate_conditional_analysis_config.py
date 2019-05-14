#!/usr/bin/env python3
import argparse
import pandas as pd
import re
from functools import partial

PHENO_TAG="{PHENO}"

def row_conf(pheno, chr, start, stop, condition_snp, null_mask, var_ratio_mask, bgen_mask, sample_file, p_thr):

    return null_mask.replace(PHENO_TAG, pheno).cat( [var_ratio_mask.replace(PHENO_TAG,pheno) ,
                                                      bgen_mask.replace(PHENO_TAG, pheno) + "bgen",
                                                      bgen_mask.replace(PHENO_TAG, pheno) + "bgi",
                                                      sample_file,
                                                      chr, start,stop, condition_snp, p_thr], sep="\t")


def row_conf_line(row, null_mask, var_ratio_mask, bgen_mask, sample_file, p_thr):

    return f'{re.sub(PHENO_TAG,row.pheno,null_mask,flags=re.IGNORECASE)}\t{re.sub(PHENO_TAG,row.pheno,var_ratio_mask,flags=re.IGNORECASE)}\t' \
        f'{re.sub(PHENO_TAG,row.pheno,bgen_mask,flags=re.IGNORECASE)}.bgen\t{re.sub(PHENO_TAG,row.pheno,bgen_mask,flags=re.IGNORECASE)}.bgi\t' \
        f'{sample_file}\t{row.chr}\t{row.start}\t{row.stop}\t{row.condition_snp}\t{p_thr}'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run conditional saige until no significant results remaining")
    parser.add_argument('locus_file', action='store', type=str, help='file with 4 columns: pheno chr start stop condition_snp')
    parser.add_argument('null_file_mask', action='store', type=str, help='mask that maps to file location of null run file. '
                                                                         ' use {PHENO} as replacement for the pheno name( e.g. gs://data/R3.grm.v1-${PHENO}.rda)')
    parser.add_argument('var_ratio_mask', action='store', type=str,
                        help='mask that maps to file location of null run variance ratio file. '
                             ' use ${PHENO} as replacement for the pheno name( e.g. gs://data/R3.grm.v1-${PHENO}.varianceRatio.txt)')
    parser.add_argument('bgen_file_mask', action='store', type=str,
                        help='mask that maps to file location of bgen file without postfix.  ')
    parser.add_argument('sample_file', action='store', type=str,
                        help='mask that maps to file location of null run files. '
                             ' use ${PHENO} as replacement for the pheno name( e.g. gs://data/R3.grm.v1-${PHENO}.rda)')

    parser.add_argument('output', action='store', type=str,
                        help='output file to be used as config to saige_conditional.wdl')

    parser.add_argument('--p_threshold', action='store', type=float, default=5e-8,
                        help='p-value threshold for signal inclusion')

    args = parser.parse_args()

    loci = pd.read_csv(args.locus_file, sep="\t")
    loci.columns = map(str.lower, loci.columns)
    r = loci.apply( partial(row_conf_line, null_mask=args.null_file_mask, var_ratio_mask=args.var_ratio_mask,
                        bgen_mask=args.bgen_file_mask, sample_file=args.sample_file, p_thr=args.p_threshold),1)
    r.to_csv(args.output, index=False)

