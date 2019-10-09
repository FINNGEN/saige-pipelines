#!/usr/bin/env python3.7
import argparse
import pandas as pd
import re
from functools import partial
import csv
from pandas.api.types import is_string_dtype
PHENO_TAG="{PHENO}"
CHR_TAG="{CHR}"

def row_conf(pheno, chr, start, stop, condition_snp, null_mask, var_ratio_mask, bgen_mask, sample_file, p_thr):

    return null_mask.replace(PHENO_TAG, pheno).cat( [var_ratio_mask.replace(PHENO_TAG,pheno) ,
                                                      bgen_mask.replace(PHENO_TAG, pheno),
                                                      bgen_mask.replace(PHENO_TAG, pheno),
                                                      sample_file,
                                                      chr, start,stop, condition_snp, p_thr], sep="\t")


def row_conf_line(row, null_mask, var_ratio_mask, bgen_mask, sample_file, p_thr, add_chr=False):
    chrom = f'chr{row.chr}' if not row.chr.startswith('chr') else row.chr
    chr_n = chrom.replace("chr","")
    snp = f'chr{row.condition_snp}' if not row.condition_snp.startswith('chr') else row.condition_snp
    return f'{re.sub(PHENO_TAG,row.pheno,null_mask,flags=re.IGNORECASE)}\t{re.sub(PHENO_TAG,row.pheno,var_ratio_mask,flags=re.IGNORECASE)}\t' \
        f'{re.sub(CHR_TAG,chr_n,bgen_mask,flags=re.IGNORECASE)}\t{re.sub(CHR_TAG,chr_n,bgen_mask,flags=re.IGNORECASE)}.bgi\t' \
        f'{sample_file}\t{chrom}\t{row.start}\t{row.stop}\t{snp}\t{p_thr}'



def from_locus(args):
    loci = pd.read_csv(args.locus_file, sep="\t")
    loci.columns = map(str.lower, loci.columns)
    r = loci.apply(partial(row_conf_line, null_mask=args.null_file_mask, var_ratio_mask=args.var_ratio_mask,
                           bgen_mask=args.bgen_file_mask, sample_file=args.sample_file, p_thr=args.p_threshold), 1)
    r.to_csv(args.output, index=False, header=True)



def from_variants(args):

    vars = pd.read_csv(args.snps_file, sep="\t", compression='gzip' if args.snps_file.endswith('gz') else 'infer',
                       dtype={args.p_col:float, args.chr_col:str})


    vars = vars[ vars[args.p_col]<args.p_threshold ].sort_values(by=[args.pheno_col,args.p_col])
    if args.exclude_regions:
        excl = args.exclude_regions.split(";")

        for e in excl:
            er = e.split(",")
            print(f'Excluding {er}')
            if len(er)!=3:
                raise Exception("Illegal exclude region definition. Define regions by ; separated list of chr,start, stop")
            vars = vars[ (vars[args.chr_col]!=er[0]) |  (vars[args.pos_col] < int(float(er[1]))) |  (vars[args.pos_col] > int(float(er[2]))) ]

    loc_width = args.locus_padding

    loci = []

    while len(vars.index) > 0:
        top = vars.iloc[0]
        loc_start = max(top[args.pos_col] - loc_width,1)
        loc_end = top[args.pos_col] + loc_width
        loci.append( top )
        ## we can sort by pheno and position and do the scanning linearly (n + n*log(n) )
        vars = vars[ (vars[args.pheno_col] != top[args.pheno_col]) | (vars[args.chr_col]!= top[args.chr_col]) | (vars[args.pos_col]<loc_start) | ( vars[args.pos_col]>loc_end ) ]

    loci_df = pd.concat(loci, axis=1).T.sort_values(by=[args.pheno_col, args.chr_col, args.pos_col] )

    merged = []

    def get_cols(s):
        return pd.Series({"pheno": s.PHENO, "chr": s[args.chr_col],
                   "start": max(s[args.pos_col] - loc_width,1), "stop": s[args.pos_col] + loc_width,
                   "condition_snp":f'{s[args.chr_col]}_{s[args.pos_col]}_{s[args.ref_col]}_{s[args.alt_col]}', "pval":s[args.p_col] })

    for i in range(0, len(loci_df.index) ):
        curr = get_cols(loci_df.iloc[i])
        if (len(merged) > 0 and merged[-1].pheno == curr.pheno and
                curr.chr == merged[-1].chr and curr.start < merged[-1].stop):
            win = curr if curr.pval < merged[-1].pval else merged[-1]
            win.start = merged[-1].start
            win.stop = curr.stop
            merged[-1] = win
        else:
            merged.append(curr)

    mergs = pd.concat(merged, axis=1).T


    r = mergs.apply(partial(row_conf_line, null_mask=args.null_file_mask, var_ratio_mask=args.var_ratio_mask,
                           bgen_mask=args.bgen_file_mask, sample_file=args.sample_file, p_thr=args.p_condition_threshold, add_chr=args.add_chr), 1)

    r.to_csv(args.output + ".merged", quoting=csv.QUOTE_NONE, index=False, header=False)
    loci_df.to_csv(args.output,sep="\t", index=False, header=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run conditional saige until no significant results remaining")
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

    parser.add_argument('--p_condition_threshold', action='store', type=float, default=5e-8,
                        help='p-value threshold for signal inclusion after first conditional analysis')

    cmd_parsers = parser.add_subparsers(title="commands", help="Sub command help")
    loci_cmd = cmd_parsers.add_parser('loci')
    loci_cmd.add_argument('locus_file', action='store', type=str,
                        help='file with 4 columns: pheno chr start stop condition_snp')
    loci_cmd.set_defaults(func=from_locus)

    vars_cmd = cmd_parsers.add_parser('vars')
    vars_cmd.add_argument('snps_file', action='store', type=str,
                          help='file giving candidate variants and phenotypes (e.g. all gw-sig variants) with at least 4 columns: PHENO chr pos ref alt p-value')
    vars_cmd.add_argument('--locus_padding', action='store', type=int, default=1500000,
                          help='Nucleotides +- from lead snp')

    vars_cmd.add_argument('--exclude_regions', action='store', type=str,
                                                help='regions to exclude semicolon separated list of chr,start,end triplets')

    vars_cmd.set_defaults(func=from_variants)

    vars_cmd.add_argument('--pheno_col', default="PHENO", type=str)
    vars_cmd.add_argument('--chr_col', default="#chrom", type=str)
    vars_cmd.add_argument('--pos_col', default="pos", type=str)
    vars_cmd.add_argument('--ref_col', default="ref", type=str)
    vars_cmd.add_argument('--alt_col', default="alt", type=str)
    vars_cmd.add_argument('--p_col', default="pval", type=str)

    vars_cmd.add_argument('--add_chr', action='store_true', help="Adds chr prefix to chromosomes if they do not exist. Useful if bgens have chr but summary stats don't.")

    args = parser.parse_args()

    if not "func" in vars(args):
        raise Exception("Command must be given")

    args.func(args)
