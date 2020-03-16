#!/usr/bin/env python3.7
import argparse
import pandas as pd
import os
import subprocess
import logging


def run_cond(null, var_ratio, bgen,bgen_idx,sample, chr, start, stop, output, condition_on, minmac, script):

    with open("tmp.region",'w') as region:
        region.write(f'{chr} {start} {stop}')
    logging.info(f'Running conditioning on { ",".join(condition_on) } for locus {chr}:{start}-{stop}')
    cmd = f'{script} --bgenFile {bgen} --bgenFileIndex {bgen_idx} --sampleFile={sample} --minMAC {minmac}' \
            f' --GMMATmodelFile {null} --varianceRatioFile {var_ratio} --SAIGEOutputFile {output}'  \
            f' --condition {",".join(condition_on)} --rangestoIncludeFile tmp.region'
    logging.info(f'Running command {cmd}')
    try:
        subprocess.check_call(cmd, shell=True)
    finally:
        os.system("rm tmp.region")


def plot_locus(result_file):

    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run conditional saige until no significant results remaining")
    parser.add_argument('null_file', action='store', type=str, help='SAIGE null file')
    parser.add_argument('variance_ratio_file', action='store', type=str, help='Result file')
    parser.add_argument('bgen_file', action='store', type=str, help='Result file')
    parser.add_argument('bgen_index', action='store', type=str, help='Result file')
    parser.add_argument('sample_file', action='store', type=str, help='Result file')
    parser.add_argument('chr', action='store', type=str, help='Result file')
    parser.add_argument('start', action='store', type=int, help='Result file')
    parser.add_argument('stop', action='store', type=int, help='Result file')
    parser.add_argument('condition_snp', action='store', type=str, help='Result file')

    parser.add_argument('--out_prefix', action='store', default="", type=str, help='')
    parser.add_argument('--max_rounds', action='store', default=10, type=int, help='Maximum number of conditioning rounds to perform')
    parser.add_argument('--p_threshold', action='store', default=5e-8, type=float, help='')
    parser.add_argument('--min_mac', action="store", default=10, type=int)
    parser.add_argument('--saige_script', action="store", default="step2_SPAtests.R")
    parser.add_argument('--plot_locus', action="store_true", default=False )

    args = parser.parse_args()

    logging.basicConfig(filename=f'{args.out_prefix+"_" if len(args.out_prefix)>0 else ""}{os.path.basename(args.null_file)}.conditional.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.INFO)

    out = f'{args.out_prefix}{"_" if len(args.out_prefix)>0 else ""}{os.path.basename(args.null_file)}_{args.condition_snp}_1.conditional'
    cond_snps = [args.condition_snp]
    run_cond(args.null_file, args.variance_ratio_file, args.bgen_file, args.bgen_index, args.sample_file,args.chr,
             args.start, args.stop, out, cond_snps, args.min_mac, args.saige_script)

    if args.plot_locus:
        plot_locus(out)
    columns = ["SNPID","BETA","SE","p.value","BETA_cond","SE_cond","p.value_cond"]
    res = pd.read_csv(out, sep=" ")
    r = res.loc[ (~res['SNPID'].isin(cond_snps)) & (res["p.value_cond"] < args.p_threshold),
                    columns].sort_values(by=["p.value_cond"])

    summary_out = f'{args.out_prefix+"_" if len(args.out_prefix)>0 else ""}{os.path.basename(args.null_file)}_{args.condition_snp}.independent.snps'

    tab="\t"
    with open( summary_out, 'w' ) as summary:
        summary.write("\t".join(columns) + "\tConditioned_on" + "\n"  )
        cond_snp=res.loc[ res['SNPID']== args.condition_snp, columns]
        summary.write(f'{tab.join( str(v) for v in cond_snp.iloc[0].values)}\t{",".join(cond_snps)}\n')
        while not r.empty and len(cond_snps)+1<=args.max_rounds:
            summary.write(f'{tab.join( str(v) for v in r.iloc[0].values)}\t{",".join(cond_snps)}\n')
            cond_snps.append(r.iloc[0].SNPID)
            out = f'{args.out_prefix + "_" if len(args.out_prefix) > 0 else ""}{os.path.basename(args.null_file)}_{args.condition_snp}_{len(cond_snps)}.conditional'

            run_cond(args.null_file, args.variance_ratio_file, args.bgen_file, args.bgen_index, args.sample_file, args.chr,
                     args.start, args.stop, out, cond_snps, args.min_mac, args.saige_script )
            if args.plot_locus:
                plot_locus(out)

            res = pd.read_csv(out, sep=" ")
            r = res.sort_values(by=["p.value_cond"]).head().loc[
                (~res['SNPID'].isin(cond_snps)) & (res["p.value_cond"] < args.p_threshold), columns]
