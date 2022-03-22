#!/usr/bin/env python3

import argparse,os.path,shlex,subprocess,sys,subprocess,shlex
from collections import defaultdict
from utils import basic_iterator,return_header,tmp_bash,file_exists,make_sure_path_exists



def extract_regions(regions_file,sumstats_file):

    out_file = sumstats_file + '.filtered_regions'
    cmd = f"tabix {sumstats_file} -h -R {regions_file} > {out_file}"
    tmp_bash(cmd)
    
    
def read_regions(regions_file):

    region_dict =defaultdict(list)

    with open(regions_file) as i:
        for line in i:
            chrom,start,end = line.strip().split()
            region_dict[chrom].append([start,end])

    return region_dict


def filter_hits(sumstats_file,regions_file,chrom_col = '#chrom',pos_col = 'pos',mlogp_col ='mlogp',ref_col = 'ref',alt_col='alt',pval_filter = 7):


    region_dict = read_regions(regions_file)
    print("region dict imported")
    subset_sumstats = sumstats_file + '.filtered_regions'
    if not os.path.isfile(subset_sumstats):
        print(f'filtering regions {subset_sumstats} ....',end = "",flush=True)
        extract_regions(regions_file,sumstats_file)
        print('done.')
    else:
        print('regions already filtered')
     
    header =return_header(subset_sumstats)
    columns = [header.index(elem) for elem in [chrom_col,pos_col,mlogp_col,ref_col,alt_col]]
    it = basic_iterator(subset_sumstats,skiprows=1,columns = columns)

    mlogs = defaultdict(lambda:pval_filter)
    hits = defaultdict(str)
    
    for elem in it:
        chrom,pos,mlogp,ref,alt = elem

        for region in region_dict[chrom]:
            start,end = region
            regionid = f"{chrom}:{start}-{end}"
            if start <= pos <= end:
                snpid = f"chr{chrom}_{pos}_{ref}_{alt}"
                if float(mlogp) > mlogs[regionid]:
                    mlogs[regionid] = float(mlogp)
                    hits[regionid] = snpid


    return hits,mlogs
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Recursive conditional analysis for regenie.")

    parser.add_argument('--pval_threshold',type = float,help ='Threshold limit (pvar -log(mpval) ',default = 7)
    parser.add_argument('--regions',type = file_exists,help ='Path to regions bed',required=True)
    parser.add_argument('--sumstats',type = file_exists,help ='Path to original sumstats',required=True)
    parser.add_argument('--chr_col', default="#chrom", type=str)
    parser.add_argument('--pos_col', default="pos", type=str)
    parser.add_argument('--ref_col', default="ref", type=str)
    parser.add_argument('--alt_col', default="alt", type=str)
    parser.add_argument('--mlogp_col', default="mlogp", type=str)
    parser.add_argument('--out',type = str,help ='Output Directory',required=True)
    parser.add_argument('--pheno',type = str,help ='Pheno column',required=True)

    args = parser.parse_args()
    make_sure_path_exists(os.path.dirname(args.out))

    hits,mlogs = filter_hits(args.sumstats,args.regions,args.chr_col,args.pos_col,args.mlogp_col,args.ref_col,args.alt_col,args.pval_threshold)
    
    out_file = os.path.join(args.out,f"{args.pheno}_sig_hits.txt")
    print(f'dumping results to {out_file}')
    with open(out_file,'wt') as o:
        for hit in hits:
            o.write('\t'.join(map(str,[args.pheno,hit,hits[hit],mlogs[hit]])) + '\n')
