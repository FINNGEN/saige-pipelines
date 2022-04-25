#!/usr/bin/env python3

import argparse,os.path,shlex,subprocess,sys,subprocess,shlex,logging
from collections import defaultdict
from utils import basic_iterator,return_header,tmp_bash,file_exists,make_sure_path_exists,log_levels,extract_int_from_string,pretty_print

sub_dict =  {str(elem):str(elem) for elem in range(1,23)}
sub_dict.update({"X":"23"})
inv_sub_dict = {v: k for k, v in sub_dict.items()}

def filter_sumstats(regions_file,sumstats_file,subset_sumstats,chrom_list):

    print(f'filtering regions {subset_sumstats}')

    # print header to file so that even if regions are empty the header is present
    tmp_bash(f"tabix -H  {sumstats_file} > {subset_sumstats}")

    # new region to be used after chrom filtering
    region_tmp = subset_sumstats + ".region.bed"   
    with open(regions_file) as i,open(region_tmp,'wt') as o:
        for line in i:
            string = line.strip().split()[0]
            for key,val in sub_dict.items():
                string.replace(key,val)
            if extract_int_from_string(string) in chrom_list:
                o.write(line)                

    #regions can be empty if no hits!
    tot_regions = sum(1 for line in open(region_tmp))
    if tot_regions:
        logging.info(f"{tot_regions} regions to be filtered.")
        cmd = f"tabix {sumstats_file}  -R {region_tmp} >> {subset_sumstats}"
        tmp_bash(cmd)
    else:
        logging.warning(f"{regions_file} empty. No hits will be returned!")

    logging.info(f"{sum(1 for line in open(subset_sumstats)) -1} hits after filtering")
    return region_tmp

def read_regions(regions_file):
    """
    Return Chrom to regions dictionary. 
    """
    region_dict =defaultdict(list)
    with open(regions_file) as i:
        for line in i:
            chrom,start,end = line.strip().split()
            region_dict[chrom].append([int(start),int(end)])
    return region_dict


def filter_hits(subset_sumstats,regions_file,chrom_col = '#chrom',pos_col = 'pos',mlogp_col ='mlogp',ref_col = 'ref',alt_col='alt',pval_filter = 7):
    """
    Function that returns the top hit for each region. It reads in the region dict that gives for each chrom all the regions. 
    Then it loops over all variants, identifies the region and then checks the mlogp. If the hit is more signifcant, the top hit dict is updated.
    """

    region_dict = read_regions(regions_file)
    
    mlogs = defaultdict(lambda:pval_filter)
    hits = defaultdict(str)
 
    for elem in basic_iterator(subset_sumstats,skiprows=1,columns =  [chrom_col,pos_col,mlogp_col,ref_col,alt_col]):
        chrom,pos,mlogp,ref,alt = elem
        pos,mlogp = int(pos),float(mlogp)

        # find region in which variant is included
        regionid = ""
        for start,end in region_dict[chrom]:
            tmpid = f"{chrom}:{start}-{end}"
            if start <= pos <= end:
                regionid = tmpid
                break

        # region should always exist since tabix filters to the regions, but better to check!
        if regionid:
            if mlogp > mlogs[regionid]:
                mlogs[regionid] = mlogp
                hits[regionid] =  f"chr{inv_sub_dict[chrom]}_{pos}_{ref}_{alt}"
        else:
            logging.warning(f"{elem} does not match any region.")


    if not len(mlogs) == len(hits):
        logging.warning(f"{len(hits)} hits don't match number of regions ({len(mlogs)})")
    return hits,mlogs

def main(args):
    
    subset_sumstats = os.path.join(args.out,f"{args.pheno}_subset_regions.txt") 
    region_tmp = filter_sumstats(args.regions,args.sumstats,subset_sumstats,args.chroms)

    print("Processing data...",end="",flush=True)
    hits,mlogs = filter_hits(subset_sumstats,region_tmp,args.chr_col,args.pos_col,args.mlogp_col,args.ref_col,args.alt_col,args.pval_threshold)    
    print('done.')
    
    out_file = os.path.join(args.out,f"{args.pheno}_sig_hits.txt")
    chroms = list(set([region.split(':')[0] for region in hits]))
    chrom_files = [out_file.replace('.txt',f"_{chrom}.txt") for chrom in chroms]
    for f in chrom_files:
        if os.path.isfile(f):
            os.remove(f)
    
    print(f'dumping results to {out_file}')
    with open(out_file,'wt') as o:
        for region in hits:
            chrom,crange = region.split(':')
            chrom_file = out_file.replace('.txt',f"_{chrom}.txt")
            o.write('\t'.join(map(str,[args.pheno,chrom,region,hits[region],mlogs[region]])) + '\n')
            with open(chrom_file,'at') as tmp:
                tmp.write('\t'.join(map(str,[args.pheno,chrom,region,hits[region],mlogs[region]])) + '\n')
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Filtering of sumstats based on regions.")

    parser.add_argument('--pval_threshold','--pt',type = float,help ='Threshold limit (pvar -log(mpval) ',default = 7)
    parser.add_argument('--regions',type = file_exists,help ='Path to regions bed',required=True)
    parser.add_argument('--sumstats',type = file_exists,help ='Path to original sumstats',required=True)
    parser.add_argument('--chr_col', default="#chrom", type=str)
    parser.add_argument('--pos_col', default="pos", type=str)
    parser.add_argument('--ref_col', default="ref", type=str)
    parser.add_argument('--alt_col', default="alt", type=str)
    parser.add_argument('--mlogp_col', default="mlogp", type=str)
    parser.add_argument('--out',type = str,help ='Output Directory',required=True)
    parser.add_argument('--pheno',type = str,help ='Phenotype',required=True)
    parser.add_argument("--chroms",nargs='+',type=str,help="List of chromosomes to include",default =list(map(str,range(1,24))))
    parser.add_argument( "-log",  "--log",  default="warning", choices = log_levels, help=(  "Provide logging level. " "Example --log debug"))

    args = parser.parse_args()
    make_sure_path_exists(os.path.dirname(args.out))

    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")

    logging.info(f"{' '.join(args.chroms)} chromosomes included.")
    pretty_print(f"{args.pheno}")
    main(args)
    
   
