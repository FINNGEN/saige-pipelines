#!/usr/bin/env python3

import argparse,os.path,shlex,subprocess,sys,subprocess,shlex,json,logging,multiprocessing,time
import numpy as np
from utils import file_exists,make_sure_path_exists,tmp_bash,pretty_print,return_open_func,log_levels,basic_iterator,return_header,check_region
from pathlib import Path
import pandas as pd
from collections import defaultdict as dd

regenie_covariates= "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm,BATCH_DS18_MIGRAINE_1_norm,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm,BATCH_DS21_SUPER_2_norm,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm"

sub_dict =  {str(elem):str(elem) for elem in range(1,23)}
sub_dict.update({"X":"23"})
inv_sub_dict = {v: k for k, v in sub_dict.items()}

def filter_pheno(args):
    """
    Extracts pheno column from big file so i don't have to reextract it every time.
    """
    header = return_header(args.pheno_file)
    columns = [1,2] + [1+header.index(elem) for elem in [args.pheno] + args.covariates.split(',')]
    tmp_pheno =   os.path.join(args.tmp_dir,f"{args.pheno}_pheno.tmp")
    if not os.path.isfile(tmp_pheno) or args.force:
        print("Creating new pheno file...",end="",flush=True)
        cmd = f"zcat -f {args.pheno_file}  | cut -f {','.join(map(str,columns))} > {tmp_pheno}"
        tmp_bash(cmd)
        logging.debug(cmd)
        print('done.')
    return tmp_pheno

def regenie_run(args,step,bgen,sample_file,pheno_file,covariates,condition_list,locus,null_file,pheno,region,log_file,threads,regenie_cmd = "regenie",params = ' ',):
    """
    Single regenie conditainal run. 
    Returns file with results after doing some renaming.
    """
   
    regenie_file = os.path.join(args.tmp_dir, f"{args.basename}_{pheno}.regenie") #default regenie output
    out_file = os.path.join(args.out_dir, f"{args.basename}_{pheno}_{locus}_{step}.conditional") #file where to redirect output
    pretty_print(f"VARIANT:{condition_list[-1]}")
    logging.info(f"generating {out_file} ...")

    if not os.path.isfile(out_file) or args.force:
        args.force = True

        #build pred
        pred_file = os.path.join(args.tmp_dir, f"{args.basename}_{args.pheno}.pred")                            
        with open(pred_file,'wt') as o: o.write(f"{pheno}\t{null_file}")
   
        # build condition file
        tmp_variant = os.path.join(args.tmp_dir, f"{args.basename}_variant.tmp")
        with open(tmp_variant,'wt') as o:
            for variant in condition_list:
                o.write(variant + '\n')

        # add sample file if passed
        sample_cmd = f" --sample {sample_file}"  if os.path.isfile(sample_file)  else ""
        
        cmd = f'{regenie_cmd} --step 2   {params} --bgen {bgen}  {sample_cmd} --out {os.path.join(args.tmp_dir,args.basename)}  --pred {pred_file} --phenoFile {pheno_file} --phenoCol {pheno} --condition-list {tmp_variant} {region}  --covarFile {pheno_file} --covarColList {covariates} --threads {threads}'
        logging.debug(cmd)
        logmode = 'wt' if step == 0 else 'a'
        with open(log_file,logmode) as o:
            start = time.time()
            ret = subprocess.call(shlex.split(cmd),stdout=o)
            logging.info(f"Script ran in {time.time() - start} seconds with {threads} cpus.")
            
        #spring cleaning if run is successful
        if not ret:
            os.replace(regenie_file,out_file)

    else:
        logging.info('file already exists')
        ret = 0

        
    return out_file,ret

def parse_sumstat_data(sumstats,pval_dict_file,column_names = ['#chrom','pos','mlogp','ref','alt','beta','sebeta']):
    """
    Reads in sumstats as dictionary for relevant fields.
    """
    sum_dict = dd(lambda : dd(float))

    print(pval_dict_file)
    if not os.path.isfile(pval_dict_file):
        logging.info("reading original pvals..")
        
        header =return_header(sumstats)
        columns = [header.index(elem) for elem in column_names]
        it = basic_iterator(sumstats,skiprows=1,columns = columns)
        for elem in it:
            chrom,pos,mlogp,ref,alt,beta,sebeta = elem
            variant ="chr" + '_'.join([inv_sub_dict[chrom],pos,ref,alt])
            for key in ["beta","sebeta","mlogp"]:
                sum_dict[variant][key] = eval(key)
            
        logging.info(f"dumping pvals..{pval_dict_file} {len(sum_dict)}")
        with open(pval_dict_file, 'wt') as fp:
            json.dump(sum_dict, fp)
        logging.info('done.')
    else:
        logging.info('reading in json...')
        with open(pval_dict_file) as json_file:
            sum_dict = json.load(json_file)
        logging.info('done.')
    logging.info(f"{len(sum_dict)} original variants in sumstats")
    return sum_dict

def check_hit(out_file,step,threshold=7):
    """
    Checks if regenie hit is significant. 
    """

    # read results as pandas df and get row with max -log10(p)
    df = pd.read_csv(out_file,sep=" ")
    # get max value,extract relevant fields and map to str
    data = list(map(str,df.iloc[df['LOG10P'].idxmax()].loc[['ID','BETA','SE','LOG10P']].values))
    variant,beta,se,pval = data
    print(f"CANDIDATE VARIANT:{variant}")
    print(f"Variant info from conditioned sumstats {map_vals_to_string(data[1:])}")
    step = step+ 1 if  float(pval) > threshold else 0
    return step,data


def map_vals_to_string(vals):
    return  str(dict(zip(['beta','sebeta','mlogp'],[round(float(elem),2) for elem in vals])))


def get_sum_dict_data(sum_dict,variant):
    keys =  ["beta","sebeta","mlogp"]
    return [sum_dict[variant][elem] for elem in keys]

def main(locus,region,args,tmp_pheno_file,sum_dict,threads):

    pretty_print(f"{locus} {region} conditional chain.",50)
    #define log_file
    log_file = args.out + f"_{args.pheno}_{locus}.log"
    print(f"Logging to {log_file}")
    
    result_file = args.out + f"_{args.pheno}_{locus}.independent.snps"
    #inital values to start looping
    condition_variant = locus
    with open(result_file,'wt') as o:
        header = ['VARIANT','BETA','SE','MLOG10P','BETA_cond','SE_cond','MLOG10P_cond','VARIANT_cond']
        o.write("\t".join(header) + '\n')
        o.write('\t'.join( [condition_variant]+ get_sum_dict_data(sum_dict,condition_variant)  +["NA"]*4    ) + '\n')

        print(f"Variant info from original FG sumstats {map_vals_to_string(get_sum_dict_data(sum_dict,condition_variant))}")
        step,condition_list = 1,[locus]
        # step is non 0 if hit is significant. Else truncate when max step is reached 
        while 0 <  step <= args.max_steps:
            out_file,ret = regenie_run(args,step,args.bgen,args.sample_file,tmp_pheno_file,args.covariates,condition_list,locus,args.null_file,args.pheno,region,log_file,threads,params = args.regenie_params)
            if ret:
                logging.error(f"RUN FAILED,check {log_file} for errors.")
                return
            # now that we have results, let's analyze them to make sure the top hit is significant
            step,data =check_hit(out_file,step,args.pval_threshold)
            if step:
                new_variant = data[0]
                o.write('\t'.join( [new_variant]+ get_sum_dict_data(sum_dict,new_variant)  + data[1:] +   [','.join(condition_list)]) + '\n')
                #update variant to conidtion on
                condition_variant = new_variant
                # update condition list
                condition_list.append(new_variant)
                print(f"Variant info from original FG sumstats {map_vals_to_string(get_sum_dict_data(sum_dict,condition_variant))}")
                print(f"Hit signficant, proceeding to condition...")
            else:
                print("Hit not significant. Ending loop.")



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Recursive conditional analysis for regenie.")

    parser.add_argument('--pval-threshold',type = float,help ='Threshold limit (-log(mpval)) ',default = 7)
    parser.add_argument('--pheno',type = str,help ='Pheno column in pheno file',required=True)
    parser.add_argument('--out',type = str,help ='Output Directory and prefix (e.g. /foo/bar/finngen)',required=True)
    parser.add_argument('--covariates',type = str,default = regenie_covariates,help='List of covariates')
    parser.add_argument('--pheno-file',type = file_exists,help ='Path to pheno file',required=True)
    parser.add_argument('--bgen',type = file_exists,help ='Path to bgen',required=True)
    parser.add_argument('--sample-file',type = file_exists,help ='Path to bgen sample file (if not in same directory as bgen)',required=False)
    parser.add_argument('--sumstats',type = file_exists,help ='Path to original sumstats',required=True)
    parser.add_argument('--regenie-params',type=str,help ='extra bgen params',default = ' --bt --bsize 200 ' )
    parser.add_argument('--null-file',type = file_exists,help ='File with null info.',required=True)
    parser.add_argument('--force',action = 'store_true',help = 'Flag for forcing re-run.')
    parser.add_argument( "-log",  "--log",  default="warning", choices = log_levels, help=(  "Provide logging level. " "Example --log debug', default='warning'"))
    parser.add_argument('--max-steps',type = int,default =10)
    parser.add_argument('--chr_col','--chr-col', default="#chrom", type=str)
    parser.add_argument('--pos_col','--pos-col', default="pos", type=str)
    parser.add_argument('--ref_col','--ref-col', default="ref", type=str)
    parser.add_argument('--alt_col','--alt-col',default="alt", type=str)
    parser.add_argument('--mlogp_col','--mlogp-col', default="mlogp", type=str)
    parser.add_argument('--beta_col','--beta-col', default="beta", type=str)
    parser.add_argument('--sebeta_col','--sebeta-col', default="sebeta", type=str)
    parser.add_argument('--threads',type = int,help ='Number of threads.',default =  multiprocessing.cpu_count())


    range_group = parser.add_mutually_exclusive_group(required=True)
    range_group.add_argument('--locus-region',type =str,nargs=2,help ='Locus & Region to filter CHR:START-END')
    range_group.add_argument('--locus-list',type = file_exists,help="File with list of locus and regions")

    args = parser.parse_args()
    

    args.out_dir,args.basename = os.path.dirname(args.out),os.path.basename(args.out)
    args.tmp_dir = os.path.join(args.out_dir,'tmp')
    make_sure_path_exists(args.out_dir)
    make_sure_path_exists(args.tmp_dir)

    # logging level
    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")
    
    # automatically look for sample file if not explicitly passed.
    if not args.sample_file:
        for sample_file in [args.bgen.replace('.bgen',elem) for elem in ['.bgen.sample','.sample']]:
            if os.path.isfile(sample_file):
                args.sample_file = sample_file
                logging.warning(f"Using {args.sample_file} as sample file.")

    # figure out whether threshold is pval or mlogp
    if args.pval_threshold <1:args.pval_threshold = -np.log10(args.pval_threshold)
    pretty_print(f"MLOGP THRESHOLD: {args.pval_threshold}")
    
    # gets original sumstats data for variants
    pval_dict_file = os.path.join(args.tmp_dir,f"{args.pheno}_pvals.json")
    columns= [args.chr_col,args.pos_col,args.mlogp_col,args.ref_col,args.alt_col,args.beta_col,args.sebeta_col]
    sum_dict = parse_sumstat_data(args.sumstats,pval_dict_file,columns)

    # formatting of args for regenie
    if args.locus_region:  region_list  = [check_region(*args.locus_region)]
    if args.locus_list:
        region_list = []
        with open(args.locus_list) as i:
            for line in i:
                region_list.append(check_region(*line.strip().split()))

    logging.info(region_list)
    logging.info(f"Using {args.threads} cpus for the regenie run.")

    # create tmp pheno file (lighter)
    tmp_pheno_file = filter_pheno(args) 
    for locus,region in region_list:
        main(locus,region,args,tmp_pheno_file,sum_dict,args.threads)
                
    
    
