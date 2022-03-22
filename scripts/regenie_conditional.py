
#!/usr/bin/env python3

import argparse,os.path,shlex,subprocess,sys,subprocess,shlex,json,logging
import numpy as np
from utils import file_exists,make_sure_path_exists,tmp_bash,pretty_print,return_open_func,log_levels,basic_iterator,return_header
import pandas as pd

regenie_covariates= "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC{1:10},IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm,BATCH_DS18_MIGRAINE_1_norm,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm,BATCH_DS21_SUPER_2_norm,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm"


def regenie_run(out_root,step,bgen,pheno_file,covariates,condition_variant,locus,null_file,pheno,region,exclusion_list,regenie_cmd = "regenie",params = '--bt --bsize 200 ',):
    """
    Single regenie conditainal run. 
    Returns file with results after doing some renaming.
    """

    regenie_file = out_root + f"_{pheno}.regenie" #default regenie output
    out_file = out_root + f"_{pheno}_{locus}_{step}.conditional" #file where to redirect output
    pretty_print(f"VARIANT:{condition_variant}")
    logging.info(f"generating {out_file} ...")
 
    if not os.path.isfile(out_file) or args.force:
        args.force = True

        #build pred
        pred_file = out_root + f"_{args.pheno}.pred"
        with open(pred_file,'wt') as o: o.write(f"{pheno}\t{null_file}")

        # build exclusion list
        tmp_exclusion = out_root + '.exclude.tmp'
        logging.info(f"Variants to exclude :{' '.join(exclusion_list)}")
        with open(tmp_exclusion,'wt') as o:
            o.write('\n'.join(exclusion_list))
    
        # build condition file
        tmp_variant = out_root + '.variant.tmp'
        logging.info(f"Variant to condition on :{condition_variant}")
        with open(tmp_variant,'wt') as o: o.write(condition_variant)
        
        cmd = f'{regenie_cmd} --step 2 {params} --bgen {bgen}  --out {out_root} --exclude {tmp_exclusion} --pred {pred_file} --phenoFile {pheno_file} --phenoCol {pheno} --condition-list {tmp_variant} {region}  --covarFile {pheno_file} --covarColList {covariates}'
        logging.debug(cmd)
        log_mode = 'wt' if step ==1 else 'a'
        with open(out_root+ f"_{args.pheno}_{locus}.log",log_mode) as o:
            subprocess.call(shlex.split(cmd),stdout=o)

        #spring cleaning
        for f in [pred_file,tmp_variant,out_root + '.log'] : os.remove(f)
        os.replace(regenie_file,out_file)

    else:
        logging.info('file already exists')     

    return out_file

def parse_sumstat_data(sumstats,pval_dict_file,column_names = ['#chrom','pos','mlogp','ref','alt','beta','sebeta']):
    """
    Reads in sumstats as dictionary for relevant fields.
    """
    sum_dict = {}

    if not os.path.isfile(pval_dict_file):
        logging.warning("reading original pvals..")
        
        header =return_header(sumstats)
        columns = [header.index(elem) for elem in column_names]
        it = basic_iterator(sumstats,skiprows=1,columns = columns)
        for elem in it:
            chrom,pos,mlogp,ref,alt,beta,sebeta = elem       
            sum_dict["chr" + '_'.join([chrom,pos,ref,alt])] = [beta,sebeta,mlogp]
        logging.warning("dumping pvals..")
        with open(pval_dict_file, 'wt') as fp:
            json.dump(sum_dict, fp)
        logging.warning('done.')
    else:
        logging.info('reading in json...')
        with open(pval_dict_file) as json_file:
            sum_dict = json.load(json_file)
        logging.info('done.')

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


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Recursive conditional analysis for regenie.")

    parser.add_argument('--pval_threshold',type = float,help ='Threshold limit (pvar -log(mpval) ',default = 6)
    parser.add_argument('--pheno',type = str,help ='Pheno column',required=True)
    parser.add_argument('--out',type = str,help ='Output Directory',required=True)
    parser.add_argument('--covariates',type = str,default = regenie_covariates,help='List of covariates')
    parser.add_argument('--pheno-file',type = file_exists,help ='Path to pheno file',required=True)
    parser.add_argument('--bgen',type = file_exists,help ='Path to bgen',required=True)
    parser.add_argument('--sumstats',type = file_exists,help ='Path to original sumstats',required=True)
    parser.add_argument('--locus',type=str,help ='Initial variant/hit',required=True)
    parser.add_argument('--regenie-params',type=str,help ='extra bgen params',default = '--bt --bsize 200 ' )
    parser.add_argument('--null-file',type = file_exists,help ='File with null info.',required=True)
    parser.add_argument('--force',action = 'store_true',help = 'Flag for forcing re-run.')
    parser.add_argument( "-log",  "--log",  default="warning", choices = log_levels, help=(  "Provide logging level. " "Example --log debug', default='warning'"))
    parser.add_argument('--max-steps',type = int,default =10)
    parser.add_argument('--chr_col', default="#chrom", type=str)
    parser.add_argument('--pos_col', default="pos", type=str)
    parser.add_argument('--ref_col', default="ref", type=str)
    parser.add_argument('--alt_col', default="alt", type=str)
    parser.add_argument('--mlogp_col', default="mlogp", type=str)
    parser.add_argument('--beta', default="beta", type=str)
    parser.add_argument('--sebeta', default="sebeta", type=str)
  

    range_group = parser.add_mutually_exclusive_group(required=True)
    range_group.add_argument('--region',type =str,help ='Region to filter CHR:START-END')
    range_group.add_argument('--chr',type = str,help ='Chromosome to filter on')

    args = parser.parse_args()
    make_sure_path_exists(os.path.dirname(args.out))

    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")
    
       
    if args.region: condition_range = f" --range {args.region}"
    if args.chr: condition_range = f" --chr {chr}"
    if args.pval_threshold <1:args.pval_threshold = -np.log10(args.pval_threshold)
    pretty_print(f"MLOGP THRESHOLD: {args.pval_threshold}")
    # gets original sumstats data for variants
    pval_dict_file = args.out + f"_{args.pheno}_pvals.json"
    columns= [args.chr_col,args.pos_col,args.mlogp_col,args.ref_col,args.alt_col,args.beta,args.sebeta]
    sum_dict = parse_sumstat_data(args.sumstats,pval_dict_file,columns)

    
    # file where to dump results as they come
    result_file = args.out + f"_{args.pheno}_{args.locus}.independent.snps"
    with open(result_file,'wt') as o:
        header = ['VARIANT','BETA','SE','MLOG10P','BETA_cond','SE_cond','MLOG10P_cond','VARIANT_cond']
        o.write("\t".join(header) + '\n')
        #inital values to start looping
        condition_variant = args.locus
        step,exclusion_list = 1,[args.locus]
        while 0 <  step <= args.max_steps:
            out_file = regenie_run(args.out,step,args.bgen,args.pheno_file,args.covariates,condition_variant,args.locus,args.null_file,args.pheno,condition_range,exclusion_list,params = args.regenie_params)
            step,data =check_hit(out_file,step,args.pval_threshold)
            if step:
                new_variant = data[0]
                #print(f"Variant info from original FG sumstats {[round(float(elem),2) for elem in sum_dict[new_variant]]}")
                print(f"Variant info from original FG sumstats {map_vals_to_string(sum_dict[new_variant])}")
                o.write('\t'.join( [new_variant]+ sum_dict[new_variant] + data[1:] +   [','.join(exclusion_list)]) + '\n')               
                condition_variant = new_variant
                exclusion_list.append(new_variant)
                print(f"Hit signficant, proceeding to condition...")
            else:
                print("Hit not significant. Ending loop.")
                


