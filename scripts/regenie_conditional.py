
#!/usr/bin/env python3

import argparse,os.path,shlex,subprocess,sys,subprocess,shlex
from utils import file_exists,make_sure_path_exists,tmp_bash,pretty_print
import pandas as pd

regenie_covariates= "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC{1:10},IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm,BATCH_DS18_MIGRAINE_1_norm,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm,BATCH_DS21_SUPER_2_norm,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm"


def regenie_run(out_root,step,bgen,pheno_file,covariates,condition_variant,original_hit,null_file,pheno,region,regenie_cmd = "regenie",params = '--bt --bsize 200 ',):
    """
    Single regenie conditainal run. 
    Returns file with results after doing some renaming.
    """

    regenie_file = out_root + f"_{pheno}.regenie" #default regenie output
    out_file = out_root + f"_{pheno}_{original_hit}_{step}.conditional" #file where to redirect output
    pretty_print(f"VARIANT:{condition_variant}")
    print(f"generating {out_file} ...",end = "", flush=True)
 
    if not os.path.isfile(out_file) or args.force:
        args.force = True

        #build pred
        pred_file = out_root + f"_{args.pheno}.pred"
        with open(pred_file,'wt') as o: o.write(f"{pheno}\t{null_file}")

        # build condition file
        tmp_variant = out_root + '.variant.tmp'
        with open(tmp_variant,'wt') as o: o.write(condition_variant)
        
        cmd = f'{regenie_cmd} --step 2 {params} --bgen {bgen}  --out {out_root} --pred {pred_file} --phenoFile {pheno_file} --phenoCol {pheno} --condition-list {tmp_variant} --range {region} --covarFile {pheno_file} --covarColList {covariates}'
        log_mode = 'wt' if step ==1 else 'a'
        with open(out_root+ f"_{args.pheno}_{original_hit}.log",log_mode) as o:
            subprocess.call(shlex.split(cmd),stdout=o)
        print('done.')

        # spring cleaning
        for f in [pred_file,tmp_variant,out_root + '.log'] : os.remove(f)
        os.replace(regenie_file,out_file)

    else:
        print('file already exists')     

    return out_file

def parse_sumstat_data(sumstats):
    """
    Reads in sumstats as dictionary for relevant fields.
    """
    sum_dict = {}
    with open(sumstats) as i:
        line= next(i)
        for entry in i:
            chrom,pos,ref,alt,pval,mlogp,beta,sebeta,*_ = entry.strip().split()
            sum_dict["chr" + '_'.join([chrom,pos,ref,alt])] = [beta,sebeta,mlogp]
    return sum_dict

def check_hit(out_file,step,threshold=7):
    """
    Checks if regenie hit is significant
    """

    # read results as pandas df and get row with max -log10(p)
    df = pd.read_csv(out_file,sep=" ")
    # get max value,extract relevant fields and map to str
    data = list(map(str,df.iloc[df['LOG10P'].idxmax()].loc[['ID','BETA','SE','LOG10P']].values))
    variant,beta,se,pval = data
    print(f"CANDIDATE VARIANT:{variant}")
    print(f"-log(pval): {pval},threshold:{threshold}")
    step = step+ 1 if  float(pval) > threshold else 0
    print(data)
    return step,data

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Recursive conditional analysis for regenie.")

    parser.add_argument('--region',type = str,help ='CHR:START-END',required=True)
    parser.add_argument('--pval_threshold',type = int,help ='-log(pval) limit ',default = 6)
    parser.add_argument('--pheno',type = str,help ='Pheno column',required=True)
    parser.add_argument('--out',type = str,help ='Output Directory',required=True)
    parser.add_argument('--covariates',type = str,default = regenie_covariates,help='List of covariates')
    parser.add_argument('--pheno-file',type = file_exists,help ='Path to pheno file',required=True)
    parser.add_argument('--bgen',type = file_exists,help ='Path to bgen',required=True)
    parser.add_argument('--sumstats',type = file_exists,help ='Path to original sumstats',required=True)
    parser.add_argument('--initial-hit',type=str,help ='Initial variant/hit',required=True)
    parser.add_argument('--null-file',type = file_exists,help ='File with null info.',required=True)
    parser.add_argument('--force',action = 'store_true',help = 'Flag for forcing re-run.')


    args = parser.parse_args()
    make_sure_path_exists(os.path.dirname(args.out))

    # gets original sumstats data for variants
    sum_dict = parse_sumstat_data(args.sumstats)
    
    # file where to dump results as they come
    result_file = args.out + f"_{args.pheno}_{args.initial_hit}.independent.snps"
    with open(result_file,'wt') as o:
        header = ['VARIANT','BETA','SE','MLOG10P','BETA_cond','SE_cond','MLOG10P_cond','VARIANT_cond']
        o.write("\t".join(header) + '\n')
        #inital values to start looping
        condition_variant = args.initial_hit
        step,results = 1,[]
        while step :
            out_file = regenie_run(args.out,step,args.bgen,args.pheno_file,args.covariates,condition_variant,args.initial_hit,args.null_file,args.pheno,args.region)
            step,data =check_hit(out_file,step,args.pval_threshold)
            if step:
                o.write('\t'.join( [data[0]]+ sum_dict[data[0]] + data[1:] +   [condition_variant]) + '\n')               
                condition_variant = data[0]
            else:
                print("Hit not significant. Ending loop.")
                
