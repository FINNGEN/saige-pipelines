
#!/usr/bin/env python3

import argparse,os.path,shlex,subprocess,sys,subprocess,shlex
from utils import file_exists,make_sure_path_exists


regenie_covariates= "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC{1:10},IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm,BATCH_DS18_MIGRAINE_1_norm,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm,BATCH_DS21_SUPER_2_norm,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm"


def regenie_run(out_root,bgen,pheno_file,covariates,condition_list,pred_file,pheno,region,regenie_cmd = "regenie",params = '--bt --bsize 200 --ref-first'):
    """
    Single regenie run.
    """

    cmd = f'{regenie_cmd} --step 2 {params} --bgen {bgen}  --out {out_root} --pred {pred_file} --phenoFile {pheno_file} --phenoCol {pheno} --condition-list {condition_list} --range {region} --covarFile {pheno_file} --covarColList {covariates}'
    subprocess.call(shlex.split(cmd))
    
    return None




if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Recursive conditional analysis for regenie.")

    parser.add_argument('--region',type = str,help ='CHR:START-END',required=True)
    parser.add_argument('--pheno',type = str,help ='Pheno column',required=True)
    parser.add_argument('--out',type = str,help ='Output Directory',required=True)
    parser.add_argument('--covariates',type = str,default = regenie_covariates,help='List of covariates')
    parser.add_argument('--pheno-file',type = file_exists,help ='Path to pheno file',required=True)
    parser.add_argument('--bgen',type = file_exists,help ='Path to bgen',required=True)
    parser.add_argument('--condition-list',type = file_exists,help ='File with snp(s) to condition on',required=True)
    parser.add_argument('--pred-file',type = file_exists,help ='File with pred info.',required=True)


    args = parser.parse_args()
    

    out_path = os.path.dirname(args.out)
    make_sure_path_exists(out_path)
    regenie_run(args.out,args.bgen,args.pheno_file,args.covariates,args.condition_list,args.pred_file,args.pheno,args.region)
