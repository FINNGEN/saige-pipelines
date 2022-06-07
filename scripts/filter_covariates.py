import pandas as pd
import numpy as np
import argparse,os
from utils import progressBar,file_exists,make_sure_path_exists

covariates = "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm,BATCH_DS18_MIGRAINE_1_norm,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm,BATCH_DS21_SUPER_2_norm,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm"


def read_pheno_df(pheno_file,phenos,test=True):

    print('reading in pheno df...',end="",flush=True)
    nrows = 100 if test else None
    df= pd.read_csv(pheno_file,nrows = nrows,sep='\t',usecols=phenos).notna().astype(int)
    print('done.')
    return df

def read_cov_df(pheno_file,covariates=covariates.split(','),test=True):

    print('reading in cov df...',end="",flush=True)
    nrows = 100 if test else None
    df= pd.read_csv(pheno_file,sep='\t',usecols=covariates,nrows =nrows).abs()
    print('done')
    return df

def return_df_count(cov_df,pheno_df):
    df = pd.DataFrame()
    columns = len(pheno_df.columns)
    # for each pheno i multply the cov data with the pheno data to get a PHENO x COV df where each entry is the count of non 0 elements in the sample-levl product. This step needs to be done regardless so i can then proceed with all the logic of filtering.
    for i,pheno in enumerate(pheno_df.columns):
        progressBar(i,columns)
        m = (cov_df.mul(pheno_df[pheno],0)>0).sum().to_frame(pheno)
        df = pd.concat([df,m],axis =1)
    progressBar(columns,columns)
    print('\ndone')
    return df


def main(args):

    # read in phenos as list of phenos regardless
    tot_phenos = []
    phenos_groups = []
    with open(args.pheno_list,'rt') as i:
        for line in i:
            phenos = line.strip().split()
            phenos_groups.append(phenos)
            tot_phenos += phenos


    outfile = os.path.join(args.out,"COV_PHENO_COUNTS.txt")
    # here i read in pheno and cov data and merge them into a single df where i calculate the product between the two
    if not os.path.isfile(outfile) or args.force:
        pheno_df = read_pheno_df(args.pheno_file,tot_phenos,args.test)
        # get cov data
        cov_df = read_cov_df(args.pheno_file,args.covariates,test=args.test)
        # do calculation
        df = return_df_count(cov_df,pheno_df)
        # dump to file
        df.to_csv(outfile,sep='\t',header=True,index=True)
    else:
        print(f"{outfile} already generated")
        df = pd.read_csv(outfile,sep='\t',index_col =0)

    print(df.head())
    outfile = os.path.join(args.out,f"pheno_covariates_{args.threshold_count}.txt")
    print(f"Dumping results to {outfile}")
    with open(outfile,'wt') as o,open(outfile.replace('.txt','.err.txt'),'wt') as tmp_err:
        for pheno_list in phenos_groups:
            pheno_name = ','.join(pheno_list)
            # return minimim value across list of phenos for each covariates
            tmp_df = df[pheno_list].min(axis =1)
            # filter out covariates that do not meet threshold count
            covs = tmp_df.index[tmp_df >= args.threshold_count].tolist()
            o.write(f"{pheno_name}\t{','.join(covs)}\n")
            # log rejected covars
            missing_covs = [elem for elem in args.covariates if elem not in covs]
            if missing_covs:tmp_err.write(f"{pheno_name}\t{'.'.join(missing_covs)}\n")
            
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Create count of valid cov entries for each pheno.")
    parser.add_argument('--pheno-file',type = file_exists,help ='Path to pheno file',required=True)
    parser.add_argument('--covariates',type = str,default = covariates,help='List of covariates')
    parser.add_argument('--pheno-list',type = file_exists,help ='File with list of phenos',required=True)
    parser.add_argument('--out',type = str,help ='Output Directory',default = os.getcwd())
    parser.add_argument('--threshold-count','--tc',type =int,help ='Minimum count for covariate',default = 10)

    parser.add_argument('--test',action = 'store_true',help = 'Flag for test run.')
    parser.add_argument('--force',action = 'store_true',help = 'Flag for forcing run run.')

    args = parser.parse_args()
    print(args)
    make_sure_path_exists(os.path.dirname(args.out))
    args.covariates = args.covariates.split(',')
    
    main(args)
