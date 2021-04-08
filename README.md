# saige-pipelines

## Running GWAS

How to run SAIGE GWAS with Cromwell  
This in an example scenario creating new phenotypes in R7 and running those

1. Create a covariate/phenotype file that contains your phenotypes. E.g. get `gs://r7_data/pheno/R7_COV_PHENO_V1.txt.gz`, add phenotypes to that (cases 1, controls 0, everyone else NA), and upload the new file to a bucket
2. Create a text file with your new phenotypes one per line, e.g.  
    my_phenos.txt
    ```
    PHENO1
    PHENO2
    ```
    and upload the file to a bucket.
3. Clone this repo `git clone https://github.com/FINNGEN/saige-pipelines`
4. Cromwell requires subworkflows be zipped: `cd saige-pipelines/wdl/gwas/ && zip saige_sub saige_sub.wdl saige_summary.wdl`
5. Change `saige.null.phenofile` in `saige.json` to the file from step 1
6. Change `saige.phenolistfile` in `saige.json` to the file from step 2  
    6.1. Use `"saige.traitType": "binary"` or `"saige.traitType": "quantitative"` depending on whether your traits are case/control or continuous  
    6.2. Use `"saige.analysisType": "additive"` or `"saige.analysisType": "recessive"`, `"saige.analysisType": "dominant"` or `"saige.analysisType": "het"` - additive being regular GWAS. 
7. Connect to Cromwell server  
    `gcloud compute ssh cromwell-fg-1 --project finngen-refinery-dev --zone europe-west1-b -- -fN -L localhost:5000:localhost:80`
8. Submit workflow  
    8.1. Using the web interface  
        8.1.1 Go to `http://0.0.0.0:5000` with your browser  
        8.1.2 Click `/api/workflows/{version}`  
        8.1.3 Choose `wdl/gwas/saige.wdl` as workflowSource  
        8.1.4 Choose the edited `wdl/gwas/saige.json` as workflowInputs  
        8.1.5 Choose `wdl/gwas/saige_sub.zip` as workflowDependencies  
        8.1.6 `Execute`  
    8.2. Or with `https://github.com/FINNGEN/CromwellInteract`
9. Use the given workflow id to look at timing diagram or to get metadata  
`http://0.0.0.0:5000/api/workflows/v1/WORKFLOW_ID/timing`
`http://0.0.0.0:5000/api/workflows/v1/WORKFLOW_ID/metadata`
10. Logs and results go under  
`gs://fg-cromwell/saige/WORKFLOW_ID`, plots `gs://fg-cromwell/saige/WORKFLOW_ID/call-test_combine/shard-*/**/*.png`, summary stats and tabix indexes `gs://fg-cromwell/saige/WORKFLOW_ID/call-test_combine/shard-*/**/*.gz*`

## Docker file creation for R5 GWAS

```
git clone https://github.com/FINNGEN/saige-pipelines
cd saige-pipelines
git clone https://github.com/weizhouUMICH/SAIGE -b finngen_r6_jk
docker build -t gcr.io/finngen-refinery-dev/saige:0.39.1-TAG -f docker/Dockerfile_SAIGE_GWAS .
```

## Conditional analysis for genomewide significant regions.

wdl/saige_conditional_full.wdl and corresponding .json scan for genomewide significant regions and then performs conditional analysis on those regions, adding significant variants as covariate and iterating on that until no significant variants are left.

If you want to run conditional analysis without scanning for gw-sig loci from results files, you can use saige_conditional.wdl/.json directly. It needs configuration file which can be greated using [scripts/generate_conditional_analysis_config.py](scripts/generate_conditional_analysis_config.py). See [scripts/generate_conditional_analysis_config_examples.sh](scripts/generate_conditional_analysis_config_examples.sh) for example commands.

## Output files

### PHENOTYPE.REGION.independent.snps files
#### Summary of top snp conditional statistics after conditioning

Columns:

- SNPID   variant
- BETA    original beta
- SE      original se
- p.value original se
- BETA_cond       beta after conditioning
- SE_cond se after conditioning
- p.value_cond    p-value after conditioning
- Conditioned_on variants used in conditioning 


### PHENOTYPE.REGION_n files
#### Summary statistics of all snps in the region after conditioning on n snps.
#### The condition snp corresponds to the line in .independent.snps file

Columns:
- CHR 
- POS 
- rsid 
- SNPID 
- Allele1 
- Allele2 
- AC_Allele2 
- AF_Allele2 
- imputationInfo 
- N 
- BETA 
- SE 
- Tstat 
- p.value  original p-value using SPA estimator (you want this p-value)
- p.value.NA  original p-value using normal approximation (don’t use!)
- Is.SPA.converge  did the model converge
- varT original t statistic
- varTstat  original variance of t-statistic
- Tstat_cond  t-statistic after conditioning
- p.value_cond  p-value after conditioning
- varT_cond  t statistic variance after conditioning
- BETA_cond  beta after conditioning
- SE_cond standard error after conditioning



