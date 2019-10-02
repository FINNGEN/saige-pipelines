# saige-pipelines

## Running GWAS

How to run SAIGE GWAS with Cromwell  
This in an example scenario creating new phenotypes in R4 and running those

1. Create a covariate/phenotype file that contains your phenotypes. E.g. get `gs://r4_data_west1/pheno/R4_COV_PHENO_V1.txt.gz`, add phenotypes to that (cases 1, controls 0, everyone else NA), and upload the new file to a bucket
2. Create a text file with your new phenotypes one per line, e.g.  
    my_phenos.txt
    ```
    PHENO1
    PHENO2
    ```
    and upload the file to a bucket.
3. Clone this repo `git clone https://github.com/FINNGEN/saige-pipelines`
4. Cromwell requires subworkflows be zipped: `cd saige-pipelines/wdl/gwas/ && zip saige_sub saige_sub.wdl`
5. Change `saige.null.phenofile` in `saige.json` to the file from step 1
6. Change `saige.phenolistfile` in `saige.json` to the file from step 2
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
10. Logs and results go to  
`gs://fg-cromwell/saige/WORKFLOW_ID/`

## Conditional analysis for genomewide significant regions.

wdl/saige_conditional_full.wdl and corresponding .json scan for genomewide significant regions and then performs conditional analysis on those regions, adding significant variants as covariate and iterating on that until no significant variants are left.

If you want to run conditional analysis without scanning for gw-sig loci from results files, you can use saige_conditional.wdl/.json directly. It needs configuration file which can be greated using [scripts/generate_conditional_analysis_config.py](scripts/generate_conditional_analysis_config.py). See [scripts/generate_conditional_analysis_config_examples.sh](scripts/generate_conditional_analysis_config_examples.sh) for example commands.



