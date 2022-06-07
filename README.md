# FINNGEN CONDITIONAL ANALYSIS

This is a wrapper pipeline of [regenie](https://rgcgithub.github.io/regenie/) for conditional analysis. The pipeline is mainly built a single python [script](scripts/regenie_conditional.py) that iteratively runs the conditional analysis of regenie until no significant hits are found anymore. The wdl is meant for release purposes and will run all hits from a list of phenos and chromosomes based on the official Finngen results. 

## regenie_conditional.py

This is the "engine" of the pipeline, that can also be used independently, so I will first explain its mechanism and inputs. 

These are the parameters:
```
optional arguments:
  -h, --help            show this help message and exit
  --pval-threshold PVAL_THRESHOLD
                        Threshold limit (-log(mpval))
  --pheno PHENO         Pheno column
  --out OUT             Output Directory and prefix (e.g. /foo/bar/finngen)
  --covariates COVARIATES
                        List of covariates
  --pheno-file PHENO_FILE
                        Path to pheno file
  --bgen BGEN           Path to bgen
  --sample-file SAMPLE_FILE
                        Path to pheno file
  --sumstats SUMSTATS   Path to original sumstats
  --regenie-params REGENIE_PARAMS
                        extra bgen params
  --null-file NULL_FILE
                        File with null info.
  --force               Flag for forcing re-run.
  -log {critical,error,warn,warning,info,debug}, --log {critical,error,warn,warning,info,debug}
                        Provide logging level. Example --log debug',
                        default='warning'
  --max-steps MAX_STEPS
  --chr_col CHR_COL, --chr-col CHR_COL
  --pos_col POS_COL, --pos-col POS_COL
  --ref_col REF_COL, --ref-col REF_COL
  --alt_col ALT_COL, --alt-col ALT_COL
  --mlogp_col MLOGP_COL, --mlogp-col MLOGP_COL
  --beta_col BETA_COL, --beta-col BETA_COL
  --sebeta_col SEBETA_COL, --sebeta-col SEBETA_COL
  --threads THREADS     Number of threads.
  --locus-region LOCUS_REGION LOCUS_REGION
                        Locus & Region to filter CHR:START-END
  --locus-list LOCUS_LIST
                        File with list of locus and regions
```

They are all quite self explanatory. By default all cpus are used and the logging level is set to `warning`. 

The null files are the `*loco.gz` outputs of regenie step 1. The last two inputs are mutually exclusive and are meant for defining the regions of choice. The vanilla mode runs just one region/locus (in any order and in regenie format, e.g. `6:34869517-37869517 chr6_35376598_G_A`). The script will automatically recognize which is the locus and which the region. Else one can pass a file with a tsv separated list of regions/locuses, one per line. The script will then run the main function for each region/locus.

Each run will iteratively condition on more and more significant variants until no hits are found under a certain threshold (`pval-threshold`, either a mlogp > 1 or a pval <1, it gets converted to mglop anyways). One can also choose to cap the iterations at a certain thershold (`max-steps`) instead. 

The outputs will be in the ```--out`` directory (generated if missing). Along with a temporary folder that contains all the necessary files, the outputs are:
- prefix*_pheno_locus.log: the stdout/err of regenie is appended to this file so all logs are available
- prefix*_pheno_locus.independent.snps: contains the chain of results 
- prefix*_pheno_locus_STEP.conditional: the regenie output of each of the [1..N] steps of the chain.

## WDL

Here I will explain the tasks and inputs of the [wdl](wdl/regenie_conditional_full.wdl)