scripts/generate_conditional_analysis_config.py gs://r3_data/saige/null/R3.grm.v1-{PHENO}.rda gs://r3_data/saige/null/R3.grm.v1-{PHENO}.varianceRatio.txt gs://r3_data/bgen/fromdatateam/{CHR}R3 gs://r3_data/pheno/R3_COV_PHENO_V1.txt.gz r3_conditional.conf vars ../R3/results/all_gw_sig_non_mhc_minimal_annot.tsv

c3392f2f-f7b1-4a32-a0d6-95aa8077ad1a
scripts/generate_conditional_analysis_config.py --chr_col chrom gs://r3_data/downsample/FG/R3.grm.v1-{PHENO}.rda gs://r3_data/downsample/FG/R3.grm.v1-{PHENO}.varianceRatio.txt gs://r3_data/bgen/fromdatateam/{CHR}R3 gs://r3_data/temp/ukb/pheno/FG_ds_cov_pheno.txt.gz FG_downsample_conditional.conf vars ../phenotype_reporting/FINNGEN_demo_downsampled_gwsig.nohla.tsv

scripts/generate_conditional_analysis_config.py  gs://r4_data/nulls_demo/R4_GRM_V1-{PHENO}.rda gs://r4_data/nulls_demo/R4_GRM_V1-{PHENO}.varianceRatio.txt gs://r4_data/bgen/R4_{CHR} gs://r4_data/pheno/R4_COV_PHENO_V1.txt.gz r4_conditional.conf vars /Users/mitja/projects/finngen/phenotype_reporting/R4_all_gw_sig_nohla.tsv

scripts/generate_conditional_analysis_config.py  gs://r3_data/downsample/UKBB/goodQC_white.British.subset.in_Relatedness_variants-{PHENO}.rda gs://r3_data/downsample/UKBB/goodQC_white.British.subset.in_Relatedness_variants-{PHENO}.varianceRatio.txt gs://r3_data/temp/ukb/geno_v3/ukb_imp_chr{CHR}_v3 gs://r3_data/temp/ukb/pheno/UKB_ds_cov_pheno.txt.gz UKBB_demo_downsample_conditional.conf vars ../phenotype_reporting/UKBB_demo_downsampled_gwsig.lifted.nohla.tsv --chr_col chrom

gs://r3_data/temp/ukb/ukb22627_imp_chr1-22_v3_s487395.sample.txt
