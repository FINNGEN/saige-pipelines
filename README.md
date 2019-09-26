# saige-pipelines


## Conditional analysis for genomewide significant regions.

wdl/saige_conditional_full.wdl and corresponding .json scan for genomewide significant regions and then performs conditional analysis on those regions, adding significant variants as covariate and iterating on that until no significant variants are left.

If you want to run conditional analysis without scanning for gw-sig loci from results files, you can use saige_conditional.wdl/.json directly. It needs configuration file which can be greated using [scripts/generate_conditional_analysis_config.py](scripts/generate_conditional_analysis_config.py). See [scripts/generate_conditional_analysis_config_examples.sh](scripts/generate_conditional_analysis_config_examples.sh) for example commands.



