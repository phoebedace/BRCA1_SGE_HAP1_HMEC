# BRCA1_SGE_HAP1_HMEC
Code for processing of SGE datasets and analyses included in 'Saturation genome editing of BRCA1 across cell types accurately resolves cancer risk'.

Please see the manuscript for a description of the datasets and analyses performed. This depository's purpose is to make available all scripts and files used to analyse SGE data for BRCA1.

The order of scripts used for SGE data analysis is:
1. Pipeline.sh - A custom pipeline for converting paired-end, Illumina-based sequencing data into variant counts, extracting editing outcomes, and annotating variants by genomic position using CADD data. Requires input files (in input_files) and custom Python scripts (in custom_python_scripts).
2. Process_pipeline_output_HAP1.Rmd or Process_pipeline_output_HMEC.Rmd – Takes the output of the Pipeline.sh script and does additional processing and creates plots. Must be run for each SGE experiment. Outputs a csv file for each SGE experiment.
3. HAP1_global.R or HMEC_global.R – Takes csv outputs from Process_pipeline_output_HAP1.Rmd or Process_pipeline_output_HMEC.Rmd and merges all data together. Applies filtering and thresholds to determine depleted variants. Data from external sources e.g. ClinVar, computational predictors, is added. Requires some input files from input_files. Outputs include data for Supplementary Tables 1–3.

Additional scripts:
BRCA1_indel_analysis_x1_5UTR.Rmd and BRCA1_indel_analysis_u1_promoter.Rmd produce a file containing function scores for indels in the promoter and 5’ UTR SGE experiments. The output of this file is imported into HAP1_global.R.
