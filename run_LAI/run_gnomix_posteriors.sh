#!/bin/bash

module load R/4.4

mkdir posterior_output
cd posterior_output

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/7v7/rfmix_format/gnomix/gnomix_posterior_filtering.R /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/GA/simplai_7/simplai_global_ancestry.txt ../../../gnomix 

