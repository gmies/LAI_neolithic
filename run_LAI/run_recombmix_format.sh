#!/bin/bash

#recombmix to rfmix format


module load R/4.4

mkdir recombmix
cd recombmix

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/7v7/rfmix_format/recombmix_format.R

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_output_to_mean.R ./output_recombmix_ /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates.

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_global.R ./output_recombmix_ 

