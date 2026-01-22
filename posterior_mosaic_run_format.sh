#!/bin/bash
#rfmix format and analysis pipeline (11/18/25 put into script then run on directories)


module load R/4.4


Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/scripts/format/posterior_mosaic_format.R

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_output_to_mean.R ./output_mosaic_ /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates.

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_global.R ./output_mosaic_ 

