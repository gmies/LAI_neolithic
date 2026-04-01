#!/bin/bash

#script to take ancestry hmm output and put into rfmix format and then calculate LA for admixed pop and global ancestry for each admixed individual 

module load R/4.4
Rscript LAI_analyses/scripts/filter_ancestry_hmm_format.R
Rscript LAI_analyses/scripts/rfmix_output_to_mean.R ./output_ancestryhmm_ /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates.
Rscript LAI_analyses/scripts/rfmix_global.R ./output_ancestryhmm_ 


