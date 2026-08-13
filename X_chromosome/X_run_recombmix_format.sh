#!/bin/bash

#recombmix to rfmix format


module load R/4.4

Rscript recombmix_format.R

Rscript ../mosaic_7v7/rfmix_output_to_mean.R ./output_recombmix_ ../mosaic_7v7/rates.

Rscript /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/remove_duplicate_cols.R output_recombmix_23.txt.gz

Rscript ../mosaic_7v7/rfmix_output_to_mean.R ./output_recombmix_23_final ../mosaic_7v7/rates.
Rscript ../mosaic_7v7/rfmix_global.R ./output_recombmix_ 

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 average_ancestry_output.txt

python  /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/scripts/analyze_rfmix_tracts.py \
--hmm-file /project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/ancestry_hmm_7v7/hap/output_hmm_run_file \
--rfmix-pattern "output_recombmix_\${chr}_final.txt.gz" \
  --output-prefix recombmix_chr23 \
  --chr-start 23 \
  --chr-end 23


#Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/replication/scripts/sliding_bins_cov_zscores.R 1 average_ancestry_output.txt

