#!/bin/bash
#rfmix format and analysis pipeline (11/18/25 put into script then run on directories)

mkdir rfmix_format
cd rfmix_format

module load R/4.4

mkdir simplai
cd simplai

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/scripts/format/simplai_format.R

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_output_to_mean.R ./output_simplai_ /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates.

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_global.R ./output_simplai_ 

cd ..


mkdir mosaic
cd mosaic

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/scripts/format/mosaic_format.R

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_output_to_mean.R ./output_mosaic_ /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates.

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_global.R ./output_mosaic_ 

cd ..

mkdir ancestry_hmm
cd ancestry_hmm

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/scripts/format/ancestry_hmm_format.R

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_output_to_mean.R ./output_ancestryhmm_ /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates.

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_global.R ./output_ancestryhmm_ 

cd ..

mkdir rfmix
cd rfmix

#format
cp ../../rfmix/output_rfmix_*0.Viterbi.txt.gz .

for f in output_rfmix_*.0.Viterbi.txt.gz; do
  new_name=$(echo "$f" | sed -E 's/\.0\.Viterbi//')
  mv "$f" "$new_name"
done

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_output_to_mean.R ./output_rfmix_ /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates.

Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/scripts/rfmix_global.R ./output_rfmix_

