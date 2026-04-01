#!/bin/bash

module load R/4.4

cp -r /project/mathilab/aaw/Mitonuclear2/data/rfmix/PopPhased/ .
cp /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/7v7/rfmix_format/rfmix/filter_posterior/RunRFMix.py .


module load bcftools/1.21  
module load plink/2.0-20210505
module load jre/1.8.0_211
module load python/3.8

for chr in {1..22}
do
#python2 /project/mathilab/aaw/Mitonuclear2/data/rfmix/RunRFMix.py PopPhased ../../../rfmix/rfmix_alleles${chr}.txt ../../../rfmix/classes.txt /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/rfmix_sep/markerLocationsChr${chr}.txt -G $input_gen -n 5  --forward-backward -o output_rfmix_${chr}

python3 RunRFMix.py PopPhased ../../../rfmix/rfmix_alleles${chr}.txt ../../../rfmix/classes.txt /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/rfmix_sep/markerLocationsChr${chr}.txt -G 35 -n 5  --forward-backward -o output_rfmix_${chr}

done

mkdir output
cd output
Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/scripts/combined_analysis/posterior_filtering.R /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates. /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/GA/simplai_7/simplai_global_ancestry.txt

#Rscript /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/posterior_filtering.R /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/mosaic/rates. /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/GA/simplai_7/simplai_global_ancestry.txt
