#!/bin/bash

#7/7/26

#run recombmix 



vcf_dir="/project/mathilab/gmies/neolithic_selection/allentoft_data/ancestry_prop/data/031825_filtering/davy_snps"
vcf_name="filt_davy.neo."

module load bcftools/1.21

module load htslib/1.21 

#module load boost/1.74.0  

module load gcc/12.3

######################################################
# RECOMB-MIX
######################################################

mkdir recombmix
cd recombmix

for chr in {1..22}
do

echo "Running chromosome ${chr}"

mkdir -p chr${chr}

 /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/RecombMix/Recomb-Mix/RecombMix \
-p ../shared_LAI_inputs/reference.chr${chr}.vcf.gz \
-q ../flare/admixed.chr${chr}.vcf.gz \
-a ../shared_LAI_inputs/recombmix.labels \
-g /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/RecombMix/Recomb-Mix/maps/grch37/genetic_map_GRCh37_chr${chr}.txt \
-o chr${chr} \
-i recombmix_chr${chr}.txt \
-f 0.01 \
-t 8

done

cd ..

