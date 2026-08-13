#!/bin/bash

#7/7/26

#run recombmix 



vcf_dir="/project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/data"
#vcf names
vcf_name="chrX.phased."

module load bcftools/1.21

module load htslib/1.21 

#module load boost/1.74.0  

module load gcc/12.3

######################################################
# RECOMB-MIX
######################################################

#mkdir recombmix
cd recombmix

chr=X

echo "Running chromosome ${chr}"

mkdir -p chr${chr}

 /project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/results/11_12_25_run_samplesizes/RecombMix/Recomb-Mix/RecombMix \
-p reference.chrX.phased.vcf.gz \
-q admixed.chrX.phased.vcf.gz \
-a ../shared_LAI_inputs/recombmix.labels \
-g ../shared_LAI_inputs/interpolated_maps/genetic_map_GRCh37_chrX.txt \
-o chrX \
-i recombmix_chrX.txt \
-f 0.01 \
-t 8


cd ..

