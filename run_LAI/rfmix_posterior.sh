#!/bin/bash

# rfmix_posterior.sh
#-------------------------------------------------------------
# Runs RFMix for chromosomes 1-22 and filters posterior 
# probabilities (>0.9) using a local R script.
#-------------------------------------------------------------
# Expected directory structure:
#   ./data/PopPhased/                        # phased input for RFMix
#   ./data/rfmix_alleles1.txt ... rfmix_alleles22.txt
#   ./data/classes.txt
#   ./data/markerLocationsChr1.txt ... markerLocationsChr22.txt
#   ./data/rates.1 ... rates.22              # genetic map / recombination rates
#   ./data/simplai_global_ancestry.txt       # sample global ancestry info
#   ./scripts/posterior_filtering.R          # R script to filter posterior
#-------------------------------------------------------------
# Usage:
#   bash rfmix_posterior.sh

# Load required modules
module load R/4.4
module load bcftools/1.21  
module load plink/2.0-20210505
module load jre/1.8.0_211
module load python/3.8

# Loop over chromosomes 1-22 and run RFMix
for chr in {1..22}
do
    echo "Running RFMix for chromosome $chr..."

    python3 RunRFMix.py data/PopPhased \
        data/rfmix_alleles${chr}.txt \
        data/classes.txt \
        data/markerLocationsChr${chr}.txt \
        -G 35 -n 5 --forward-backward -o output_rfmix_${chr}
done

# Create output folder for filtered files
mkdir -p output
cd output

# Run posterior probability filtering using local R script
Rscript ../scripts/posterior_filtering.R \
    ../data/rates. \
    ../data/simplai_global_ancestry.txt
