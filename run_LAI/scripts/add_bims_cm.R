# add_bims_cm.R
#-------------------------------------------------------
# Purpose: Read in per-chromosome BIM files, 
#          interpolate genetic map distances (cM), 
#          and output a combined BIM file with cM positions.
# Input: Directory with BIM files for chromosomes 1-22
# Output: final_allentoft_bim.bim in working directory
#-------------------------------------------------------

library(tidyverse)

# Get input directory from command-line argument
args <- commandArgs(TRUE)
input_file_dir <- args[1]

# Initialize empty tibble for concatenated BIM data
large_bim_file <- tibble()

# Loop over chromosomes 1-22
for (chr in 1:22) {
  
  # Read BIM file for this chromosome
  file_path_bim <- paste0(input_file_dir, "full_dataset_chr", chr, ".bim")
  chr_bim <- read_delim(file_path_bim, delim = "\t", col_names = FALSE)
  
  # Read chromosome-specific genetic map (GRCh37)
  file_path <- paste0("./data/genetic_map_GRCh37_chr", chr, ".txt")
  chr_specific_recombination <- read_delim(file_path, delim = " ", col_names = TRUE)
  colnames(chr_specific_recombination) <- c("bp", "map", "cm")
  
  # Interpolate genetic map (cM) for each SNP position
  map_fun <- approxfun(chr_specific_recombination$bp, chr_specific_recombination$cm)
  chr_bim <- chr_bim %>%
    mutate(X3 = map_fun(X4))
  
  # Combine this chromosome's BIM with the full dataset
  large_bim_file <- bind_rows(large_bim_file, chr_bim)
}

# Write combined BIM file with cM positions
write.table(large_bim_file, "final_allentoft_bim.bim", quote = FALSE, row.names = FALSE, col.names = FALSE)
