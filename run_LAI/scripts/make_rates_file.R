# This script generates recombination rate files for MOSAIC
# by interpolating SNP positions onto a genetic map.
# Inputs:
#   1. Directory containing chromosome SNP files (snpfile.1 … snpfile.22)
# Outputs:
#   - Recombination rate files: rates.1 … rates.22
#     Each file contains:
#       :sites: <num_snps>
#       <space-separated SNP positions>
#       <space-separated cM positions>

library(tidyverse)
library(dplyr)
library(ggplot2)

args <- commandArgs(TRUE)
input_file_dir <- args[1]

for (chr in 1:22) {
  
  # Read chromosome SNP file
  file_path_bim <- paste0(input_file_dir, "snpfile.", chr)
  chr_bim <- read_delim(file_path_bim, delim = "\t", col_names = FALSE)
  colnames(chr_bim) <- c("snp", "chr", "dist", "pos", "a1", "a2")
  
  # Read chromosome-specific genetic map
  map_file <- paste0("../data/shapeit_maps/genetic_map_GRCh37_chr", chr, ".txt")
  chr_map <- read_delim(map_file, delim = " ", col_names = TRUE)
  colnames(chr_map) <- c("bp", "map", "cm")
  
  # Interpolate cM values for SNP positions
  map_fun <- approxfun(chr_map$bp, chr_map$cm)
  chr_bim <- chr_bim %>% mutate(cm = map_fun(pos))
  
  # Prepare output content
  output_content <- c(
    paste(":sites:", nrow(chr_bim)),
    paste(chr_bim$pos, collapse = " "),
    paste(chr_bim$cm, collapse = " ")
  )
  
  # Write output to file
  output_file <- paste0("rates.", chr)
  writeLines(output_content, con = output_file, sep = "\n", useBytes = TRUE)
}
