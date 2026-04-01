#!/usr/bin/env Rscript

# rfmix_output_to_mean.R
# -------------------------
# Computes mean ancestry proportions per SNP across all individuals
# from RFMix Viterbi output files.
#
# Inputs:
# 1) Directory containing RFMix Viterbi outputs: output_rfmix_<chr>.0.Viterbi.txt
# 2) SNP position files (one per chromosome), e.g., rates.1, rates.2, ...
#
# Outputs:
# - average_ancestry_output.txt : table with columns CHR, POS, HGall (mean ancestry)
#
# Usage:
# Rscript ./scripts/rfmix_output_to_mean.R <rfmix_output_dir> <pos_file_prefix>
# Example:
# Rscript ./scripts/rfmix_output_to_mean.R ./data/output_rfmix ./data/rates.

library(dplyr)
library(data.table)

args <- commandArgs(TRUE)

if(length(args) != 2) {
  stop("Usage: Rscript rfmix_average_ancestry.R <rfmix_output_dir> <pos_file_prefix>")
}

rfmix_dir      <- args[1]  # Directory containing RFMix Viterbi files
pos_file_prefix <- args[2] # Prefix for chromosome position files (e.g., rates.)

all_chromosomes_data <- list()

# Loop through chromosomes 1–22
for(chr in 1:22) {
  
  # Construct input file paths
  chr_viterbi_file <- paste0(rfmix_dir, chr, ".0.Viterbi.txt")
  chr_pos_file     <- paste0(pos_file_prefix, chr)
  
  # Read RFMix output
  chr_viterbi <- fread(chr_viterbi_file, header = FALSE, sep = ' ')
  
  # Read SNP positions (skip first line if header info present)
  chr_map <- fread(chr_pos_file, skip = 1)
  
  # Initialize data frame for output
  combined_chr_data <- data.frame(matrix(ncol = 3, nrow = ncol(chr_map)))
  colnames(combined_chr_data) <- c("CHR", "POS", "HGall")
  combined_chr_data$CHR <- chr
  combined_chr_data$POS <- unlist(chr_map[1,])
  
  # Compute mean ancestry for each SNP
  combined_chr_data <- combined_chr_data %>%
    mutate(
      HGall = rowSums(chr_viterbi == 2, na.rm = TRUE) /
              (rowSums(chr_viterbi == 1, na.rm = TRUE) + rowSums(chr_viterbi == 2, na.rm = TRUE))
    )
  
  # Store chromosome-level data
  all_chromosomes_data[[chr]] <- combined_chr_data
}

# Combine all chromosomes
final_output <- bind_rows(all_chromosomes_data)

# Write final output
write.table(final_output, file = "average_ancestry_output.txt", quote = FALSE, row.names = FALSE)
