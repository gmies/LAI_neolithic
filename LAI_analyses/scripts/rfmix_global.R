#!/usr/bin/env Rscript

# rfmix_global.R
# -------------------------
# Computes per-sample global ancestry proportions from RFMix Viterbi output files.
# Each sample has two columns per SNP; this script calculates the proportion of 2s (ancestry calls)
# and outputs a single mean ancestry per sample.
#
# Inputs:
# 1) Prefix for per-chromosome RFMix output files (e.g., "allentoft_12_10_24_rfmix_")
#     The script will look for files named: <prefix>1.txt.gz, <prefix>2.txt.gz, ..., <prefix>22.txt.gz
# 2) File containing sample names: simplai_global_ancestry.txt
#
# Outputs:
# - output_global_ancestry.txt : table with columns: sample, global_ancestry
#
# Usage:
# Rscript ./scripts/rfmix_global.R <file_prefix>
# Example:
# Rscript ./scripts/rfmix_global.R ./data/allentoft_12_10_24_rfmix_

library(dplyr)

args <- commandArgs(TRUE)

if(length(args) != 1){
  stop("Usage: Rscript rfmix_output_to_mean.R <rfmix_file_prefix>")
}

prefix <- args[1]

# Path to the simplai_global_ancestry.txt file (contains sample names)
simplai_file <- "./data/simplai_global_ancestry.txt"

# Read sample names
simplai_data <- read.table(simplai_file, header = TRUE, sep = " ")
sample_names <- simplai_data$sample  # Assumes first column has sample names

# Initialize empty data frame for combined per-chromosome data
all_data <- data.frame()

# Loop through chromosomes 1-22
for(chr in 1:22){
  
  # Build path to chromosome file
  chr_file <- paste0(prefix, chr, ".txt.gz")
  
  # Read the per-chromosome Viterbi output
  chr_data <- read.table(chr_file, header = FALSE, sep = " ")
  
  # Combine chromosome data into all_data
  all_data <- rbind(all_data, chr_data)
}

# Remove unwanted column if exists
if("V353" %in% colnames(all_data)){
  all_data <- all_data %>% select(-V353)
}

# Ensure number of columns is divisible by 2
num_columns <- ncol(all_data)
if(num_columns %% 2 != 0){
  stop("Number of columns not divisible by 2. Check input files.")
}

# Calculate mean ancestry (proportion of 2s) per sample
column_twos_counts <- numeric()
for(col_idx in seq(1, num_columns, by = 2)){
  col1 <- all_data[[col_idx]]
  col2 <- all_data[[col_idx + 1]]
  
  mean_count <- mean(c(
    sum(col1 == 2, na.rm = TRUE) / nrow(all_data),
    sum(col2 == 2, na.rm = TRUE) / nrow(all_data)
  ))
  
  column_twos_counts <- c(column_twos_counts, mean_count)
}

# Verify sample count matches number of columns
if(length(column_twos_counts) != length(sample_names)){
  stop("Mismatch between number of samples and calculated global ancestry values")
}

# Create final output
results_df <- data.frame(sample = sample_names, global_ancestry = column_twos_counts)

# Write output file
write.table(results_df, "output_global_ancestry.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
