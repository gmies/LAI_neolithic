# This script reads MOSAIC output (local ancestry) and genotype RData files
# and computes the mean ancestry per SNP across individuals.
# Inputs:
#   1. Local ancestry RData file (localanc)
#   2. Genotype RData file (g.loc)
#   3. Output filename
# Output:
#   - Table with columns: CHR, POS, HGall (mean ancestry per SNP)

library(tidyverse)

args <- commandArgs(TRUE)
localanc_file <- args[1]
geno_file <- args[2]
output_file <- args[3]

# Load MOSAIC output and genotype data
load(localanc_file)  # expects 'localanc' object
load(geno_file)      # expects 'g.loc' object

# Initialize empty dataframe with desired columns
final_dataframe <- data.frame(
  CHR = integer(),
  POS = numeric(),
  HGall = numeric(),
  stringsAsFactors = FALSE
)

# Merge positions from genotype data across all chromosomes
for (i in 1:22) {
  chr_geno <- g.loc[[i]]
  chr_df <- data.frame(CHR = i, POS = chr_geno)
  
  if (nrow(final_dataframe) == 0) {
    final_dataframe <- chr_df
  } else {
    final_dataframe <- merge(final_dataframe, chr_df, by = c("CHR", "POS"), all = TRUE)
  }
}

# Compute mean ancestry per SNP for each chromosome
for (chr in 1:22) {
  chr_lai <- localanc[[chr]]
  subset_array <- chr_lai[2, , ]
  chr_means <- apply(subset_array, 2, mean)
  chr_df <- data.frame(CHR = chr, HGall = chr_means)
  
  # Update final dataframe
  if (nrow(chr_df) == nrow(final_dataframe[final_dataframe$CHR == chr, ])) {
    final_dataframe$HGall[final_dataframe$CHR == chr] <- chr_df$HGall
  } else {
    warning(paste("Row mismatch for chromosome", chr))
  }
}

# Print dimensions and overall mean for sanity check
print(nrow(final_dataframe))
print(mean(final_dataframe$HGall))

# Write output table
write.table(final_dataframe, output_file, quote = FALSE, row.names = FALSE, col.names = TRUE)
