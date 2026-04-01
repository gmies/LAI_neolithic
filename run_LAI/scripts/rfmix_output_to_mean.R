# This script reads RFMix Viterbi output and calculates average ancestry per SNP for each chromosome.
# Inputs:
#   1. Directory with RFMix Viterbi output files (1.0.Viterbi.txt for each chr)
#   2. Directory with map files for positions (one file per chromosome)
# Output:
#   - average_ancestry_output.txt: CHR, POS, HGall (mean ancestry per SNP)

library(dplyr)
library(data.table)

args <- commandArgs(TRUE)
dir_file_path <- args[1]  # Directory containing RFMix output
txt_file_path <- args[2]  # Directory containing map files

all_chromosomes_data <- list()

for (i in 1:22) {
  
  # Read chromosome-specific Viterbi output
  chr_txt_file <- paste0(dir_file_path, i, ".0.Viterbi.txt")
  chr_txt <- fread(file = chr_txt_file, header = FALSE, sep = ' ')
  
  # Read chromosome-specific map file (skip first line)
  chr_map_file <- paste0(txt_file_path, i)
  chr_map <- fread(chr_map_file, skip = 1)
  
  # Initialize dataframe for this chromosome
  chr_data <- data.frame(
    CHR = i,
    POS = unlist(chr_map[1,]),
    HGall = NA_real_
  )
  
  # Compute mean ancestry per SNP
  chr_data$HGall <- rowSums(chr_txt == 2, na.rm = TRUE) / 
                    (rowSums(chr_txt == 1, na.rm = TRUE) + rowSums(chr_txt == 2, na.rm = TRUE))
  
  # Store chromosome data
  all_chromosomes_data[[i]] <- chr_data
}

# Combine all chromosomes into one dataframe
final_output <- bind_rows(all_chromosomes_data)

# Write output
write.table(final_output, file = "average_ancestry_output.txt", quote = FALSE, row.names = FALSE)
