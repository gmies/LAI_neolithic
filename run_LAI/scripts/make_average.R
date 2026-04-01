# This script processes FLARE ancestry output by individual to compute
# the proportion of minor ancestry per SNP across all individuals.
# Inputs:
#   1. FLARE by-individual output file 
#   2. Desired output filename
# Output:
#   - Data frame with columns: CHR, POS, HGall (mean minor ancestry proportion)

library(dplyr)
library(tidyverse)
library(stringr)

args <- commandArgs(TRUE)
flare_input_file <- args[1]
output_file_name <- args[2]

# Read FLARE output
flare_input <- read.table(
  flare_input_file,
  header = FALSE,
  colClasses = c("numeric", "numeric", "character", "character")
)

# Extract columns representing haplotype ancestry
col3 <- as.character(flare_input$V3)
col4 <- as.character(flare_input$V4)

# Count ones (minor ancestry) in both columns per row
count_ones_col3 <- str_count(col3, "1")
count_ones_col4 <- str_count(col4, "1")
total_ones_per_row <- count_ones_col3 + count_ones_col4

# Count total number of alleles (0s + 1s) per row
total_entries_col3 <- nchar(col3)
total_entries_col4 <- nchar(col4)
total_entries_per_row <- total_entries_col3 + total_entries_col4

# Compute minor ancestry proportion per SNP
ratio <- total_ones_per_row / total_entries_per_row

# Prepare output data frame
flare_output <- data.frame(
  CHR = flare_input$V1,
  POS = flare_input$V2,
  HGall = ratio
)

# Save output to file
write.table(flare_output, output_file_name, quote = FALSE, row.names = FALSE, col.names = TRUE)
