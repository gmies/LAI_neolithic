# This script formats simpLAI output for downstream analysis:
# correlation analysis, recombination analysis, and Z-score outlier detection.
# Inputs: 
#   1. simpLAI output file (recombination version)
#   2. Number of chromosomes per admixed individual
# Output:
#   - simpLAI_final_output.txt with columns: CHR, POS, HGall

library(tidyverse)
library(dplyr)

args <- commandArgs(TRUE)

lai_file <- args[1]
num_chr_admixed <- as.numeric(args[2])

# Read simpLAI output
lai <- read.table(lai_file, header = TRUE)

# Adjust and normalize values
data <- lai - 1
data$beg <- data$beg + 1

# Compute sum per row and adjust
data$row_sum <- rowSums(data) - data$beg

# Compute mean ancestry proportion for each SNP
data$mean <- 1 - (data$row_sum / num_chr_admixed)
print(mean(data$mean))

# Assign chromosome numbers based on start positions
data$chr <- 1
for (i in 2:length(data$beg)) {
  if (data$beg[i] < data$beg[i - 1]) {
    data$chr[i] <- data$chr[i - 1] + 1
  } else {
    data$chr[i] <- data$chr[i - 1]
  }
}

# Create final output
lai_output <- data.frame(
  CHR = data$chr,
  POS = data$beg,
  HGall = data$mean
)

# Save output to file
write.table(lai_output, "simpLAI_final_output.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
