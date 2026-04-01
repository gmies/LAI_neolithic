# This script reads ancestry_hmm posterior files and computes the average ancestry per SNP.
# Inputs:
#   1. Directory containing ancestry_hmm posterior files (*.posterior*)
#   2. Output filename
# Output:
#   - Table with columns: CHR, POS, HGall (average ancestry per SNP)

library(tidyverse)
library(dplyr)

args <- commandArgs(TRUE)
path_directory <- args[1]
output_file <- args[2]

# Change to the directory containing posterior files
setwd(path_directory)

# List posterior files
files <- list.files(pattern = "*.posterior*")
file_count <- length(files)

# Initialize dataframe with chromosome and position from the first file
first_file <- read.table(files[1], header = TRUE)
output_data <- data.frame(CHR = first_file[, 1], POS = first_file[, 2])

# Initialize sum vector
total_sum <- NULL

# Iterate over files to compute average ancestry
for (file in files) {
  data <- read.table(file, header = TRUE)
  
  if (ncol(data) >= 4) {
    # Extract columns 3 and 4 and compute weighted sum
    calculated <- data[[3]] + 0.5 * data[[4]]
    
    # Accumulate sum
    if (is.null(total_sum)) {
      total_sum <- calculated
    } else {
      total_sum <- total_sum + calculated
    }
  } else {
    warning(paste("File", file, "does not have at least 4 columns. Skipping."))
  }
}

# Compute average across files
average_column <- total_sum / file_count
output_data$HGall <- average_column

# Write output table
write.table(output_data, output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")
