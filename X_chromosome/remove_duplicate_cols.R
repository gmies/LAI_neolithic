#!/usr/bin/env Rscript


#script created 2/9/26
#to take rfmix format with duplicated males and put in format with no duplicates (back in males are haploid)

# Load required library
suppressMessages(library(data.table))

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the input gzipped file as an argument.")
}

input_file <- args[1]

# Read the gzipped file using fread
cat("Reading file:", input_file, "\n")
x <- fread(paste("zcat", input_file), header = FALSE)

# Find duplicated columns
dup_cols <- which(duplicated(as.list(x)))
cat("Number of duplicated columns:", length(dup_cols), "\n")
if (length(dup_cols) > 0) {
  cat("Duplicated columns indices:", dup_cols, "\n")
}

# Remove duplicate columns
x_unique <- x[, !dup_cols, with = FALSE]

# Create output file name
output_file <- sub("(\\.txt)?\\.gz$", "_final\\2.txt.gz", input_file, perl = TRUE)
if (output_file == input_file) {
  output_file <- sub("\\.gz$", "_final.txt.gz", input_file)
}

# Write gzipped output
cat("Writing output to:", output_file, "\n")
fwrite(x_unique, output_file, sep = " ", quote = FALSE, col.name = FALSE)

cat("Done!\n")
