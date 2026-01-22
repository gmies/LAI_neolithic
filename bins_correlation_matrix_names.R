#9/4/25 updated - with diagnostics

#9/11/24 write script to take different number for correlation:

#8/12/24 run correlation matrix of ancestry hmm results:

library(tidyverse)
library(corrplot)

args <- commandArgs(TRUE)

num_datasets <- as.numeric(args[1])
dataset_paths <- args[2:(1 + num_datasets)]
dataset_names <- args[(2 + num_datasets):(1 + 2 * num_datasets)]

# Get bin size (in kb) as last argument
bin_size_kb <- as.numeric(args[length(args)])
bin_size_bp <- bin_size_kb * 1000  # Convert to base pairs

cat("Bin size input:", bin_size_kb, "kb =", bin_size_bp, "bp\n")

bin_positions <- function(pos, bin_size) {
  return(floor(pos / bin_size) * bin_size)
}

read_and_process <- function(file_path, dataset_id, bin_size) {
  data <- read.table(file_path, header = TRUE)
  data$dataset <- dataset_id
  
  # Add diagnostics
  cat("Dataset", dataset_id, ":\n")
  cat("  Original data points:", nrow(data), "\n")
  cat("  Position range:", min(data$POS), "-", max(data$POS), "\n")
  
  data <- data %>%
    mutate(BIN = bin_positions(POS, bin_size)) %>%
    group_by(CHR, BIN) %>%
    summarise(AVG = mean(HGall), .groups = 'drop')
  
  cat("  After binning:", nrow(data), "bins\n")
  cat("  Unique chromosomes:", length(unique(data$CHR)), "\n")
  
  return(data)
}

# Read and process all datasets
dataset_list <- vector("list", num_datasets)
for (i in 1:num_datasets) {
  dataset_list[[i]] <- read_and_process(dataset_paths[i], dataset_names[i], bin_size_bp)
}

# Merge all datasets
merged_bins <- reduce(dataset_list, 
                       function(x, y) merge(x, y, by = c("CHR", "BIN"), all = TRUE))

# Rename columns for clarity
colnames(merged_bins) <- c("CHR", "BIN", paste0("AVG_", dataset_names))

# More diagnostics
cat("\nMerged data:\n")
cat("Total bins:", nrow(merged_bins), "\n")
cat("Complete cases:", sum(complete.cases(merged_bins)), "\n")

# Show some sample data
cat("\nFirst few rows of merged data:\n")
print(head(merged_bins))

# Compute correlation matrix
correlation_matrix <- cor(merged_bins[, grep("^AVG_", colnames(merged_bins))], use = "complete.obs", method = "pearson")

cat("\nCorrelation matrix:\n")
print(correlation_matrix)


# Define output file name based on bin size
cor_matrix_filename <- paste0("LA_correlation_matrix_binsize_", bin_size_kb, "kb.txt")

# Write correlation matrix to the file
write.table(
  correlation_matrix,
  file = cor_matrix_filename,
  sep = "\t",
  quote = FALSE,
  col.names = NA  # So that row names are preserved
)


# Use the input kb value for filename (not the bp value)
pdf(paste0("LA_correlation_matrix_binsize_", bin_size_kb, "kb.pdf"), width = 10, height = 10)

# Plot correlation matrix
corrplot(correlation_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", # Add correlation coefficients
         diag = FALSE, # Don't show diagonal
         title = paste("Correlation Matrix - Bin Size:", bin_size_kb, "kb"))

dev.off()

cat("Plot saved as: LA_correlation_matrix_binsize_", bin_size_kb, "kb.pdf\n")
