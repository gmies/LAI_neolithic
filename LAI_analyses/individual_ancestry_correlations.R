#!/usr/bin/env Rscript

#started 9/4/25

# Script to calculate per-individual local ancestry correlations across methods
# Usage: Rscript individual_ancestry_correlation.R <num_methods> <method1_dir> <method1_prefix> <method2_dir> <method2_prefix> ... <method1_name> <method2_name> ... <bin_size_kb>

library(tidyverse)
library(corrplot)
library(RColorBrewer)

args <- commandArgs(TRUE)

# Parse arguments
num_methods <- as.numeric(args[1])
method_dirs <- args[2:(1 + num_methods)]
method_prefixes <- args[(2 + num_methods):(1 + 2 * num_methods)]
method_names <- args[(2 + 2 * num_methods):(1 + 3 * num_methods)]
bin_size_kb <- as.numeric(args[length(args)])
bin_size_bp <- bin_size_kb * 1000

cat("Processing", num_methods, "methods with bin size:", bin_size_kb, "kb\n")
cat("Methods:", paste(method_names, collapse = ", "), "\n")

# Load position reference and global ancestry data
pos_file <- "/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/tract_lengths/tractor/mosaic/average_ancestry_output.txt"
global_file <- "/project/mathilab/gmies/neolithic_selection/qpadm/qpadm_by_indiv_082925/allen_imputed/global_ancestry.txt"

pos_data <- read.table(pos_file, header = TRUE)
global_data <- read.table(global_file, header = TRUE)

cat("Loaded", nrow(pos_data), "SNP positions\n")
cat("Loaded", nrow(global_data), "individuals\n")

# Function to bin positions
bin_positions <- function(pos, bin_size) {
  return(floor(pos / bin_size) * bin_size)
}

# Function to read ancestry data for all chromosomes
read_ancestry_method <- function(method_dir, method_prefix, method_name) {
  cat("Reading method:", method_name, "\n")
  
  all_chr_data <- list()
  
  for (chr in 1:22) {
    file_path <- file.path(method_dir, paste0(method_prefix, chr, ".ancestry_only.txt.gz"))
    
    if (!file.exists(file_path)) {
      cat("Warning: File not found:", file_path, "\n")
      next
    }
    
    # Read the ancestry data
    chr_data <- read.table(file_path, header = FALSE)
    
    # Get positions for this chromosome
    chr_pos <- pos_data[pos_data$CHR == chr, ]
    
    if (nrow(chr_data) != nrow(chr_pos)) {
      cat("Warning: Position mismatch for chr", chr, "- expected", nrow(chr_pos), "got", nrow(chr_data), "\n")
      # Take minimum to avoid index errors
      min_rows <- min(nrow(chr_data), nrow(chr_pos))
      chr_data <- chr_data[1:min_rows, ]
      chr_pos <- chr_pos[1:min_rows, ]
    }
    
    # Add position information
    chr_data$CHR <- chr
    chr_data$POS <- chr_pos$POS
    
    all_chr_data[[chr]] <- chr_data
  }
  
  # Combine all chromosomes
  combined_data <- do.call(rbind, all_chr_data)
  cat("  Total SNPs loaded:", nrow(combined_data), "\n")
  
  return(combined_data)
}

# Read all methods
method_data <- list()
for (i in 1:num_methods) {
  method_data[[i]] <- read_ancestry_method(method_dirs[i], method_prefixes[i], method_names[i])
}

# Calculate per-individual ancestry and correlations
num_individuals <- nrow(global_data)
num_haplotypes <- ncol(method_data[[1]]) - 2  # Subtract CHR and POS columns

cat("Number of individuals:", num_individuals, "\n")
cat("Number of haplotypes:", num_haplotypes, "\n")

if (num_haplotypes != num_individuals * 2) {
  cat("Warning: Expected", num_individuals * 2, "haplotypes but found", num_haplotypes, "\n")
}

# Function to calculate individual ancestry averages with binning
calculate_individual_ancestry <- function(ancestry_data, individual_idx, bin_size) {
  # Get the two haplotypes for this individual (columns are 0-indexed after CHR, POS)
  hap1_col <- (individual_idx - 1) * 2 + 1   # +2 for CHR, POS columns
  hap2_col <- (individual_idx - 1) * 2 + 2 
  
  # Create data frame with positions and ancestry
  ind_data <- data.frame(
    CHR = ancestry_data$CHR,
    POS = ancestry_data$POS,
    HAP1 = ancestry_data[, hap1_col],
    HAP2 = ancestry_data[, hap2_col]
  )
  
  # Calculate average ancestry per SNP (HG ancestry = 2, so we want proportion of 2s)
  ind_data$AVG_ANCESTRY <- (ind_data$HAP1 == 2) * 0.5 + (ind_data$HAP2 == 2) * 0.5
  
  # Bin the data
  ind_data$BIN <- bin_positions(ind_data$POS, bin_size)
  
  # Average within bins
  binned_data <- ind_data %>%
    group_by(CHR, BIN) %>%
    summarise(AVG = mean(AVG_ANCESTRY, na.rm = TRUE), .groups = 'drop')
  
  return(binned_data)
}

# Calculate correlations for each individual
correlation_results <- data.frame(
  Individual = global_data$sample,
  Global_Ancestry = global_data$global_ancestry
)

# Add correlation columns for each method pair
method_pairs <- combn(method_names, 2, simplify = FALSE)
pair_names <- sapply(method_pairs, function(x) paste(x[1], x[2], sep = "_vs_"))

for (pair_name in pair_names) {
  correlation_results[[pair_name]] <- NA
}

# Add average correlation column
correlation_results$Average_Correlation <- NA

cat("Calculating correlations for each individual...\n")

for (i in 1:num_individuals) {
  cat("Processing individual", i, "/", num_individuals, "(", global_data$sample[i], ")\n")
  
  # Get binned ancestry data for this individual from each method
  individual_data <- list()
  for (j in 1:num_methods) {
    individual_data[[j]] <- calculate_individual_ancestry(method_data[[j]], i, bin_size_bp)
  }
  
  # Merge all methods for this individual
  merged_individual <- reduce(individual_data, 
                              function(x, y) merge(x, y, by = c("CHR", "BIN"), all = TRUE,
                                                   suffixes = c("", paste0("_", length(individual_data)))))
  
  # Rename columns
  colnames(merged_individual) <- c("CHR", "BIN", paste0("AVG_", method_names))
  
  # Calculate correlations between methods for this individual
  ancestry_cols <- grep("^AVG_", colnames(merged_individual))
  if (length(ancestry_cols) >= 2 && sum(complete.cases(merged_individual[, ancestry_cols])) > 1) {
    
    method_correlations <- cor(merged_individual[, ancestry_cols], use = "complete.obs", method = "pearson")
    
    # Extract pairwise correlations
    pair_correlations <- numeric(length(pair_names))
    for (k in 1:length(method_pairs)) {
      method1_idx <- which(method_names == method_pairs[[k]][1])
      method2_idx <- which(method_names == method_pairs[[k]][2])
      pair_correlations[k] <- method_correlations[method1_idx, method2_idx]
    }
    
    # Store correlations
    correlation_results[i, pair_names] <- pair_correlations
    
    # Calculate average correlation (mean of upper triangle excluding diagonal)
    correlation_results$Average_Correlation[i] <- mean(method_correlations[upper.tri(method_correlations)], na.rm = TRUE)
  }
}

# Save correlation results
write.csv(correlation_results, paste0("individual_correlations_binsize_", bin_size_kb, "kb.csv"), row.names = FALSE)
cat("Correlation results saved to: individual_correlations_binsize_", bin_size_kb, "kb.csv\n")

# Create plots
# Plot 1: Correlations by method pairs
cat("Creating correlation plots...\n")

plot_data_long <- correlation_results %>%
  pivot_longer(cols = all_of(pair_names), 
               names_to = "Method_Pair", 
               values_to = "Correlation") %>%
  filter(!is.na(Correlation))

pdf(paste0("individual_correlations_by_method_binsize_", bin_size_kb, "kb.pdf"), width = 12, height = 8)

p1 <- ggplot(plot_data_long, aes(x = Global_Ancestry, y = Correlation, color = Method_Pair)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.6) +
  labs(x = "Global Ancestry", 
       y = "Local Ancestry Correlation",
       title = paste("Individual Local Ancestry Correlations by Method Pair (Bin Size:", bin_size_kb, "kb)"),
       color = "Method Pair") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 3))

print(p1)
dev.off()

# Plot 2: Average correlations
pdf(paste0("individual_average_correlations_binsize_", bin_size_kb, "kb.pdf"), width = 10, height = 8)

p2 <- ggplot(correlation_results %>% filter(!is.na(Average_Correlation)), 
             aes(x = Global_Ancestry, y = Average_Correlation)) +
  geom_point(alpha = 0.7, size = 2, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.3) +
  labs(x = "Global Ancestry", 
       y = "Average Local Ancestry Correlation",
       title = paste("Individual Average Local Ancestry Correlations (Bin Size:", bin_size_kb, "kb)"),
       subtitle = "Average correlation across all method pairs") +
  theme_minimal()

print(p2)
dev.off()

cat("Plots saved:\n")
cat("  - individual_correlations_by_method_binsize_", bin_size_kb, "kb.pdf\n")
cat("  - individual_average_correlations_binsize_", bin_size_kb, "kb.pdf\n")

# Print summary statistics
cat("\nSummary Statistics:\n")
cat("Number of individuals processed:", sum(!is.na(correlation_results$Average_Correlation)), "/", num_individuals, "\n")
cat("Average correlation range:", round(min(correlation_results$Average_Correlation, na.rm = TRUE), 3), 
    "to", round(max(correlation_results$Average_Correlation, na.rm = TRUE), 3), "\n")
cat("Mean average correlation:", round(mean(correlation_results$Average_Correlation, na.rm = TRUE), 3), "\n")

cat("Analysis complete!\n")
