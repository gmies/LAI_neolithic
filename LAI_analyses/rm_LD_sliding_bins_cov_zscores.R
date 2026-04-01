#sliding bins of 51 snps 
#5/27/25

#bins of 51 snps
#1/24/25

library(tidyverse)

library(corrplot)


library(MASS)
library('dplyr')
library('qqman')
library('ggplot2')


args <- commandArgs(TRUE)



# Number of datasets

num_datasets <- as.numeric(args[1])



# Dataset paths and names

dataset_paths <- args[2:(1 + num_datasets)]

dataset_names <- args[(2 + num_datasets):(1 + 2 * num_datasets)]



###


# Function to compute Z-scores
compute_zscore <- function(data) {
  mean_val <- mean(data$HGall, na.rm = TRUE)
  sd_val <- sd(data$HGall, na.rm = TRUE)
  data$new_z <- (data$HGall - mean_val) / sd_val
  return(data)
}

# Function to read and process a dataset

#and remove long range LD from price paper


# Load long range LD regions once outside the function (global)
ld_regions <- read.table("/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/selection_analysis/same_format_042125/071425_removed/long_range_LD_regions_hg19_5mb.txt", header = FALSE,
                         col.names = c("CHR", "start", "end", "label"))

read_and_process_dataset <- function(file_path, dataset_name) {
  data_input <- read.table(file_path, header = TRUE)
  colnames(data_input) <- c("CHR", "POS", "HGall")
  
  # Filter out positions within long-range LD regions
  filtered_data <- data_input %>% 
    rowwise() %>%
    filter(!any(CHR == ld_regions$CHR & POS >= ld_regions$start & POS <= ld_regions$end)) %>%
    ungroup()
  
  filtered_data <- filtered_data[order(filtered_data$CHR, filtered_data$POS), ]  
  
  binned_data_list <- list()
  
  for (chr in unique(filtered_data$CHR)) {
    chr_data <- filtered_data[filtered_data$CHR == chr, ]
    
    n <- nrow(chr_data)
    window_size <- 51
    sliding_means <- numeric(n)
    
    for (i in 1:n) {
      window_end <- min(i + window_size - 1, n)
      sliding_means[i] <- mean(chr_data$HGall[i:window_end], na.rm = TRUE)
    }
    
    chr_binned_data <- data.frame(
      CHR = chr_data$CHR,
      POS = chr_data$POS,
      HGall = sliding_means,
      dataset = dataset_name
    )
    
    binned_data_list[[as.character(chr)]] <- chr_binned_data
  }
  
  binned_data <- do.call(rbind, binned_data_list)
  
  data <- compute_zscore(binned_data)
  return(data)
}



# Read and process all datasets
all_datasets <- list()
z_scores_matrix <- list()

for (i in 1:num_datasets) {
  dataset_data <- read_and_process_dataset(dataset_paths[i], dataset_names[i])
  all_datasets[[i]] <- dataset_data
  z_scores_matrix[[i]] <- dataset_data$new_z  # Collect Z-scores
}

# Combine Z-scores into a matrix
z_scores_matrix <- do.call(cbind, z_scores_matrix)  # Combine all Z-scores into a matrix

# Compute the covariance matrix of Z-scores
cov_matrix <- cov(z_scores_matrix, use = "pairwise.complete.obs")  # Covariance matrix between Z-scores

# Display covariance matrix (optional, for debugging/insight)
#print(cov_matrix)

# Calculate weighted Z-scores using the formula:
# Z' = sum(Z_i) / sqrt(sum_ij R_ij), where R_ij is the covariance matrix

# For each SNP, calculate the weighted Z-score
combined_data <- bind_rows(all_datasets) %>%
  group_by(CHR, POS) %>%
  summarise(AVG_Z = sum(new_z, na.rm = TRUE), .groups = 'drop')

# Compute the denominator: sqrt(sum_ij R_ij)
denominator <- sqrt(sum(cov_matrix))

# Compute the final weighted Z-scores
combined_data$weighted_z <- combined_data$AVG_Z / denominator

# Create a combined plot for the weighted Z-scores across all datasets
combined_data <- combined_data %>%
  arrange(CHR, POS) %>%
  mutate(index = 1:n()) %>%
  mutate(p_value = 2 * (1 - pnorm(abs(weighted_z))))

# Save combined Z-scores to file
write.table(combined_data, "combined_weighted_z_scores.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)



combined_data$color <- ifelse(as.numeric(combined_data$CHR) %% 2 == 0, "red", "blue")



pdf("combined_avg_zscore_plot.pdf", width = 9, height = 6)



ggplot(combined_data, aes(x = index, y = weighted_z, color = color)) +

  geom_point(size = 0.5) +  # Add points

  labs(x = "Chromosome", y = "Average Z-score") +  # Axis labels

  theme_minimal() +  # Minimal theme 

  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at y = 0

geom_hline(yintercept = 3, color = "blue") +
geom_hline(yintercept = 3.5, color = "blue") +
geom_hline(yintercept = -3, color = "blue") +
geom_hline(yintercept = -3.5, color = "blue") +

  scale_color_identity() +  # Use colors defined in the color column

  theme(axis.text.x = element_blank(),  # Remove x-axis text

        axis.ticks.x = element_blank())



dev.off()



qq_data <- combined_data %>%
  dplyr::select(p_value = p_value, chromosome = CHR, POS) %>%
  dplyr::mutate(chromosome = as.numeric(chromosome), POS = as.numeric(POS), p_value = as.numeric(p_value))

pdf("combined_avg_qq.pdf")

qq(qq_data$p_value)

dev.off()
