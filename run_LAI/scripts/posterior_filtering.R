# Script to filter RFMix posterior probabilities for Neolithic data
# Computes per-SNP and per-individual average ancestry probabilities
# Outputs:
# 1) average_ancestry_output.txt  (per-SNP)
# 2) rfmix_global_ancestry.txt    (per-individual)

library(dplyr)
library(data.table)

args <- commandArgs(TRUE)

txt_file_path <- args[1]      # Directory containing position (.map/.txt) files per chromosome
simplai_file <- args[2]       # File containing sample names (simplai_global_ancestry.txt)

# Initialize empty dataframe to store per-individual posterior probabilities
all_data <- data.frame()

# Read sample names from simplai file
simplai_data <- read.table(simplai_file, header = TRUE, sep = " ")
sample_names <- simplai_data$sample

# List to store per-chromosome results
all_chromosomes_data <- list()

# Loop through chromosomes 1-22
for (i in 1:22) {
  
  # Read RFMix ForwardBackward output for this chromosome
  chr_specific_txt_name <- paste0("../output_rfmix_", i, ".0.ForwardBackward.txt")
  chr_specific_txt <- data.table::fread(chr_specific_txt_name, header = FALSE, sep = ' ')
  
  # Copy for clarity
  chr_specific_txt_correct_indiv <- chr_specific_txt
  
  # Read position file for this chromosome
  chr_specific_map_name <- paste0(txt_file_path, i)
  chr_specific_map <- fread(chr_specific_map_name, skip = 1)
  
  # Initialize dataframe for posterior probabilities
  post_dataframe <- data.frame(matrix(ncol = (ncol(chr_specific_txt)/2), nrow = nrow(chr_specific_txt)))
  
  # Take every other column (odd columns) as relevant probabilities
  data_odd <- chr_specific_txt[, seq(1, ncol(chr_specific_txt), by = 2), with = FALSE]
  
  # Assign 1 for high-confidence (>0.9) and 2 for low-confidence (<0.1)
  post_dataframe[data_odd > 0.9] <- 1
  post_dataframe[data_odd < 0.1] <- 2    
    
  # Compute per-SNP ancestry proportions
  combined_chr_specific_data <- data.frame(matrix(ncol = 3, nrow = ncol(chr_specific_map)))
  colnames(combined_chr_specific_data) <- c('CHR', 'POS', 'HGall')
  combined_chr_specific_data$CHR <- i
  combined_chr_specific_data$POS <- unlist(chr_specific_map[1,])
  
  combined_chr_specific_data <- combined_chr_specific_data %>%
    mutate(
      HGall = rowSums(post_dataframe == 2, na.rm = TRUE) / 
              (rowSums(post_dataframe == 1, na.rm = TRUE) + rowSums(post_dataframe == 2, na.rm = TRUE))
    )
  
  # Store per-chromosome results
  all_chromosomes_data[[i]] <- combined_chr_specific_data
  
  # Append per-individual data for global ancestry calculation
  all_data <- rbind(all_data, post_dataframe)
}

# Combine per-SNP results and write output
final_output <- bind_rows(all_chromosomes_data)
write.table(final_output, file = "average_ancestry_output.txt", quote = FALSE, row.names = FALSE)

# Compute global ancestry per individual
num_columns <- ncol(all_data)
if (num_columns %% 2 != 0) stop("Number of columns not divisible by 2")

column_twos_counts <- numeric()

for (col_idx in seq(1, num_columns, by = 2)) {
  col1 <- all_data[[col_idx]]
  col2 <- all_data[[col_idx + 1]]
  
  count_twos_col1 <- sum(col1 == 2, na.rm = TRUE) / nrow(all_data)
  count_twos_col2 <- sum(col2 == 2, na.rm = TRUE) / nrow(all_data)
  
  column_twos_counts <- c(column_twos_counts, mean(c(count_twos_col1, count_twos_col2)))
}

if (length(column_twos_counts) != length(sample_names)) stop("Mismatch in number of samples")

results_df <- data.frame(sample = sample_names, global_ancestry = column_twos_counts)
write.table(results_df, "rfmix_global_ancestry.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
