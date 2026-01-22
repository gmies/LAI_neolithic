
#3/26/25 to filter rfmix on neolithic data for posterior probabilites 


library(dplyr)
library(data.table)

args <- commandArgs(TRUE)


#arg 1 is pos file
txt_file_path <- args[1]
#ex: ../mosaic/rates.


#sample name list:
simplai_file <- args[2]
#ex: ../simplai_7/simplai_global_ancestry.txt

all_data <- data.frame()

# Read in the simplai_global_ancestry.txt file to get the sample names
simplai_data <- read.table(simplai_file, header = TRUE, sep = " ")
sample_names <- simplai_data$sample




#read in files that have output_rfmix_*ForwardBackward.txt 



f = 0

all_chromosomes_data <- list()


for (i in 1:22) {
  
  #read in file for each chr:
  chr_specific_txt_name <- paste0("../output_rfmix_", i, ".0.ForwardBackward.txt")
  
  
  chr_specific_txt <- data.table::fread(file = chr_specific_txt_name, header = FALSE, sep = ' ')
  
  
  #renaming (to hold name from old script)
  chr_specific_txt_correct_indiv <- chr_specific_txt
  
  
  #for naming pos:
    chr_specific_map_name <- paste0(txt_file_path, i)
    chr_specific_map <- fread(chr_specific_map_name, skip = 1)
  
  
  
  #make dataframe based on posterior probabilites:
  
    post_dataframe <- data.frame(matrix(ncol = (ncol(chr_specific_txt)/2), nrow = nrow(chr_specific_txt)))
                                

data_odd <-chr_specific_txt[, seq(1, ncol(chr_specific_txt), by = 2), with = FALSE]



post_dataframe[data_odd > 0.9] <- 1

post_dataframe[data_odd < 0.1] <- 2    
    
  
  
  
  
  #for by snp:
  
  combined_chr_specific_data <- data.frame(matrix(ncol = 3, nrow = ncol(chr_specific_map)))
  
  x <- c('CHR', 'POS', 'HGall')
  colnames(combined_chr_specific_data) <- x
  
  combined_chr_specific_data$CHR <- i
  combined_chr_specific_data$POS <- unlist(chr_specific_map[1,])
  
  
 #trying this instead of for loop that is used in previous script:
  
  
  combined_chr_specific_data <- combined_chr_specific_data %>%
  mutate(
    HGall = rowSums(post_dataframe == 2, na.rm = TRUE) / (rowSums(post_dataframe == 1, na.rm = TRUE) + rowSums(post_dataframe == 2, na.rm = TRUE))
  )
  

  
  #for each chromosome, make a new dataframe for just that chromosome
  chr_specific_ancestry <- combined_chr_specific_data
  
  # Append the chromosome data to the list
  all_chromosomes_data[[i]] <- combined_chr_specific_data  
  
  
  
  
  
  #for by individual:
  
  
  all_data <- rbind(all_data, post_dataframe)
  
  
  
  
  
}




#finish by snp:

final_output <- bind_rows(all_chromosomes_data)

output_name <- paste0("average_ancestry_output.txt")

write.table(final_output, file = output_name, quote = FALSE, row.names = FALSE)




#finish by indiv:



# Check if the number of columns is divisible by 2
num_columns <- ncol(all_data)
if (num_columns %% 2 != 0) {
  stop("The number of columns is not divisible by 2. Please check your input data.")
}

# Initialize an empty vector to store counts of twos for every two columns
column_twos_counts <- numeric()

# Loop through every two columns at a time and calculate the count of 2s for each pair
for (col_idx in seq(1, num_columns, by = 2)) {
  # Take the two columns (i.e., data from the two samples)
  col1 <- all_data[[col_idx]]
  col2 <- all_data[[col_idx + 1]]
  
  # Count the number of 2s in the first and second column (i.e., the number of 2s for each sample)
  count_twos_col1 <- sum(col1 == 2, na.rm = TRUE) / nrow(all_data)
  count_twos_col2 <- sum(col2 == 2, na.rm = TRUE) / nrow(all_data)
  


  # Calculate the mean count of 2s for the pair of columns
  mean_count_twos <- mean(c(count_twos_col1, count_twos_col2))
  
  # Append the mean value to the column_twos_counts vector
  column_twos_counts <- c(column_twos_counts, mean_count_twos)
}

# Check if the length of column_twos_counts matches the number of samples
if (length(column_twos_counts) != length(sample_names)) {
  stop("The number of column twos counts does not match the number of samples. Please check your input data.")
}

# Now create the final data frame with sample names and calculated global ancestry (count of 2s)
results_df <- data.frame(sample = sample_names, global_ancestry = column_twos_counts)

# Write the output to a file
write.table(results_df, "rfmix_global_ancestry.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)




