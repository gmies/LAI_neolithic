#11/18/24 script adapted from below to run mean output for rfmix output based on arg input of files


#script to run average ancestry for RFMix calls from PMBB AA
#10/4/24
#script adapted from the one I used to run for Goldberg Cabo Verde populations: /project/mathilab/gmies/prelim_results/correlations/Goldberg_results/scripts/2_27_24_rec_ancestry_correlation.R



library(dplyr)
library(data.table)

args <- commandArgs(TRUE)

#arg 1 is the directory where to find rfmix output files
dir_file_path <- args[1]

#arg 2 is pos file
txt_file_path <- args[2]


f = 0

all_chromosomes_data <- list()

i=23
  
  #read in file for each chr:
#  chr_specific_txt_name <- paste0(dir_file_path, i, ".txt.gz")
  chr_specific_txt_name <- paste0(dir_file_path, ".txt.gz")  
  chr_specific_map_name <- paste0(txt_file_path, i)
  
  
  chr_specific_txt <- data.table::fread(file = chr_specific_txt_name, header = FALSE, sep = ' ')
                                          
 # chr_specific_map <- read.table(chr_specific_map_name, header = FALSE, sep = ' ', skip = 1)
  
 chr_specific_map <- fread(chr_specific_map_name, skip = 1)

  #renaming (to hold name from old script)
  chr_specific_txt_correct_indiv <- chr_specific_txt
  

  
  #make new dataframe with the three columns of:  chr (from the input), bp (from the .map first column), average of the rows for that snp in a third column
  combined_chr_specific_data <- data.frame(matrix(ncol = 3, nrow = ncol(chr_specific_map)))
  
  x <- c('CHR', 'POS', 'HGall')
  colnames(combined_chr_specific_data) <- x
  
  combined_chr_specific_data$CHR <- i
  combined_chr_specific_data$POS <- unlist(chr_specific_map[1,])
  
  
 #trying this instead of for loop that is used in previous script:
  
  
  combined_chr_specific_data <- combined_chr_specific_data %>%
  mutate(
    HGall = rowSums(chr_specific_txt_correct_indiv == 2, na.rm = TRUE) / (rowSums(chr_specific_txt_correct_indiv == 1, na.rm = TRUE) + rowSums(chr_specific_txt_correct_indiv == 2, na.rm = TRUE))
  )

    
  #want to output these files at some point but unsure where that will be / how to do
  
#for each chromosome, make a new dataframe for just that chromosome
  chr_specific_ancestry <- combined_chr_specific_data
  
  # Append the chromosome data to the list
  all_chromosomes_data[[i]] <- combined_chr_specific_data  

#}

final_output <- bind_rows(all_chromosomes_data)

output_name <- paste0("average_ancestry_output.txt")

write.table(final_output, file = output_name, quote = FALSE, row.names = FALSE)




