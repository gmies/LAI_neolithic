# Define file paths
individuals_file <- "/project/mathilab/gmies/neolithic_selection/allentoft_data/031825_rerun_filter_allen_LAI/GA/simplai_7/mneo_samples"
posterior_dir <- "../../../ancestry_hmm"
output_dir <- "."  # Replace with your desired output directory

# Read the list of individuals
individuals <- readLines(individuals_file)

# Loop through each chromosome
for (chr in 1:22) {
  cat("Processing chromosome", chr, "\n")
  
  # Initialize a list to store data for this chromosome
  chr_data <- list()
  
  # Loop through each individual
  for (ind in individuals) {
    # Construct the file path for the individual's .posterior file
posterior_file <- file.path(posterior_dir, paste0(ind, ".posterior.gz"))
    
    # Check if the file exists
    if (file.exists(posterior_file)) {
      # Read the .posterior file
      posterior_data <- read.table(posterior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      
      # Filter rows where chrom == chr
      chr_rows <- posterior_data[posterior_data$chrom == chr, ]
      
      # Extract the last three columns
posterior_values <- chr_rows[, 3:5]
      
      # Determine the ancestry state for each row
      ancestry_states <- apply(posterior_values, 1, function(row) {
        max_val <- max(row)

        if (max_val <= 0.9) {
          return("NA NA")
        } else if (max_val == row[1]) {
          return("2 2")
        } else if (max_val == row[2]) {
          return("1 2")
        } else {
          return("1 1")
        }
      })
      
      # Add the ancestry states to the data
      chr_data[[ind]] <- ancestry_states
    } else {
      cat("Warning: File not found:", posterior_file, "\n")
    }
  }
  
  # Combine the data for this chromosome
  chr_data_combined <- do.call(cbind, chr_data)
  
  # Write the combined data to a file
  output_file <- file.path(output_dir, paste0("output_ancestryhmm_", chr, ".txt.gz"))
  write.table(chr_data_combined, output_file, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat("Finished processing chromosome", chr, "\n")
}
