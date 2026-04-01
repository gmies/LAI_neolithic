# combine_indivs.R
#-------------------------------------------------------
# Purpose: Combine per-chromosome allele count files 
#          for each admixed individual into a single file
# Input: Number of admixed individuals (from command-line)
#        Files named as {chr}_{indiv}_mneo.frq.counts in working directory
# Output: One combined file per individual: {indiv}mneo.counts
#-------------------------------------------------------

args <- commandArgs(TRUE)
num_indivs <- as.numeric(args[1])  # Number of admixed individuals

# Loop over each individual
for (indiv in 1:num_indivs) {
    
    indiv_concat <- NULL  # Initialize concatenated data for this individual
    
    # Loop over chromosomes 1-22
    for (chr in 1:22) {
        file_name <- paste0(chr, "_", indiv, "_mneo.frq.counts")
        chr_data <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE)
        
        # Concatenate chromosome data
        if (is.null(indiv_concat)) {
            indiv_concat <- chr_data
        } else {
            indiv_concat <- rbind(indiv_concat, chr_data)
        }
    }
    
    # Save concatenated data for this individual
    write.table(
        indiv_concat, 
        paste0(indiv, "mneo.counts"), 
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
    )
}
