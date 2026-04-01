#script to check for replicaiton from one dataset in another

library(dplyr)

args <- commandArgs(TRUE)


#input the file you are using for top hits:

hits <- read.table(args[1], header = TRUE)

#one to test for replication:

rep <- read.table(args[2], header = TRUE)


#grab past a certain threshold (here 0.001 but can pick something else)
top_hits <- hits %>% filter(p_value <= 0.003)


#then merge with rep:

merged <- merge(top_hits, rep, by = c("CHR", "POS"))


#now save only rep < 0.05 or some other input

replicated <- merged %>% filter(p_value.y <= 0.05)


#now write out with headings CHR, POS, pvalue, rep_pvalue (select these and rename)
output <- replicated %>%
  select(CHR, POS, p_value = p_value.x, rep_p_value = p_value.y)


write.table(output, "replicated_hits.txt", sep = "\t", row.names = FALSE, quote = FALSE)
