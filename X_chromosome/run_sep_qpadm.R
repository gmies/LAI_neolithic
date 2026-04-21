
library(admixtools)
library(tidyverse)

args <- commandArgs(TRUE)

# Your full population set
left <- c("NEO", "HG")

right <- c(
 "Karitiana.DG", "Mbuti.DG", "Papuan.DG", "Han.DG",
  "Russia_Ust_Ishim_HG_published.DG", "Belgium_UP_GoyetQ116_1_published_all",
  "Russia_MA1_HG.SG", "Russia_Kostenki14.SG", "Italy_North_Villabruna_HG",
  "Czech_Vestonice16", "Russia_Kostenki14", "Spain_ElMiron", "Russia_AfontovaGora3"
)

#right <- c("Papuan.DG", "Han.DG")

#diff for replication 
ind_df <- read_table(args[1], col_names = FALSE)
#ind_df <- read_table("/project/mathilab/gmies/neolithic_selection/qpadm/0100825_allentoft_qpadm/output_data/qpadm_input.ind", col_names = FALSE)
colnames(ind_df) <- c("name", "sex", "group")
targets <- ind_df %>% filter(sex == "U") %>% pull(name)
length(targets)

# Loop over targets and collect results
results_list <- map(targets, function(target) {
  message("Running qpAdm for: ", target)
  pops_for_target <- unique(c(target, left, right))
  
  # Recompute f2_blocks per target to minimize RAM usage
  f2_blocks <- f2_from_geno("qpadm_input", pops = pops_for_target, maxmiss = 1, minmaf = 0, maxmaf = 1, adjust_pseudohaploid = TRUE, poly_only = NULL, auto_only = FALSE)
  
  tryCatch({
    res <- qpadm(f2_blocks, left = left, right = right, target = target)
    
    res$weights %>%
      mutate(individual = target) %>%
      select(individual, left, weight)
    
  }, error = function(e) {
    message("Failed on: ", target, " — ", e$message)
    tibble(individual = target, left = NA, weight = NA)
  })
})

# Combine all into a single tibble
results_df <- bind_rows(results_list)

# Optionally, pivot to wide format
final_df <- results_df %>%
  pivot_wider(names_from = left, values_from = weight)

# Save to CSV
write_csv(final_df, "qpadm_individual_ancestry_proportions.csv")

hg_ancestry <- results_df %>%
  filter(left == "HG") %>%
  select(individual, global_ancestry = weight)
  
  write.table(hg_ancestry, "global_ancestry.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)


