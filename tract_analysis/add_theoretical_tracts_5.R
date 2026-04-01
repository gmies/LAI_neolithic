#plot to run density distributions together with theoretical distributions:

library(tidyverse)
library(ggplot2)

# Parse command-line arguments
args <- commandArgs(TRUE)

# The first argument is the number of datasets
num_datasets <- as.numeric(args[1])

#second arg is pop name
pop <- args[2]

# The next 'num_datasets' arguments are the paths to the data files
dataset_paths <- args[3:(2 + num_datasets)]

# The following 'num_datasets' arguments are the dataset names (labels for the sources)
dataset_names <- args[(3 + num_datasets):(2 + 2 * num_datasets)]

# Decide which ancestry to filter based on population
target_ancestry <- ifelse(pop == "hg", "ancestry_2",
                   ifelse(pop == "farmer", "ancestry_1", NA))

if (is.na(target_ancestry)) {
  stop("Invalid population label. Use 'hg' or 'farmer'.")
}

# Function to read, filter, and label each dataset
read_and_filter_dataset <- function(file_path, source_label, ancestry_label) {
  data <- read.table(file_path, header = TRUE, sep = "\t") %>%
    filter(ancestry == ancestry_label) %>%
    mutate(Source = source_label)
  return(data)
}

# Read and process all datasets
dataset_list <- vector("list", num_datasets)
for (i in 1:num_datasets) {
  dataset_list[[i]] <- read_and_filter_dataset(dataset_paths[i], dataset_names[i], target_ancestry)
}

# Combine all filtered datasets
combined_data <- bind_rows(dataset_list)

# Check if data is empty (e.g., wrong ancestry present in input)
if (nrow(combined_data) == 0) {
  stop("No data found for the selected ancestry in provided files.")
}

write.table(combined_data, file = paste0(pop, "_combined_tract_data.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Log-spaced bins from 1 to 1000 cM
bin_edges <- 10^seq(log10(0.5), log10(1000), length.out = 60)

# Filter out zeros and negative (tract length must be >0 for log binning)
filtered_data <- combined_data %>%
  filter(tract_length_cm > 0)

# Assign each tract to a bin
filtered_data <- filtered_data %>%
  mutate(bin = cut(tract_length_cm, breaks = bin_edges, include.lowest = TRUE, right = FALSE))

# Create a lookup table for bin centers to avoid regex issues
bin_centers <- data.frame(
  bin = levels(cut(1, breaks = bin_edges, include.lowest = TRUE, right = FALSE)),
  bin_start = bin_edges[-length(bin_edges)],
  bin_end = bin_edges[-1]
) %>%
  mutate(bin_center = (bin_start + bin_end) / 2)


# Generate theoretical distributions for T=10, 35, 50
T_values <- c(10, 20, 30, 40, 50, 60, 100, 200)
p_farmer <- 0.6
p_hg <- 0.4

# Summarize counts per Source and bin
binned_summary <- filtered_data %>%
  group_by(Source, bin) %>%
  summarise(count = n(), .groups = 'drop') %>%
  left_join(bin_centers, by = "bin") %>%
  filter(!is.na(bin_center))

# Calculate total tracts per Source for empirical datasets
total_counts_per_source <- binned_summary %>%
  group_by(Source) %>%
  summarise(total_count = sum(count))

# Add proportion column for empirical data by dividing counts by total count per Source
binned_summary <- binned_summary %>%
  left_join(total_counts_per_source, by = "Source") %>%
  mutate(proportion = count / total_count)


# Define scaling factor based on population proportion
p_pop <- ifelse(pop == "hg", p_hg, 
                ifelse(pop == "farmer", p_farmer, NA))

if (is.na(p_pop)) {
  stop("Population must be 'hg' or 'farmer' for theoretical scaling.")
}


# Generate theoretical data for T values
theoretical_data <- data.frame()

for (T in T_values) {
#  lambda_rate <- (T / 100)*p_pop
 lambda_rate <- (T*(1-p_pop)) / 100
  
  theory_bin_data <- bin_centers %>%
    mutate(
      bin_probability = (exp(-lambda_rate * bin_start) - exp(-lambda_rate * bin_end)),
      Source = paste0("Theory T=", T)
    ) %>%
    filter(bin_probability > 1e-6)  # Filter out tiny probs for clarity
  
  theoretical_data <- rbind(theoretical_data, theory_bin_data)
}

# Combine empirical (using proportions) and theoretical data (bin_probability)
plot_data <- rbind(
  binned_summary %>% select(Source, bin_center, proportion) %>% rename(value = proportion),
  theoretical_data %>% select(Source, bin_center, bin_probability) %>% rename(value = bin_probability)
)

# Identify empirical and theoretical sources
empirical_sources <- unique(binned_summary$Source)
theoretical_sources <- paste0("Theory T=", T_values)

# Plot
pdf(paste0(pop, "_tract_length_with_theory.pdf"), width = 10, height = 6)

p <- ggplot(plot_data, aes(x = bin_center, y = value, color = Source)) +
  geom_point(data = plot_data %>% filter(!grepl("Theory", Source)),
             alpha = 0.7, size = 2) +
  geom_line(data = plot_data %>% filter(grepl("Theory", Source)),
            aes(group = Source), linewidth = 1.5, alpha = 0.9) +
  scale_x_log10(name = "Tract Length (cM)", 
                breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
                labels = c("1", "2", "5", "10", "20", "50", "100", "200", "500", "1000")) +
  scale_y_log10(name = "Proportion of tracts (log scale)") +
  scale_color_manual(values = c(
    setNames(c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")[seq_along(empirical_sources)], empirical_sources),
  "Theory T=10" = "black",
  "Theory T=20" = "gray15",
  "Theory T=30" = "gray25",
  "Theory T=40" = "gray35",
  "Theory T=50" = "gray45",
  "Theory T=60" = "gray55",
"Theory T=100" = "gray65",
"Theory T=200" = "gray75"
  )) +
  ggtitle(paste("Ancestry Tract Length Distribution -", 
                ifelse(pop == "hg", "Hunter-Gatherer", "Farmer"), "Ancestry")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(
    linetype = c(rep("blank", length(empirical_sources)), rep("solid", length(theoretical_sources))),
    shape = c(rep(16, length(empirical_sources)), rep(NA, length(theoretical_sources)))
  )))

print(p)
dev.off()


# Print summary statistics
cat("\nSummary Statistics:\n")
cat(sprintf("Total empirical tracts: %d\n", total_empirical_tracts))
for (T in T_values) {
  lambda_rate <- T / 100
  mean_length <- 1 / lambda_rate
  percentile_90 <- -log(0.1) / lambda_rate
  percentile_95 <- -log(0.05) / lambda_rate
  
  cat(sprintf("\nT = %d generations:\n", T))
  cat(sprintf("  Mean tract length: %.2f cM\n", mean_length))
  cat(sprintf("  90%% of tracts shorter than: %.2f cM\n", percentile_90))
  cat(sprintf("  95%% of tracts shorter than: %.2f cM\n", percentile_95))
}

cat(sprintf("\nPlot saved as '%s_tract_length_with_theory.pdf'\n", pop))
