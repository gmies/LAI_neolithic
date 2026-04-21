
# LOAD LIBRARIES
library(tidyverse)
library(cowplot)


# COMMAND LINE ARGUMENTS

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <sex_file> <method_file1> <method_file2> ...")
}

sex_file <- args[1]
method_files <- args[-1]

 
# READ SEX FILE
 
sex_df <- read.table(sex_file, header = FALSE, stringsAsFactors = FALSE)
colnames(sex_df) <- c("sample", "ploidy", "ancestry_ref")

sex_df <- sex_df %>%
  mutate(
    sex = factor(ploidy, levels = c(1, 2), labels = c("Male", "Female"))
  )

 
# READ ALL METHOD FILES
 
read_method <- function(file) {
  method_name <- file %>%
    basename() %>%
    gsub("^chr_x_", "", .) %>%
    gsub("\\.txt$", "", .)

  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  colnames(df) <- c("sample", "prop")

  df %>%
    mutate(method = method_name)
}

method_df <- map_dfr(method_files, read_method)

 
# ADD QPADM
 
qpadm_df <- read.table("/project/mathilab/gmies/neolithic_selection/X_chromosome/LAI/qpadm/global_ancestry.txt",
                       header = TRUE, stringsAsFactors = FALSE)
colnames(qpadm_df) <- c("sample", "prop")

qpadm_df <- qpadm_df %>%
  mutate(method = "X_qpadm")

method_df <- bind_rows(method_df, qpadm_df)

 
# MERGE WITH SEX
 
combined_df <- method_df %>%
  inner_join(sex_df %>% select(sample, sex), by = "sample")

 
# T-TESTS PER METHOD
 
ttest_df <- combined_df %>%
  group_by(method) %>%
  summarise(
    p_value = tryCatch(
      t.test(prop ~ sex)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = ifelse(is.na(p_value),
                     "p = NA",
                     paste0("p = ", signif(p_value, 3)))
  )

combined_df <- combined_df %>%
  left_join(ttest_df, by = "method")

 
# PLOTTING FUNCTION
 
plot_methods <- function(df, output_file) {

  p <- ggplot(df, aes(x = sex, y = prop, color = sex)) +
    geom_jitter(width = 0.2, size = 1.8, alpha = 0.8) +
    geom_boxplot(outlier.shape = NA, alpha = 0.3) +
    facet_wrap(~ method, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c("Male" = "steelblue",
                                  "Female" = "orange")) +
    geom_text(
      data = distinct(df, method, p_label),
      aes(x = 1.5, y = Inf, label = p_label),
      inherit.aes = FALSE,
      vjust = 1.2,
      size = 3.2
    ) +
    labs(
      x = "Sex",
      y = "X ancestry proportion",
      color = "Sex"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(output_file, p, width = 9, height = 7, dpi = 300)
}

 
# OUTPUT 1: WITH FLARE
 
plot_methods(
  combined_df,
  output_file = "X_methods_with_flare.pdf"
)

 
# OUTPUT 2: WITHOUT FLARE
 
plot_methods(
  combined_df %>% filter(method != "flare"),
  output_file = "X_methods_no_flare.pdf"
)

message("Plots saved:")
message("  - X_methods_with_flare.pdf")
message("  - X_methods_no_flare.pdf")
