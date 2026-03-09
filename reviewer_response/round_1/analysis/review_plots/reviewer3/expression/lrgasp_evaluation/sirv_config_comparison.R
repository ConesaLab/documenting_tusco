#!/usr/bin/env Rscript
###############################################################################
# SIRV Baseline Comparison and Statistical Analysis
#
# Compare TUSCO configuration FN rates against SIRV (spike-in) FN rates
# using MEDIAN FN rates across pipelines for comparison.
#
# Rationale: Median provides a robust central tendency measure that is
# less affected by outliers than mean, with statistical meaning.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# Set working directory to script location
script_dir <- getwd()

# Define paths relative to repository root
repo_root <- file.path(script_dir, "..", "..", "..", "..")

# Source figure utilities for TUSCO theme and colors
fig_utils_path <- file.path(repo_root, "scripts", "figure_utils.R")
if (file.exists(fig_utils_path)) {
  source(fig_utils_path)
} else {
  # Fallback theme if figure_utils.R not found
  theme_tusco <- function(base_size = 7, base_family = "Helvetica") {
    theme_classic(base_size = base_size, base_family = base_family) +
      theme(
        plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
        axis.title = element_text(size = base_size),
        axis.text = element_text(size = base_size - 1),
        legend.title = element_text(size = base_size, face = "bold"),
        legend.text = element_text(size = base_size - 1),
        strip.text = element_text(size = base_size, face = "bold"),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25)
      )
  }
  TUSCO_COLORS <- list(
    human_tusco = "#a8d5a0",
    mouse_tusco = "#1b9e77",
    sirvs = "#cab2d6"
  )
}

###############################################################################
# 1. Load SIRV FN Data from Figure 3a Tables
###############################################################################

cat("\n=== Loading SIRV FN data from figure3a tables ===\n")

# Human SIRV data
sirv_human_file <- file.path(repo_root, "figs/figure-03/tables/fig-3a-human.tsv")
sirv_human <- read_tsv(sirv_human_file, show_col_types = FALSE) %>%
  filter(Type == "SIRVs", final_label == "Missing", record_type == "bar_distribution") %>%
  select(pipeline, sirv_fn_count = count, sirv_fn_pct = percentage) %>%
  mutate(
    species = "hsa",
    pipeline = toupper(pipeline)
  )

# Mouse SIRV data
sirv_mouse_file <- file.path(repo_root, "figs/figure-03/tables/fig-3a-mouse.tsv")
sirv_mouse <- read_tsv(sirv_mouse_file, show_col_types = FALSE) %>%
  filter(Type == "SIRVs", final_label == "Missing", record_type == "bar_distribution") %>%
  select(pipeline, sirv_fn_count = count, sirv_fn_pct = percentage) %>%
  mutate(
    species = "mmu",
    pipeline = toupper(pipeline)
  )

# Combine SIRV data
sirv_data <- bind_rows(sirv_human, sirv_mouse)

cat("SIRV FN data loaded:\n")
print(sirv_data)

###############################################################################
# 2. Load TUSCO Config Evaluation Results
###############################################################################

cat("\n=== Loading TUSCO config evaluation results ===\n")

eval_file <- file.path(script_dir, "evaluation_results.tsv")
eval_df <- read_tsv(eval_file, show_col_types = FALSE) %>%
  mutate(
    fn_pct = 100 * n_fn / n_genes,
    config_type = case_when(
      str_detect(config_name, "strict") ~ "strict",
      str_detect(config_name, "current") ~ "current",
      str_detect(config_name, "more_lenient") ~ "more_lenient",
      str_detect(config_name, "lenient") ~ "lenient",
      TRUE ~ "unknown"
    ),
    pipeline = toupper(pipeline)
  )

###############################################################################
# 3. Join TUSCO Configs with SIRV Baseline
###############################################################################

cat("\n=== Joining TUSCO configs with SIRV baseline ===\n")

combined_df <- eval_df %>%
  left_join(sirv_data, by = c("species", "pipeline"))

###############################################################################
# 4. Calculate MEDIAN Statistics per Config
###############################################################################

cat("\n=== Calculating median statistics ===\n")

# Calculate median, Q1, Q3 for each config
config_stats <- combined_df %>%
  group_by(config_name, species, config_type) %>%
  summarise(
    n_genes = first(n_genes),
    median_tusco_fn = median(fn_pct, na.rm = TRUE),
    q1_tusco_fn = quantile(fn_pct, 0.25, na.rm = TRUE),
    q3_tusco_fn = quantile(fn_pct, 0.75, na.rm = TRUE),
    median_sirv_fn = median(sirv_fn_pct, na.rm = TRUE),
    q1_sirv_fn = quantile(sirv_fn_pct, 0.25, na.rm = TRUE),
    q3_sirv_fn = quantile(sirv_fn_pct, 0.75, na.rm = TRUE),
    n_pipelines = n(),
    .groups = "drop"
  ) %>%
  mutate(
    fn_difference = median_tusco_fn - median_sirv_fn,
    species_label = ifelse(species == "hsa", "Human", "Mouse"),
    config_label = case_when(
      config_type == "strict" ~ "Strict",
      config_type == "current" ~ "Current",
      config_type == "lenient" ~ "Lenient",
      config_type == "more_lenient" ~ "More Lenient"
    )
  )

cat("\nMedian Statistics by Config:\n")
print(config_stats %>%
        select(species_label, config_label, n_genes, median_tusco_fn, median_sirv_fn, fn_difference) %>%
        arrange(species_label, factor(config_label, levels = c("Strict", "Current", "Lenient", "More Lenient"))))

###############################################################################
# 5. Create Visualization with TUSCO Theme
###############################################################################

cat("\n=== Creating visualization ===\n")

# Order configs
config_stats$config_label <- factor(
  config_stats$config_label,
  levels = c("Strict", "Current", "Lenient", "More Lenient")
)

# Define colors using TUSCO palette
config_colors <- c(
  "Strict" = "#2166ac",
  "Current" = "#4daf4a",
  "Lenient" = "#ff7f00",
  "More Lenient" = "#e41a1c"
)

# Create bar plot with IQR error bars
p <- ggplot(config_stats, aes(x = config_label, y = median_tusco_fn, fill = config_label)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = q1_tusco_fn, ymax = q3_tusco_fn),
    width = 0.2,
    linewidth = 0.3
  ) +
  geom_hline(
    aes(yintercept = median_sirv_fn),
    linetype = "dashed",
    color = TUSCO_COLORS$sirvs,
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", median_tusco_fn), y = q3_tusco_fn + 3),
    size = 2,
    fontface = "bold"
  ) +
  geom_text(
    aes(label = paste0("n=", n_genes), y = -4),
    size = 1.8,
    color = "gray40"
  ) +
  facet_wrap(~species_label, scales = "free_y") +
  scale_fill_manual(values = config_colors) +
  labs(
    x = "TUSCO Configuration",
    y = "False Negative Rate (%)",
    title = "TUSCO Config FN Rate vs SIRV Baseline",
    subtitle = "Median across 6 pipelines; error bars = IQR; dashed line = SIRV median"
  ) +
  coord_cartesian(ylim = c(-8, max(config_stats$q3_tusco_fn) + 10)) +
  theme_tusco(base_size = 7) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 5, hjust = 0.5, color = "gray40"),
    strip.background = element_rect(fill = "grey95", linewidth = 0.25)
  )

###############################################################################
# 6. Save Output Files
###############################################################################

cat("\n=== Saving output files ===\n")

# Save PDF (Nature single column width)
pdf_file <- file.path(script_dir, "sirv_config_comparison.pdf")
ggsave(pdf_file, p, width = 3.35, height = 3.0, units = "in", device = "pdf", dpi = 300)
cat("Saved plot:", pdf_file, "\n")

# Save TSV with comparison data
tsv_file <- file.path(script_dir, "sirv_config_comparison.tsv")

output_data <- config_stats %>%
  select(
    config_name,
    species,
    species_label,
    config_type,
    config_label,
    n_genes,
    median_tusco_fn_pct = median_tusco_fn,
    q1_tusco_fn_pct = q1_tusco_fn,
    q3_tusco_fn_pct = q3_tusco_fn,
    median_sirv_fn_pct = median_sirv_fn,
    fn_difference,
    n_pipelines
  )

write_tsv(output_data, tsv_file)
cat("Saved data:", tsv_file, "\n")

###############################################################################
# 7. Summary and Interpretation
###############################################################################

cat("\n=== SUMMARY ===\n")
cat("\nMethodology: Using MEDIAN FN rate across all 6 pipelines for comparison\n")
cat("- Provides robust central tendency less affected by outliers\n")
cat("- Error bars show interquartile range (IQR: 25th-75th percentile)\n")

cat("\n=== Final Comparison Table ===\n")
final_table <- config_stats %>%
  select(species_label, config_label, n_genes,
         median_tusco_fn, median_sirv_fn, fn_difference) %>%
  mutate(across(where(is.numeric) & !matches("n_genes"), ~round(., 1))) %>%
  arrange(species_label, config_label)

print(final_table)

cat("\n=== Key Findings ===\n")
cat("\nThe 'Current' threshold balances:\n")
cat("1. Gene count: Provides sufficient genes for robust benchmarking\n")
cat("2. FN rate: Maintains acceptable detection rate across platforms\n")
cat("3. Expression: Ensures genes are detectable (medium-high expression)\n")
cat("\nLowering thresholds (Lenient/More Lenient) sharply increases FN rate,\n")
cat("compromising the reliability of the benchmark.\n")

cat("\nAnalysis complete!\n")
