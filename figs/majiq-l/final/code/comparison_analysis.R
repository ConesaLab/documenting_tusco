#!/usr/bin/env Rscript
# ==============================================================================
# TUSCO vs MAJIQ-L Benchmarking Comparison Analysis
# ==============================================================================
# This script compares benchmarking metrics between TUSCO and MAJIQ-L approaches
# across 6 long-read sequencing pipelines.
#
# Metrics computed:
# - Sensitivity = TP / (TP + FN)
# - FN Rate = FN / (TP + FN)
# - FDR = 100 × (N - TP) / N  (N = observed transcripts/junctions)
#
# MAJIQ-L calculations include BOTH known and novel junctions.
#
# Output:
# - Three-panel figure comparing Sensitivity, FN rate, and FDR
# - Comparison metrics table
# ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Setup paths
# ------------------------------------------------------------------------------
# Get script directory from command line args or use current directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()

# Repository root
repo_root <- normalizePath(file.path(script_dir, "../../../.."), mustWork = FALSE)

# Output directories
output_base <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
plot_dir <- file.path(output_base, "plots")
table_dir <- file.path(output_base, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Input files
tusco_file <- file.path(repo_root, "figs/figure-03/tables/figure3b-human.tsv")
majiq_file <- file.path(repo_root, "figs/majiq-l/tables/junction_category_totals.tsv")

cat("Script directory:", script_dir, "\n")
cat("Repository root:", repo_root, "\n")
cat("TUSCO data file:", tusco_file, "\n")
cat("MAJIQ-L data file:", majiq_file, "\n")

# ------------------------------------------------------------------------------
# Load and process TUSCO data
# ------------------------------------------------------------------------------
cat("\n=== Loading TUSCO data ===\n")

tusco_raw <- read_tsv(tusco_file, show_col_types = FALSE)

# Filter for TUSCO type and raw records only
# Get both counts and percentages
tusco_data <- tusco_raw %>%
  filter(Type == "TUSCO", record_type == "raw") %>%
  select(big_category, count, percentage, pipeline) %>%
  pivot_wider(
    names_from = big_category,
    values_from = c(count, percentage),
    names_glue = "{big_category}_{.value}"
  ) %>%
  mutate(
    # N = observed transcripts (TP + PTP + FP, excluding FN which are not observed)
    N_observed = TP_count + PTP_count + FP_count,
    # Compute sensitivity: TP / (TP + FN)
    sensitivity = TP_percentage / (TP_percentage + FN_percentage) * 100,
    fn_rate = FN_percentage,
    # FDR = 100 × (N - TP) / N
    fdr = 100 * (N_observed - TP_count) / N_observed,
    method = "TUSCO"
  ) %>%
  select(pipeline, sensitivity, fn_rate, fdr, method,
         TP_count, PTP_count, FP_count, FN_count = FN_percentage, N_observed)

# Create pipeline labels
pipeline_labels <- c(
  "WTC11_cdna_ont" = "cDNA-ONT",
  "WTC11_cdna_ont_ls" = "cDNA-ONT (LS)",
  "WTC11_cdna_pacbio" = "cDNA-PacBio",
  "WTC11_cdna_pacbio_ls" = "cDNA-PacBio (LS)",
  "WTC11_drna_ont" = "dRNA-ONT",
  "WTC11_drna_ont_ls" = "dRNA-ONT (LS)"
)

tusco_data <- tusco_data %>%
  mutate(pipeline_label = pipeline_labels[pipeline])

cat("TUSCO data loaded:\n")
print(tusco_data %>% select(pipeline_label, TP_count, N_observed, sensitivity, fn_rate, fdr))

# ------------------------------------------------------------------------------
# Load and process MAJIQ-L data
# ------------------------------------------------------------------------------
cat("\n=== Loading MAJIQ-L data ===\n")

majiq_raw <- read_tsv(majiq_file, show_col_types = FALSE)

# Aggregate across all gene categories for each pipeline
# Include BOTH known and novel junctions
majiq_data <- majiq_raw %>%
  group_by(pipeline, pipeline_label) %>%
  summarise(
    # Known junctions (in annotation)
    TP_known = sum(All, na.rm = TRUE),
    FP_known = sum(LR_Annotation, na.rm = TRUE),  # LR found, no SR support
    FN_known = sum(MAJIQ_Annotation, na.rm = TRUE),  # SR found, LR missed
    # Novel junctions (not in annotation)
    TP_novel = sum(Both_denovo, na.rm = TRUE),
    FP_novel = sum(LR_denovo, na.rm = TRUE),
    FN_novel = sum(MAJIQ_denovo, na.rm = TRUE),
    # Total LR junctions (N for FDR)
    N_observed = sum(total_LR_junctions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Combined TP and FN (both known and novel)
    TP_total = TP_known + TP_novel,
    FN_total = FN_known + FN_novel,
    # Compute sensitivity for ALL junctions (known + novel): TP / (TP + FN)
    sensitivity = TP_total / (TP_total + FN_total) * 100,
    fn_rate = FN_total / (TP_total + FN_total) * 100,
    # FDR = 100 × (N - TP) / N
    fdr = 100 * (N_observed - TP_total) / N_observed,
    method = "MAJIQ-L"
  )

cat("MAJIQ-L data loaded (including both known and novel):\n")
print(majiq_data %>% select(pipeline_label, TP_total, FN_total, N_observed, sensitivity, fn_rate, fdr))

# ------------------------------------------------------------------------------
# Create MAJIQ-L benchmarking table (for LaTeX document)
# ------------------------------------------------------------------------------
cat("\n=== Creating MAJIQ-L benchmarking table ===\n")

majiq_table <- majiq_data %>%
  select(pipeline_label, TP_known, FP_known, FN_known, TP_novel, FP_novel, FN_novel) %>%
  arrange(factor(pipeline_label, levels = c(
    "cDNA-ONT", "cDNA-ONT (LS)",
    "cDNA-PacBio", "cDNA-PacBio (LS)",
    "dRNA-ONT", "dRNA-ONT (LS)"
  )))

write_tsv(majiq_table, file.path(table_dir, "majiq_benchmarking.tsv"))
cat("Saved MAJIQ-L benchmarking table\n")
print(majiq_table)

# ------------------------------------------------------------------------------
# Combine data for comparison
# ------------------------------------------------------------------------------
cat("\n=== Combining data for comparison ===\n")

# Prepare TUSCO for joining
tusco_for_join <- tusco_data %>%
  select(pipeline, pipeline_label, sensitivity, fn_rate, fdr, method)

# Prepare MAJIQ-L for joining
majiq_for_join <- majiq_data %>%
  select(pipeline, pipeline_label, sensitivity, fn_rate, fdr, method)

# Combine
comparison_data <- bind_rows(tusco_for_join, majiq_for_join)

# Create comparison table (wide format)
comparison_wide <- comparison_data %>%
  select(pipeline_label, method, sensitivity, fn_rate, fdr) %>%
  pivot_wider(
    names_from = method,
    values_from = c(sensitivity, fn_rate, fdr),
    names_glue = "{.value}_{method}"
  ) %>%
  arrange(factor(pipeline_label, levels = c(
    "cDNA-PacBio", "cDNA-PacBio (LS)",
    "dRNA-ONT", "dRNA-ONT (LS)",
    "cDNA-ONT", "cDNA-ONT (LS)"
  )))

write_tsv(comparison_wide, file.path(table_dir, "comparison_metrics.tsv"))
cat("Saved comparison metrics table\n")
print(comparison_wide)

# ------------------------------------------------------------------------------
# Define theme and colors
# ------------------------------------------------------------------------------
theme_comparison <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10)
    )
}

METHOD_COLORS <- c(
  "TUSCO" = "#a8d5a0",    # Light green (project standard)
  "MAJIQ-L" = "#B2182B"   # Red
)

# Order pipelines by TUSCO sensitivity (best to worst)
pipeline_order <- comparison_wide %>%
  arrange(desc(sensitivity_TUSCO)) %>%
  pull(pipeline_label)

comparison_data <- comparison_data %>%
  mutate(pipeline_label = factor(pipeline_label, levels = pipeline_order))

# ------------------------------------------------------------------------------
# Create figure: Three-panel dot plot with connecting lines
# ------------------------------------------------------------------------------
cat("\n=== Creating comparison figure ===\n")

# Prepare data for plotting with lines
plot_data_wide <- comparison_data %>%
  select(pipeline_label, method, sensitivity, fn_rate, fdr) %>%
  pivot_wider(
    names_from = method,
    values_from = c(sensitivity, fn_rate, fdr)
  )

# Panel A: Sensitivity
p_sensitivity <- ggplot() +
  # Lines connecting TUSCO and MAJIQ-L for each pipeline
  geom_segment(
    data = plot_data_wide,
    aes(x = pipeline_label, xend = pipeline_label,
        y = sensitivity_TUSCO, yend = `sensitivity_MAJIQ-L`),
    color = "gray60", linewidth = 0.8, alpha = 0.6
  ) +
  # Points for each method
  geom_point(
    data = comparison_data,
    aes(x = pipeline_label, y = sensitivity, color = method, shape = method),
    size = 4, alpha = 0.9
  ) +
  scale_color_manual(values = METHOD_COLORS, name = "Method") +
  scale_shape_manual(values = c("TUSCO" = 16, "MAJIQ-L" = 17), name = "Method") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Sensitivity",
    x = NULL,
    y = "Sensitivity (%)"
  ) +
  theme_comparison() +
  theme(legend.position = "none")

# Panel B: FN Rate
p_fn_rate <- ggplot() +
  # Lines connecting TUSCO and MAJIQ-L for each pipeline
  geom_segment(
    data = plot_data_wide,
    aes(x = pipeline_label, xend = pipeline_label,
        y = fn_rate_TUSCO, yend = `fn_rate_MAJIQ-L`),
    color = "gray60", linewidth = 0.8, alpha = 0.6
  ) +
  # Points for each method
  geom_point(
    data = comparison_data,
    aes(x = pipeline_label, y = fn_rate, color = method, shape = method),
    size = 4, alpha = 0.9
  ) +
  scale_color_manual(values = METHOD_COLORS, name = "Method") +
  scale_shape_manual(values = c("TUSCO" = 16, "MAJIQ-L" = 17), name = "Method") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "False Negative Rate",
    x = NULL,
    y = "FN Rate (%)"
  ) +
  theme_comparison() +
  theme(legend.position = "none")

# Panel C: FDR
p_fdr <- ggplot() +
  # Lines connecting TUSCO and MAJIQ-L for each pipeline
  geom_segment(
    data = plot_data_wide,
    aes(x = pipeline_label, xend = pipeline_label,
        y = fdr_TUSCO, yend = `fdr_MAJIQ-L`),
    color = "gray60", linewidth = 0.8, alpha = 0.6
  ) +
  # Points for each method
  geom_point(
    data = comparison_data,
    aes(x = pipeline_label, y = fdr, color = method, shape = method),
    size = 4, alpha = 0.9
  ) +
  scale_color_manual(values = METHOD_COLORS, name = "Method") +
  scale_shape_manual(values = c("TUSCO" = 16, "MAJIQ-L" = 17), name = "Method") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "False Discovery Rate",
    x = NULL,
    y = "FDR (%)"
  ) +
  theme_comparison() +
  theme(legend.position = "none")

# Combine panels with shared legend
p_combined <- (p_sensitivity | p_fn_rate | p_fdr) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.margin = margin(5, 5, 5, 5)
    )
  ) &
  theme(legend.position = "bottom")

# Save figure
ggsave(
  file.path(plot_dir, "fig-tusco-vs-majiql.pdf"),
  p_combined,
  width = 12, height = 5,
  device = "pdf"
)
cat("Saved figure:", file.path(plot_dir, "fig-tusco-vs-majiql.pdf"), "\n")

# Also save as PNG for quick viewing
ggsave(
  file.path(plot_dir, "fig-tusco-vs-majiql.png"),
  p_combined,
  width = 12, height = 5,
  dpi = 300
)
cat("Saved PNG:", file.path(plot_dir, "fig-tusco-vs-majiql.png"), "\n")

# ------------------------------------------------------------------------------
# Print summary
# ------------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("\nPipeline Ranking (by TUSCO Sensitivity):\n")
for (i in seq_along(pipeline_order)) {
  cat(sprintf("  %d. %s\n", i, pipeline_order[i]))
}

cat("\nKey interpretation:\n")
cat("  - Only cDNA-PacBio ranks best in both methods\n")
cat("  - TUSCO reflects RNA integrity (sample quality)\n")
cat("  - MAJIQ-L reflects sequencing depth\n")
cat("  - dRNA-ONT: High TUSCO Sn (good RNA) but low MAJIQ-L Sn (low depth)\n")
cat("  - cDNA-ONT: Low TUSCO Sn (many PTP = poor TSS/TTS)\n")

cat("\nKey Formulas:\n")
cat("  Sensitivity = TP / (TP + FN)\n")
cat("  FN Rate = FN / (TP + FN)\n")
cat("  FDR = 100 × (N - TP) / N, where N = observed transcripts/junctions\n")

cat("\n=== Analysis complete ===\n")
