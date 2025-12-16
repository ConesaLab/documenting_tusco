#!/usr/bin/env Rscript

###############################################################
# Supplementary Figure - Sensitivity Analysis
# Four-panel figure showing parameter threshold sensitivity
#
# Panel A: Gene Set Size vs Threshold
# Panel B: FP and FN vs Threshold (Key panel)
# Panel C: FP Breakdown by Category
# Panel D: PTP Count vs TSS Window
#
# Usage: Rscript sensitivity_analysis.R [output_dir] [width] [height]
###############################################################

# =============================================================================
# 1. SETUP AND LIBRARY LOADING
# =============================================================================

# Source utilities library
if (file.exists("../../../scripts/figure_utils.R")) {
  source("../../../scripts/figure_utils.R")
} else {
  stop("Cannot find figure_utils.R - please run from the code directory")
}

# Build unified context
ctx <- figure_context(defaults = list(
  out_dir = "..",
  width = 7.09,  # Nature double column width
  height = 8.0   # Tall figure for 4 panels
))
params <- ctx$params

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(cowplot)
  library(scales)
})

# Suppress warnings for global variables
utils::globalVariables(c(
  "config_name", "species", "sweep_type", "n_genes", "n_fp_total", "n_fn",
  "n_ptp", "n_tp", "threshold_value", "threshold_label", "metric", "value",
  "fp_category", "count", "threshold_order"
))

# =============================================================================
# 2. OUTPUT DIRECTORIES
# =============================================================================

plot_dir <- ctx$plot_dir
table_dir <- ctx$table_dir

message("Output directories:")
message("  Plots: ", plot_dir)
message("  Tables: ", table_dir)

# =============================================================================
# 3. DATA LOADING
# =============================================================================

# Find the evaluation results file
data_candidates <- c(
  file.path(ctx$repo_root, "reviewer_response", "sensitivity_runs", "summary", "evaluation_results.tsv"),
  file.path(ctx$figure_dir, "data", "evaluation_results.tsv"),
  file.path(ctx$repo_root, "data", "processed", "sensitivity", "evaluation_results.tsv")
)

eval_file <- NULL
for (f in data_candidates) {
  if (file.exists(f)) {
    eval_file <- f
    break
  }
}

if (is.null(eval_file)) {
  stop("Cannot find evaluation_results.tsv. Please run evaluate_sensitivity_runs.py first.\n",
       "Searched in:\n  ", paste(data_candidates, collapse = "\n  "))
}

message("Loading evaluation data from: ", eval_file)
eval_data <- read_tsv(eval_file, show_col_types = FALSE)
message("  Loaded ", nrow(eval_data), " rows")

# =============================================================================
# 4. DATA PROCESSING
# =============================================================================

# Parse threshold values from config names
parse_threshold <- function(config_name, sweep_type) {
  # Extract numeric threshold from config name patterns
  # e.g., "hsa_splice_0.01" -> 0.01, "hsa_tss_300" -> 300

  if (sweep_type == "splice") {
    # Pattern: species_splice_value (e.g., hsa_splice_strict, hsa_splice_lenient1)
    if (grepl("strict", config_name)) return(0.001)
    if (grepl("default", config_name)) {
      if (grepl("^hsa", config_name)) return(0.01)
      if (grepl("^mmu", config_name)) return(0.05)
    }
    if (grepl("lenient1", config_name)) {
      if (grepl("^hsa", config_name)) return(0.05)
      if (grepl("^mmu", config_name)) return(0.15)
    }
    if (grepl("lenient2", config_name)) {
      if (grepl("^hsa", config_name)) return(0.10)
      if (grepl("^mmu", config_name)) return(0.30)
    }
    if (grepl("lenient3", config_name)) {
      if (grepl("^hsa", config_name)) return(0.25)
      if (grepl("^mmu", config_name)) return(0.50)
    }
  }

  if (sweep_type == "tss") {
    if (grepl("off", config_name)) return(0)
    if (grepl("lenient", config_name)) return(50)
    if (grepl("default", config_name)) return(300)
    if (grepl("strict", config_name)) return(1000)
  }

  if (sweep_type == "expression") {
    if (grepl("lenient2", config_name)) return(0.1)
    if (grepl("lenient1", config_name)) return(0.5)
    if (grepl("default", config_name)) return(1.0)
    if (grepl("strict", config_name)) return(5.0)
  }

  return(NA_real_)
}

# Get threshold label for display
get_threshold_label <- function(config_name, sweep_type) {
  if (sweep_type == "splice") {
    if (grepl("strict", config_name)) return("Strict")
    if (grepl("default", config_name)) return("Default")
    if (grepl("lenient1", config_name)) return("Lenient")
    if (grepl("lenient2", config_name)) return("Very Lenient")
    if (grepl("lenient3", config_name)) return("Extreme")
  }

  if (sweep_type == "tss") {
    if (grepl("off", config_name)) return("Off (0 bp)")
    if (grepl("lenient", config_name)) return("50 bp")
    if (grepl("default", config_name)) return("300 bp")
    if (grepl("strict", config_name)) return("1000 bp")
  }

  if (sweep_type == "expression") {
    if (grepl("lenient2", config_name)) return("0.1 RPKM")
    if (grepl("lenient1", config_name)) return("0.5 RPKM")
    if (grepl("default", config_name)) return("1.0 RPKM")
    if (grepl("strict", config_name)) return("5.0 RPKM")
  }

  return(config_name)
}

# Get threshold order for plotting
get_threshold_order <- function(config_name, sweep_type) {
  # Returns order for x-axis (lenient to strict, left to right)
  if (sweep_type == "splice") {
    if (grepl("strict", config_name)) return(5)
    if (grepl("default", config_name)) return(4)
    if (grepl("lenient1", config_name)) return(3)
    if (grepl("lenient2", config_name)) return(2)
    if (grepl("lenient3", config_name)) return(1)
  }

  if (sweep_type == "tss") {
    # TSS: larger window = stricter = right side
    if (grepl("off", config_name)) return(1)
    if (grepl("lenient", config_name)) return(2)
    if (grepl("default", config_name)) return(3)
    if (grepl("strict", config_name)) return(4)
  }

  if (sweep_type == "expression") {
    # Expression: higher threshold = stricter = right side
    if (grepl("lenient2", config_name)) return(1)
    if (grepl("lenient1", config_name)) return(2)
    if (grepl("default", config_name)) return(3)
    if (grepl("strict", config_name)) return(4)
  }

  return(NA_integer_)
}

# Add parsed columns to data
eval_data <- eval_data %>%
  rowwise() %>%
  mutate(
    threshold_value = parse_threshold(config_name, sweep_type),
    threshold_label = get_threshold_label(config_name, sweep_type),
    threshold_order = get_threshold_order(config_name, sweep_type)
  ) %>%
  ungroup()

# Aggregate across pipelines (average per config)
summary_data <- eval_data %>%
  group_by(config_name, species, sweep_type, threshold_value, threshold_label, threshold_order) %>%
  summarise(
    n_genes = mean(n_genes, na.rm = TRUE),
    n_tp = mean(n_tp, na.rm = TRUE),
    n_ptp = mean(n_ptp, na.rm = TRUE),
    n_fn = mean(n_fn, na.rm = TRUE),
    n_fp_total = mean(n_fp_total, na.rm = TRUE),
    n_fp_nic = mean(n_fp_nic, na.rm = TRUE),
    n_fp_nnc = mean(n_fp_nnc, na.rm = TRUE),
    n_fp_genic = mean(n_fp_genic, na.rm = TRUE),
    n_fp_antisense = mean(n_fp_antisense, na.rm = TRUE),
    n_fp_intergenic = mean(n_fp_intergenic, na.rm = TRUE),
    sensitivity = mean(sensitivity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    .groups = "drop"
  )

# Species labels
species_labels <- c("hsa" = "Human", "mmu" = "Mouse")
sweep_labels <- c("splice" = "Splice Threshold", "tss" = "TSS Window", "expression" = "Expression Threshold")

# =============================================================================
# 5. COLOR PALETTE
# =============================================================================

# Use TUSCO colors from figure_utils
sweep_colors <- c(
  "splice" = "#2ca02c",      # Green
  "tss" = "#1f77b4",         # Blue
  "expression" = "#d62728"   # Red
)

species_colors <- c(
  "hsa" = TUSCO_COLORS$human_tusco,
  "mmu" = TUSCO_COLORS$mouse_tusco
)

metric_colors <- c(
  "FP" = "#d62728",
  "FN" = "#1f77b4",
  "TP" = "#2ca02c",
  "PTP" = "#ff7f0e"
)

fp_category_colors <- c(
  "NIC" = "#e41a1c",
  "NNC" = "#377eb8",
  "Genic" = "#4daf4a",
  "Antisense" = "#984ea3",
  "Intergenic" = "#ff7f00"
)

# =============================================================================
# 6. PANEL A: Gene Set Size vs Threshold
# =============================================================================

create_panel_a <- function(data, sweep = "splice") {
  sweep_data <- data %>%
    filter(sweep_type == sweep) %>%
    arrange(threshold_order)

  if (nrow(sweep_data) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No data for", sweep)))
  }

  # Create ordered factor for x-axis
  sweep_data <- sweep_data %>%
    mutate(threshold_label = factor(threshold_label, levels = unique(threshold_label[order(threshold_order)])))

  ggplot(sweep_data, aes(x = threshold_label, y = n_genes, color = species, group = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = species_colors, labels = species_labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = sweep_labels[sweep],
      x = NULL,
      y = "Number of Genes",
      color = "Species"
    ) +
    theme_tusco() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      legend.position = "bottom",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold")
    )
}

# =============================================================================
# 7. PANEL B: FP and FN vs Threshold (KEY PANEL)
# =============================================================================

create_panel_b <- function(data, sweep = "splice") {
  sweep_data <- data %>%
    filter(sweep_type == sweep) %>%
    arrange(threshold_order) %>%
    select(config_name, species, threshold_label, threshold_order, n_fp_total, n_fn) %>%
    pivot_longer(
      cols = c(n_fp_total, n_fn),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(metric, "n_fp_total" = "FP", "n_fn" = "FN"),
      threshold_label = factor(threshold_label, levels = unique(threshold_label[order(threshold_order)]))
    )

  if (nrow(sweep_data) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No data for", sweep)))
  }

  ggplot(sweep_data, aes(x = threshold_label, y = value, color = metric, group = interaction(species, metric))) +
    geom_line(aes(linetype = species), linewidth = 0.8) +
    geom_point(aes(shape = species), size = 2) +
    scale_color_manual(values = metric_colors[c("FP", "FN")]) +
    scale_linetype_manual(values = c("hsa" = "solid", "mmu" = "dashed"), labels = species_labels) +
    scale_shape_manual(values = c("hsa" = 16, "mmu" = 17), labels = species_labels) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~species, scales = "free_y", labeller = labeller(species = species_labels)) +
    labs(
      title = paste(sweep_labels[sweep], "- FP and FN"),
      subtitle = "Both FP and FN decrease as threshold becomes stricter (left to right)",
      x = "Threshold (Lenient \u2192 Strict)",
      y = "Count",
      color = "Metric",
      linetype = "Species",
      shape = "Species"
    ) +
    theme_tusco() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6, color = "gray40"),
      strip.text = element_text(size = 7, face = "bold")
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = "none",
      shape = "none"
    )
}

# =============================================================================
# 8. PANEL C: FP Breakdown by Category
# =============================================================================

create_panel_c <- function(data, sweep = "splice") {
  sweep_data <- data %>%
    filter(sweep_type == sweep) %>%
    arrange(threshold_order) %>%
    select(config_name, species, threshold_label, threshold_order,
           n_fp_nic, n_fp_nnc, n_fp_genic, n_fp_antisense, n_fp_intergenic) %>%
    pivot_longer(
      cols = starts_with("n_fp_"),
      names_to = "fp_category",
      values_to = "count"
    ) %>%
    mutate(
      fp_category = recode(fp_category,
        "n_fp_nic" = "NIC",
        "n_fp_nnc" = "NNC",
        "n_fp_genic" = "Genic",
        "n_fp_antisense" = "Antisense",
        "n_fp_intergenic" = "Intergenic"
      ),
      fp_category = factor(fp_category, levels = c("NIC", "NNC", "Genic", "Antisense", "Intergenic")),
      threshold_label = factor(threshold_label, levels = unique(threshold_label[order(threshold_order)]))
    )

  if (nrow(sweep_data) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No data for", sweep)))
  }

  ggplot(sweep_data, aes(x = threshold_label, y = count, fill = fp_category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_manual(values = fp_category_colors) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~species, scales = "free_y", labeller = labeller(species = species_labels)) +
    labs(
      title = paste(sweep_labels[sweep], "- FP by Category"),
      x = "Threshold (Lenient \u2192 Strict)",
      y = "FP Count",
      fill = "FP Category"
    ) +
    theme_tusco() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      legend.position = "bottom",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      strip.text = element_text(size = 7, face = "bold")
    ) +
    guides(fill = guide_legend(nrow = 1))
}

# =============================================================================
# 9. PANEL D: PTP Count vs TSS Window
# =============================================================================

create_panel_d <- function(data) {
  tss_data <- data %>%
    filter(sweep_type == "tss") %>%
    arrange(threshold_order) %>%
    mutate(threshold_label = factor(threshold_label, levels = unique(threshold_label[order(threshold_order)])))

  if (nrow(tss_data) == 0) {
    return(ggplot() + theme_void() + labs(title = "No TSS data available"))
  }

  ggplot(tss_data, aes(x = threshold_label, y = n_ptp, color = species, group = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = species_colors, labels = species_labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "PTP Count vs TSS Window Size",
      subtitle = "Stricter TSS windows (larger) convert more PTPs to FPs",
      x = "TSS Window (Lenient \u2192 Strict)",
      y = "PTP Count",
      color = "Species"
    ) +
    theme_tusco() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      legend.position = "bottom",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6, color = "gray40")
    )
}

# =============================================================================
# 10. CREATE COMBINED FIGURE
# =============================================================================

message("\nCreating panels...")

# Check available sweep types
available_sweeps <- unique(summary_data$sweep_type)
message("  Available sweeps: ", paste(available_sweeps, collapse = ", "))

# Create panels based on available data
panels <- list()

# Panel A: Gene set sizes for all sweeps
if ("splice" %in% available_sweeps) {
  panels$a1 <- create_panel_a(summary_data, "splice")
} else {
  panels$a1 <- ggplot() + theme_void() + labs(title = "Splice sweep not available")
}

if ("expression" %in% available_sweeps) {
  panels$a2 <- create_panel_a(summary_data, "expression")
} else {
  panels$a2 <- ggplot() + theme_void() + labs(title = "Expression sweep not available")
}

# Panel B: FP/FN for splice sweep (key panel)
if ("splice" %in% available_sweeps) {
  panels$b <- create_panel_b(summary_data, "splice")
} else {
  panels$b <- ggplot() + theme_void() + labs(title = "FP/FN - No splice data")
}

# Panel C: FP breakdown
if ("splice" %in% available_sweeps) {
  panels$c <- create_panel_c(summary_data, "splice")
} else {
  panels$c <- ggplot() + theme_void() + labs(title = "FP breakdown - No data")
}

# Panel D: TSS PTP
panels$d <- create_panel_d(summary_data)

# Combine panels
message("  Combining panels...")

# Top row: Gene set size panels
row_a <- plot_grid(
  panels$a1, panels$a2,
  ncol = 2, align = "hv",
  labels = c("A", "B"),
  label_size = 10, label_fontface = "bold"
)

# Middle row: FP/FN panel (full width)
row_b <- plot_grid(
  panels$b,
  ncol = 1,
  labels = c("C"),
  label_size = 10, label_fontface = "bold"
)

# Bottom row: FP breakdown and PTP
row_c <- plot_grid(
  panels$c, panels$d,
  ncol = 2, align = "hv",
  labels = c("D", "E"),
  label_size = 10, label_fontface = "bold"
)

# Final combined figure
combined <- plot_grid(
  row_a, row_b, row_c,
  ncol = 1,
  rel_heights = c(0.8, 1.0, 0.9)
)

# =============================================================================
# 11. SAVE OUTPUTS
# =============================================================================

message("\nSaving outputs...")

# Save plot
out_pdf <- file.path(plot_dir, "fig-sensitivity.pdf")
ggsave(out_pdf, combined, width = params$width, height = params$height, units = "in", device = "pdf")
message("  Saved plot: ", out_pdf)

# Save PNG version
out_png <- file.path(plot_dir, "fig-sensitivity.png")
ggsave(out_png, combined, width = params$width, height = params$height, units = "in", dpi = 300)
message("  Saved plot: ", out_png)

# Save summary table
summary_out <- file.path(table_dir, "sensitivity_summary.tsv")
write_tsv(summary_data, summary_out)
message("  Saved table: ", summary_out)

# Save detailed evaluation data
eval_out <- file.path(table_dir, "sensitivity_evaluation.tsv")
write_tsv(eval_data, eval_out)
message("  Saved table: ", eval_out)

# =============================================================================
# 12. GENERATE KEY FINDINGS SUMMARY
# =============================================================================

message("\n=== KEY FINDINGS ===")

# For each sweep type, show the range of genes and FP/FN
for (sweep in c("splice", "tss", "expression")) {
  sweep_subset <- summary_data %>% filter(sweep_type == sweep)
  if (nrow(sweep_subset) == 0) next

  message(sprintf("\n%s Sweep:", sweep_labels[sweep]))

  for (sp in c("hsa", "mmu")) {
    sp_data <- sweep_subset %>% filter(species == sp)
    if (nrow(sp_data) == 0) next

    min_genes <- min(sp_data$n_genes, na.rm = TRUE)
    max_genes <- max(sp_data$n_genes, na.rm = TRUE)
    min_fp <- min(sp_data$n_fp_total, na.rm = TRUE)
    max_fp <- max(sp_data$n_fp_total, na.rm = TRUE)
    min_fn <- min(sp_data$n_fn, na.rm = TRUE)
    max_fn <- max(sp_data$n_fn, na.rm = TRUE)

    message(sprintf("  %s: %d-%d genes, FP: %.0f-%.0f, FN: %.0f-%.0f",
                    species_labels[sp], min_genes, max_genes, min_fp, max_fp, min_fn, max_fn))
  }
}

message("\n=== COMPLETE ===")
