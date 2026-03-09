#!/usr/bin/env Rscript
# ==============================================================================
# TUSCO vs MAJIQ-L Benchmarking Comparison Analysis
# Supplementary Figure S4
# ==============================================================================

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
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()
output_base <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
plot_dir <- file.path(output_base, "plots")
table_dir <- file.path(output_base, "tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("Script directory:", script_dir, "\n")
cat("Output directory:", plot_dir, "\n")

# ------------------------------------------------------------------------------
# Load pre-computed comparison data
# ------------------------------------------------------------------------------
cat("\n=== Loading comparison data ===\n")

comparison_file <- file.path(table_dir, "comparison_metrics.tsv")
if (!file.exists(comparison_file)) {
  stop("Comparison data file not found: ", comparison_file)
}

df <- read_tsv(comparison_file, show_col_types = FALSE)
cat("Loaded", nrow(df), "pipelines\n")
print(df)

# Create long format data for plotting
tusco_data <- df %>%
  select(pipeline_label, 
         sensitivity = sensitivity_TUSCO, 
         fn_rate = fn_rate_TUSCO, 
         fdr = fdr_TUSCO) %>%
  mutate(method = "TUSCO")

majiq_data <- df %>%
  select(pipeline_label, 
         sensitivity = `sensitivity_MAJIQ-L`, 
         fn_rate = `fn_rate_MAJIQ-L`, 
         fdr = `fdr_MAJIQ-L`) %>%
  mutate(method = "MAJIQ-L")

comparison_data <- bind_rows(tusco_data, majiq_data)

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

METHOD_COLORS <- c("TUSCO" = "#a8d5a0", "MAJIQ-L" = "#B2182B")

# Order pipelines by TUSCO sensitivity (best to worst)
pipeline_order <- df %>%
  arrange(desc(sensitivity_TUSCO)) %>%
  pull(pipeline_label)

comparison_data <- comparison_data %>%
  mutate(
    pipeline_label = factor(pipeline_label, levels = pipeline_order),
    method = factor(method, levels = c("TUSCO", "MAJIQ-L"))
  )

df <- df %>%
  mutate(pipeline_label = factor(pipeline_label, levels = pipeline_order))

# ------------------------------------------------------------------------------
# Create figure: Three-panel dot plot with connecting lines
# ------------------------------------------------------------------------------
cat("\n=== Creating comparison figure ===\n")

# Panel A: Sensitivity
p_sensitivity <- ggplot() +
  geom_segment(
    data = df,
    aes(x = pipeline_label, xend = pipeline_label,
        y = sensitivity_TUSCO, yend = `sensitivity_MAJIQ-L`),
    color = "gray60", linewidth = 0.8, alpha = 0.6
  ) +
  geom_point(
    data = comparison_data,
    aes(x = pipeline_label, y = sensitivity, color = method, shape = method),
    size = 4, alpha = 0.9
  ) +
  scale_color_manual(values = METHOD_COLORS, name = "Method",
                     labels = c("TUSCO" = "TUSCO evaluation", "MAJIQ-L" = "MAJIQ-L")) +
  scale_shape_manual(values = c("TUSCO" = 16, "MAJIQ-L" = 17), name = "Method",
                     labels = c("TUSCO" = "TUSCO evaluation", "MAJIQ-L" = "MAJIQ-L")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(title = "Sensitivity", x = NULL, y = "Sensitivity (%)") +
  theme_comparison() +
  theme(legend.position = "none")

# Panel B: FN Rate
p_fn_rate <- ggplot() +
  geom_segment(
    data = df,
    aes(x = pipeline_label, xend = pipeline_label,
        y = fn_rate_TUSCO, yend = `fn_rate_MAJIQ-L`),
    color = "gray60", linewidth = 0.8, alpha = 0.6
  ) +
  geom_point(
    data = comparison_data,
    aes(x = pipeline_label, y = fn_rate, color = method, shape = method),
    size = 4, alpha = 0.9
  ) +
  scale_color_manual(values = METHOD_COLORS, name = "Method",
                     labels = c("TUSCO" = "TUSCO evaluation", "MAJIQ-L" = "MAJIQ-L")) +
  scale_shape_manual(values = c("TUSCO" = 16, "MAJIQ-L" = 17), name = "Method",
                     labels = c("TUSCO" = "TUSCO evaluation", "MAJIQ-L" = "MAJIQ-L")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(title = "False Negative Rate", x = NULL, y = "FN Rate (%)") +
  theme_comparison() +
  theme(legend.position = "none")

# Panel C: FDR
p_fdr <- ggplot() +
  geom_segment(
    data = df,
    aes(x = pipeline_label, xend = pipeline_label,
        y = fdr_TUSCO, yend = `fdr_MAJIQ-L`),
    color = "gray60", linewidth = 0.8, alpha = 0.6
  ) +
  geom_point(
    data = comparison_data,
    aes(x = pipeline_label, y = fdr, color = method, shape = method),
    size = 4, alpha = 0.9
  ) +
  scale_color_manual(values = METHOD_COLORS, name = "Method",
                     labels = c("TUSCO" = "TUSCO evaluation", "MAJIQ-L" = "MAJIQ-L")) +
  scale_shape_manual(values = c("TUSCO" = 16, "MAJIQ-L" = 17), name = "Method",
                     labels = c("TUSCO" = "TUSCO evaluation", "MAJIQ-L" = "MAJIQ-L")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(title = "False Discovery Rate", x = NULL, y = "FDR (%)") +
  theme_comparison() +
  theme(legend.position = "none")

# Combine panels with shared legend
p_combined <- (p_sensitivity | p_fn_rate | p_fdr) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5))) &
  theme(legend.position = "bottom")

# Save figure
ggsave(file.path(plot_dir, "fig-s4.pdf"), p_combined, width = 12, height = 5, device = "pdf")
cat("Saved figure:", file.path(plot_dir, "fig-s4.pdf"), "\n")

ggsave(file.path(plot_dir, "fig-s4.png"), p_combined, width = 12, height = 5, dpi = 300)
cat("Saved PNG:", file.path(plot_dir, "fig-s4.png"), "\n")

cat("\n=== Analysis complete ===\n")
