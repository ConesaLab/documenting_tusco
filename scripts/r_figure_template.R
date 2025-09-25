#!/usr/bin/env Rscript

# TUSCO Figure Script Template
# This template provides a standardized structure for all figure generation scripts
# Usage: Rscript figure_script.R [output_dir] [width] [height]

# =============================================================================
# 1. COMMAND LINE ARGUMENT PARSING
# =============================================================================

# Source utilities library
if (file.exists("../../scripts/figure_utils.R")) {
  source("../../scripts/figure_utils.R")
} else if (file.exists("../../../scripts/figure_utils.R")) {
  source("../../../scripts/figure_utils.R")
} else {
  # Fallback - define minimal functions inline
  parse_figure_args <- function(defaults = list(out_dir = "..", width = 3.35, height = 4.0)) {
    args <- commandArgs(trailingOnly = TRUE)
    result <- defaults
    if (length(args) > 0) result$out_dir <- args[1]
    if (length(args) > 1) result$width <- as.numeric(args[2])
    if (length(args) > 2) result$height <- as.numeric(args[3])
    return(result)
  }
  
  create_output_dirs <- function(base_dir = ".") {
    plot_dir <- file.path(base_dir, "plots")
    table_dir <- file.path(base_dir, "tables")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    list(plot_dir = plot_dir, table_dir = table_dir)
  }
  
  save_figure <- function(plot_obj, fig_id, plot_dir = "plots", table_dir = "tables", 
                         width = 3.35, height = 4.0, data = NULL, data_suffix = NULL) {
    plot_file <- file.path(plot_dir, paste0(fig_id, ".pdf"))
    ggsave(plot_file, plot_obj, width = width, height = height, 
           units = "in", device = "pdf", dpi = 300)
    message("Saved plot: ", plot_file)
    if (!is.null(data)) {
      data_name <- if (is.null(data_suffix)) fig_id else paste0(fig_id, "-", data_suffix)
      data_file <- file.path(table_dir, paste0(data_name, ".tsv"))
      readr::write_tsv(data, data_file)
      message("Saved data: ", data_file)
    }
    invisible(plot_file)
  }
}

# Parse arguments
params <- parse_figure_args(defaults = list(
  out_dir = "..",           # Output to parent directory (figure-XX/)
  width = 3.35,            # Nature single column width
  height = 4.0             # Default height
))

# =============================================================================
# 2. PACKAGE LOADING
# =============================================================================

message("Loading required packages...")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  # Add other required packages here
})

# =============================================================================
# 3. OUTPUT DIRECTORY SETUP
# =============================================================================

# Create output directories
dirs <- create_output_dirs(params$out_dir)
plot_dir <- dirs$plot_dir
table_dir <- dirs$table_dir

message("Output directories:")
message("  Plots: ", plot_dir)
message("  Tables: ", table_dir)

# =============================================================================
# 4. DATA INPUT AND PATH RESOLUTION
# =============================================================================

# Define data file paths using resolve_data_path if available, or manual resolution
if (exists("resolve_data_path")) {
  # Use utility function
  example_data_file <- resolve_data_path("example", "data.tsv")
} else {
  # Manual resolution fallback
  data_candidates <- c(
    "../../data/raw/example/data.tsv",
    "../../../data/raw/example/data.tsv", 
    "../../data/processed/example/data.tsv"
  )
  example_data_file <- data_candidates[file.exists(data_candidates)][1]
  if (is.na(example_data_file)) stop("Required data file not found")
}

# Read data files
# example_data <- read_tsv(example_data_file, show_col_types = FALSE)

# =============================================================================
# 5. DATA PROCESSING AND ANALYSIS
# =============================================================================

# Add your data processing code here
message("Processing data...")

# Example: Create sample plot data
sample_data <- data.frame(
  x = 1:10,
  y = rnorm(10),
  group = rep(c("A", "B"), each = 5)
)

# =============================================================================
# 6. PLOT GENERATION
# =============================================================================

# Define consistent theme
theme_figure <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25)
  )

# Create main plot
main_plot <- ggplot(sample_data, aes(x = x, y = y, color = group)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Example Figure",
    x = "X Value", 
    y = "Y Value",
    color = "Group"
  ) +
  theme_figure

# =============================================================================
# 7. OUTPUT GENERATION
# =============================================================================

# Save figure and data
save_figure(
  plot_obj = main_plot,
  fig_id = "fig-example",  # Change this to match your figure
  plot_dir = plot_dir,
  table_dir = table_dir,
  width = params$width,
  height = params$height,
  data = sample_data,
  data_suffix = "raw-data"
)

# Optional: Save additional outputs
summary_data <- sample_data %>%
  group_by(group) %>%
  summarize(
    mean_y = mean(y),
    sd_y = sd(y),
    .groups = "drop"
  )

save_figure(
  plot_obj = NULL,  # No plot, just data
  fig_id = "fig-example", 
  plot_dir = plot_dir,
  table_dir = table_dir,
  data = summary_data,
  data_suffix = "summary"
)

message("Figure generation completed successfully!")

# =============================================================================
# 8. SESSION INFO (OPTIONAL)
# =============================================================================

if (interactive()) {
  message("Session info:")
  sessionInfo()
}