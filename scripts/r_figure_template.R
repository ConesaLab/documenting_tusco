#!/usr/bin/env Rscript

# TUSCO Figure Script Template
# This template provides a standardized structure for all figure generation scripts
# Usage: Rscript figure_script.R [output_dir] [width] [height]

# =============================================================================
# 1. COMMAND LINE ARGUMENT PARSING
# =============================================================================

# Robust repo root discovery - works from any working directory
find_repo_root_for_utils <- function(start = getwd(), limit = 10) {
  cur <- tryCatch(normalizePath(start, winslash = "/", mustWork = FALSE), error = function(e) start)
  for (i in seq_len(limit)) {
    if (file.exists(file.path(cur, "scripts", "figure_utils.R"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not find repo root containing scripts/figure_utils.R starting from: ", start)
}

# Detect script directory (handles Rscript and source())
script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- if (length(script_path) == 1 && nzchar(script_path)) {
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
repo_root <- find_repo_root_for_utils(script_dir)
source(file.path(repo_root, "scripts", "figure_utils.R"))

# Parse arguments
params <- parse_figure_args(defaults = list(
  out_dir = "..", # Output to parent directory (figure-XX/)
  width = 3.35, # Nature single column width
  height = 4.0 # Default height
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

# Optional: Download data from S3 if missing
if (!file.exists(example_data_file) && exists("download_s3_parallel")) {
  message("Data file not found: ", example_data_file)
  message("Attempting download from S3...")

  # Configure your S3 credentials and path
  # download_s3_parallel(
  #   bucket = "your-bucket-name",
  #   key = "path/to/data.tsv",
  #   outfile = example_data_file,
  #   ak = Sys.getenv("AWS_ACCESS_KEY_ID"),
  #   sk = Sys.getenv("AWS_SECRET_ACCESS_KEY")
  # )
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
  fig_id = "fig-example", # Change this to match your figure
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
  plot_obj = NULL, # No plot, just data
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
