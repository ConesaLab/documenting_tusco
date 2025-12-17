#!/usr/bin/env Rscript
# alphagenome_expression_distribution.R
# Reviewer response figure: AlphaGenome expression distribution
# Shows TUSCO genes among top expressers in single-isoform gene pool

###############################################################
# 0. Setup and Dependencies
###############################################################

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# Get script directory and navigate to repo root
get_script_dir <- function() {
  # Try multiple methods to determine script location
  # Method 1: commandArgs
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  # Method 2: sys.frame (for source())
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile)))
    }
  }
  # Method 3: fallback to current directory
  return(getwd())
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)
setwd(repo_root)

# Source common utilities
source("scripts/figure_utils.R")

# Output directory
out_dir <- "reviewer_response/review_plots/alphagenome"
tables_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################
# 1. Define File Paths
###############################################################

# AlphaGenome JSON files
alphagenome_human_json <- "data/hsa/alphagenome_human.json"
alphagenome_mouse_json <- "data/mmu/alphagenome_mouse.json"

# TUSCO gene lists (final selected)
tusco_human_tsv <- "data/processed/tusco/tusco_human.tsv"
tusco_mouse_tsv <- "data/processed/tusco/tusco_mouse.tsv"

###############################################################
# 2. Helper Functions
###############################################################

#' Read TUSCO gene IDs from TSV (no headers, gene_id in column 1)
#' @param tsv_path Path to TSV file
#' @return Character vector of Ensembl gene IDs (version stripped)
read_gene_ids <- function(tsv_path) {
  lines <- readLines(tsv_path)
  # Remove comment lines
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) return(character(0))

  dt <- fread(text = paste(lines, collapse = "\n"),
              sep = "\t", header = FALSE)
  gene_ids <- dt[[1]]
  # Strip version suffix
  gene_ids <- gsub("\\..*$", "", gene_ids)
  unique(gene_ids)
}

#' Parse AlphaGenome JSON and compute median expression per gene
#' @param json_path Path to AlphaGenome JSON file
#' @return data.table with gene_id and median_expression columns
parse_alphagenome_expression <- function(json_path) {
  message("Parsing AlphaGenome JSON: ", json_path)

  records <- fromJSON(json_path, simplifyVector = FALSE)
  n <- length(records)

  gene_id <- character(n)
  median_expression <- numeric(n)
  n_tissues <- integer(n)

  for (i in seq_len(n)) {
    rec <- records[[i]]

    # Extract and clean gene ID
    gid <- rec$gene_id
    if (is.null(gid)) gid <- NA_character_
    gene_id[i] <- gsub("\\..*$", "", gid)

    # Extract tissue expression values
    expr_dict <- rec$tissue_expression
    if (!is.null(expr_dict) && length(expr_dict) > 0) {
      expr_vals <- as.numeric(unlist(expr_dict, use.names = FALSE))
      expr_vals <- expr_vals[is.finite(expr_vals)]

      median_expression[i] <- if (length(expr_vals) > 0) {
        median(expr_vals, na.rm = TRUE)
      } else {
        NA_real_
      }
      n_tissues[i] <- length(expr_vals)
    } else {
      median_expression[i] <- NA_real_
      n_tissues[i] <- 0L
    }
  }

  data.table(
    gene_id = gene_id,
    median_expression = median_expression,
    n_tissues = n_tissues
  )
}

###############################################################
# 3. Load and Process Data
###############################################################

message("Loading gene lists...")

# Load TUSCO gene IDs
tusco_human <- read_gene_ids(tusco_human_tsv)
tusco_mouse <- read_gene_ids(tusco_mouse_tsv)

message("  Human TUSCO genes: ", length(tusco_human))
message("  Mouse TUSCO genes: ", length(tusco_mouse))

# Parse AlphaGenome expression data
human_expr <- parse_alphagenome_expression(alphagenome_human_json)
mouse_expr <- parse_alphagenome_expression(alphagenome_mouse_json)

message("  Human AlphaGenome genes: ", nrow(human_expr))
message("  Mouse AlphaGenome genes: ", nrow(mouse_expr))

###############################################################
# 4. Annotate and Combine Data
###############################################################

# Annotate human data
human_expr[, `:=`(
  species = "Human",
  is_tusco = gene_id %in% tusco_human
)]

# Annotate mouse data
mouse_expr[, `:=`(
  species = "Mouse",
  is_tusco = gene_id %in% tusco_mouse
)]

# Combine species
combined_dt <- rbindlist(list(human_expr, mouse_expr), use.names = TRUE)

# Filter to valid expression values
plot_dt <- combined_dt[
  !is.na(median_expression) &
  median_expression > 0
]

# Create group labels
plot_dt[, group := ifelse(is_tusco, "TUSCO", "Other single-isoform")]
plot_dt[, group := factor(group, levels = c("Other single-isoform", "TUSCO"))]
plot_dt[, species := factor(species, levels = c("Human", "Mouse"))]

# Log10 transform expression
LOG_OFFSET <- 1e-6
plot_dt[, log10_expression := log10(median_expression + LOG_OFFSET)]

message("\nData summary:")
message("  Total genes for plotting: ", nrow(plot_dt))
print(plot_dt[, .N, by = .(species, group)])

###############################################################
# 5. Calculate Summary Statistics
###############################################################

# Percentile ranks of TUSCO genes within single-isoform pool
calc_percentile_rank <- function(dt, sp) {
  sp_dt <- dt[species == sp]
  all_expr <- sp_dt$median_expression
  tusco_expr <- sp_dt[is_tusco == TRUE, median_expression]

  # Calculate percentile for each TUSCO gene
  percentiles <- sapply(tusco_expr, function(x) {
    ecdf(all_expr)(x) * 100
  })

  data.table(
    species = sp,
    n_tusco_genes = length(tusco_expr),
    n_total_genes = length(all_expr),
    tusco_median_percentile = median(percentiles, na.rm = TRUE),
    tusco_min_percentile = min(percentiles, na.rm = TRUE),
    tusco_max_percentile = max(percentiles, na.rm = TRUE),
    tusco_q25_percentile = quantile(percentiles, 0.25, na.rm = TRUE),
    tusco_q75_percentile = quantile(percentiles, 0.75, na.rm = TRUE)
  )
}

percentile_stats <- rbindlist(list(
  calc_percentile_rank(plot_dt, "Human"),
  calc_percentile_rank(plot_dt, "Mouse")
))

# Expression summary by group
expression_summary <- plot_dt[, .(
  n_genes = .N,
  median_rpkm = median(median_expression, na.rm = TRUE),
  mean_rpkm = mean(median_expression, na.rm = TRUE),
  sd_rpkm = sd(median_expression, na.rm = TRUE),
  q25_rpkm = quantile(median_expression, 0.25, na.rm = TRUE),
  q75_rpkm = quantile(median_expression, 0.75, na.rm = TRUE),
  min_rpkm = min(median_expression, na.rm = TRUE),
  max_rpkm = max(median_expression, na.rm = TRUE),
  median_log10 = median(log10_expression, na.rm = TRUE)
), by = .(species, group)]

message("\nExpression summary:")
print(expression_summary)

message("\nTUSCO percentile ranks:")
print(percentile_stats)

###############################################################
# 6. Create Visualization
###############################################################

# Define colors
plot_colors <- c(
  "Other single-isoform" = "#D1D3D4",  # Light grey
  "TUSCO" = "#1b9e77"                   # TUSCO green
)

# Get TUSCO median for vertical lines
tusco_medians <- expression_summary[group == "TUSCO", .(species, median_log10)]

# Main density plot with TUSCO overlay
p <- ggplot() +
  # Background density for all genes (grey)
  geom_density(
    data = plot_dt[is_tusco == FALSE],
    aes(x = log10_expression),
    fill = "#D1D3D4",
    color = "#999999",
    alpha = 0.7,
    linewidth = 0.3
  ) +
  # Overlaid density for TUSCO genes (green)
  geom_density(
    data = plot_dt[is_tusco == TRUE],
    aes(x = log10_expression),
    fill = TUSCO_COLORS$human_tusco,
    color = "#1b9e77",
    alpha = 0.7,
    linewidth = 0.5
  ) +
  # Rug marks for TUSCO genes
  geom_rug(
    data = plot_dt[is_tusco == TRUE],
    aes(x = log10_expression),
    color = "#1b9e77",
    alpha = 0.9,
    length = unit(0.04, "npc"),
    sides = "b",
    linewidth = 0.5
  ) +
  # Vertical lines at TUSCO median expression
  geom_vline(
    data = tusco_medians,
    aes(xintercept = median_log10),
    color = "#1b9e77",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  facet_wrap(~species, nrow = 1, scales = "fixed") +
  scale_x_continuous(
    breaks = seq(-2, 2, by = 1),
    labels = math_format(10^.x)
  ) +
  labs(
    x = "Median RPKM across tissues",
    y = "Density"
  ) +
  theme_tusco(base_size = 7) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(1, "lines"),
    axis.title.x = element_text(margin = margin(t = 5))
  )

# Add percentile annotation
annotation_df <- percentile_stats[, .(
  species = factor(species, levels = c("Human", "Mouse")),
  label = sprintf("TUSCO: %.0f%% percentile\n(n=%d)",
                  tusco_median_percentile, n_tusco_genes)
)]

p <- p +
  geom_text(
    data = annotation_df,
    aes(x = -1.5, y = Inf, label = label),
    hjust = 0, vjust = 1.2,
    size = 2.2,
    color = "#1b9e77",
    fontface = "bold"
  )

###############################################################
# 7. Alternative: Histogram with density overlay
###############################################################

p_hist <- ggplot() +
  # Histogram for non-TUSCO genes
  geom_histogram(
    data = plot_dt[is_tusco == FALSE],
    aes(x = log10_expression, y = after_stat(density)),
    fill = "#D1D3D4",
    color = "white",
    bins = 40,
    alpha = 0.8
  ) +
  # Density curve for TUSCO genes
  geom_density(
    data = plot_dt[is_tusco == TRUE],
    aes(x = log10_expression),
    fill = TUSCO_COLORS$human_tusco,
    color = "#1b9e77",
    alpha = 0.6,
    linewidth = 0.7
  ) +
  # Rug marks for TUSCO genes
  geom_rug(
    data = plot_dt[is_tusco == TRUE],
    aes(x = log10_expression),
    color = "#1b9e77",
    alpha = 0.9,
    length = unit(0.04, "npc"),
    sides = "b",
    linewidth = 0.5
  ) +
  # Vertical lines at TUSCO median
  geom_vline(
    data = tusco_medians,
    aes(xintercept = median_log10),
    color = "#1b9e77",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  facet_wrap(~species, nrow = 1, scales = "fixed") +
  scale_x_continuous(
    breaks = seq(-2, 2, by = 1),
    labels = math_format(10^.x)
  ) +
  labs(
    x = "Median RPKM across tissues",
    y = "Density"
  ) +
  theme_tusco(base_size = 7) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(1, "lines"),
    axis.title.x = element_text(margin = margin(t = 5))
  ) +
  geom_text(
    data = annotation_df,
    aes(x = -1.5, y = Inf, label = label),
    hjust = 0, vjust = 1.2,
    size = 2.2,
    color = "#1b9e77",
    fontface = "bold"
  )

###############################################################
# 8. Save Outputs
###############################################################

# Save density plot
ggsave(
  file.path(out_dir, "alphagenome_expression_density.pdf"),
  p,
  width = 6,
  height = 2.5,
  units = "in",
  device = cairo_pdf
)
message("\nSaved: ", file.path(out_dir, "alphagenome_expression_density.pdf"))

# Save histogram plot
ggsave(
  file.path(out_dir, "alphagenome_expression_histogram.pdf"),
  p_hist,
  width = 6,
  height = 2.5,
  units = "in",
  device = cairo_pdf
)
message("Saved: ", file.path(out_dir, "alphagenome_expression_histogram.pdf"))

# Save summary statistics
fwrite(
  expression_summary,
  file.path(tables_dir, "alphagenome_expression_summary.tsv"),
  sep = "\t"
)
message("Saved: ", file.path(tables_dir, "alphagenome_expression_summary.tsv"))

# Save percentile statistics
fwrite(
  percentile_stats,
  file.path(tables_dir, "alphagenome_tusco_percentile_ranks.tsv"),
  sep = "\t"
)
message("Saved: ", file.path(tables_dir, "alphagenome_tusco_percentile_ranks.tsv"))

# Save full plot data
fwrite(
  plot_dt[, .(species, gene_id, group, is_tusco, median_expression,
              log10_expression, n_tissues)],
  file.path(tables_dir, "alphagenome_expression_data.tsv"),
  sep = "\t"
)
message("Saved: ", file.path(tables_dir, "alphagenome_expression_data.tsv"))

message("\nDone!")
