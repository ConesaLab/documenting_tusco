#!/usr/bin/env Rscript
#####################################
#####       TUSCO Report       ######
#####################################
### Author: Tianyuan Liu
### Last Modified: 11/07/2024
###
### Mouse version using ComplexUpset for the FN UpSet plot
### (Adapted from human_FN_validation.R)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do not change working directory; restrict outputs to local ./plot and ./tsv
# Mitigate OpenMP SHM issues in restricted environments
Sys.setenv(OMP_NUM_THREADS = "1", OMP_PROC_BIND = "FALSE", OMP_WAIT_POLICY = "PASSIVE", KMP_INIT_AT_FORK = "0")

# Helper to resolve preferred paths: prefer new data layout, fallback to legacy
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && base::file.exists(p)) return(p)
    if (is_dir && base::dir.exists(p)) return(p)
  }
  return(candidates[[1]])
}

suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  # ComplexUpset is optional; we fall back if unavailable
  # library(ComplexUpset)
  library(ggplot2)
})
has_ComplexUpset <- requireNamespace("ComplexUpset", quietly = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Define file paths/folders
#    (EDIT paths to match your local environment)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Example argument usage (uncomment if you want to take from commandArgs)
# args <- commandArgs(trailingOnly = TRUE)
# class_file <- args[1]
# tusco_file <- args[2] # Changed from bugsi_file
# transcript_gtf_file <- args[3] # This was bugsi_file in mouse, now tusco_file is TSV
# utilities_path <- args[4]

# Data roots (new structure under data/, fallback to legacy figs/data)
parent_dir <- resolve_path(c(
  base::file.path("..", "..", "..", "data", "raw", "lrgasp", "mouse"),
  base::file.path("..", "data", "lrgasp", "mouse")
), is_dir = TRUE)

classification_folders <- c(
  "ES_cdna_ont",
  "ES_cdna_ont_ls",
  "ES_cdna_pacbio",
  "ES_cdna_pacbio_ls",
  "ES_drna_ont",
  "ES_drna_ont_ls"
)

classification_filename <- "_classification.txt"

# TUSCO TSV file for mouse (processed output)
tusco_file <- resolve_path(c(
  base::file.path("..", "..", "..", "data", "processed", "tusco", "mmu", "tusco_mouse.tsv"),
  base::file.path("..", "data", "tusco", "tusco_mouse.tsv")
))

# Short-read quant
short_read_quant <- resolve_path(c(
  base::file.path("..", "..", "..", "data", "raw", "lrgasp", "short_read_quant", "mouse_quant.genes.sf"),
  base::file.path("..", "data", "lrgasp", "short_read_quant", "mouse_quant.genes.sf")
))

# Output folders (project-relative under figs/figure-03)
plot_dir <- base::file.path("..", "plots")
tsv_dir  <- base::file.path("..", "tables")
if (!base::dir.exists(plot_dir)) base::dir.create(plot_dir, recursive = TRUE)
if (!base::dir.exists(tsv_dir))  base::dir.create(tsv_dir,  recursive = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Helper: read TUSCO TSV, get TUSCO gene set (no version)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_tusco_genes <- function(tsv_path) {
  if (!base::file.exists(tsv_path)) {
    base::stop("TUSCO TSV file not found at: ", tsv_path)
  }
  base::message("Reading TUSCO TSV from: ", tsv_path)

  # Robust read, similar to human script
  tusco_df <- readr::read_delim(
    tsv_path,
    delim = "\\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"), # Assuming same cols as human
    col_types = readr::cols(.default = "c"),
    trim_ws = TRUE,
    comment = "#" # Allow comments in TUSCO TSV if any
  )
  if (base::ncol(tusco_df) == 1 && base::nrow(tusco_df) > 0) { # Fallback if not tab-delimited and not empty
    base::message("Falling back to space-delimited reading for TUSCO TSV.")
    # Re-read header to get actual column names if first attempt didn't get them
    # This part might need adjustment based on actual mouse tusco file format if it has headers and is space delimited.
    # For now, assuming the col_names are standard or it's truly headerless.
    tusco_df <- readr::read_table2(
      tsv_path,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = readr::cols(.default = "c"),
      comment = "#"
    )
  }

  if (!"ensembl" %in% base::colnames(tusco_df)) {
    base::stop("TUSCO TSV must contain 'ensembl' column. Current columns: ", paste(colnames(tusco_df), collapse=", "))
  }
  
  tusco_genes_raw <- tusco_df %>%
    dplyr::mutate(ensembl_id = stringr::str_remove(ensembl, "\\..+$")) %>% # remove version from ensembl ID
    dplyr::distinct(ensembl_id) %>%
    dplyr::pull(ensembl_id)
  
  # Filter out invalid gene IDs like "#", empty strings, or NA
  tusco_genes_filtered <- tusco_genes_raw[!tusco_genes_raw %in% c("#", "") & !is.na(tusco_genes_raw) & tusco_genes_raw != "ensembl"] # Added "ensembl" in case header is read as data
  
  if (base::length(tusco_genes_filtered) == 0) {
    base::warning("No TUSCO genes found after filtering. Check TSV content or parsing logic.")
  }
  
  return(tusco_genes_filtered)
}

tusco_genes <- read_tusco_genes(tusco_file) # Changed from bugsi_genes
base::message("Total TUSCO genes: ", base::length(tusco_genes))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Define short-read FNs
#    A TUSCO gene is found if TPM > 0 in short_read_quant
#    => FN = TUSCO set - found set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!base::file.exists(short_read_quant)) {
  base::stop("Short-read quant file not found at: ", short_read_quant)
}
base::message("Reading short-read quant from: ", short_read_quant)

short_df <- readr::read_tsv(short_read_quant, comment = "#") %>%
  dplyr::mutate(gene_id = stringr::str_remove(Name, "\\..+$")) # remove version

if (!base::all(c("Name", "TPM") %in% base::colnames(short_df))) {
  base::stop("short_df is missing expected columns (Name or TPM). Check the file format.")
}

short_found <- short_df %>%
  dplyr::filter(TPM > 0) %>%
  dplyr::distinct(gene_id) %>%
  dplyr::pull(gene_id)

short_tusco_found <- base::intersect(tusco_genes, short_found) # Changed from short_bugsi_found
short_tusco_fn <- base::setdiff(tusco_genes, short_tusco_found) # Changed from short_bugsi_fn
base::message("Short-read FN (TUSCO): ", base::length(short_tusco_fn))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Define a function to read a classification file,
#    then identify which TUSCO genes are "found" by FSM/ISM.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_tusco_fsm_ism_found <- function(class_file, tusco_set) { # Changed from get_bugsi_fsm_ism_found
  if (!base::file.exists(class_file)) {
    base::stop("Classification file not found: ", class_file)
  }
  
  base::message("Reading classification file: ", class_file)
  df <- readr::read_tsv(class_file, col_types = readr::cols(.default = "c"))
  
  required_cols <- c("isoform", "chrom", "strand", "length", "exons", 
                     "structural_category", "associated_gene")
  
  if (base::ncol(df) < 7) {
    base::stop("Expected at least 7 columns in classification file, got: ", base::ncol(df))
  }
  
  base::colnames(df)[1:7] <- required_cols
  
  if (!base::all(c("structural_category", "associated_gene") %in% base::colnames(df))) {
    base::stop("Column naming error: structural_category or associated_gene not found.")
  }
  
  df_clean <- df %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(
      df %>%
        dplyr::filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_gene, sep = "_")
    )
  
  df_clean <- df_clean %>%
    dplyr::mutate(associated_gene = stringr::str_remove(associated_gene, "\\.\\d+$"))
  
  found_genes <- df_clean %>%
    dplyr::filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    dplyr::distinct(associated_gene) %>%
    dplyr::pull(associated_gene)
  
  tusco_found <- base::intersect(tusco_set, found_genes) # Changed from bugsi_found
  
  base::message("  -> Found in classification file: ", base::length(tusco_found),
                " (FSM/ISM) out of ", base::length(tusco_set), " total TUSCO genes.") # Changed from BUGSI
  
  return(tusco_found)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) For each classification folder, define FN = TUSCO - found
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fn_list <- base::list()

for (folder in classification_folders) {
  class_file <- base::file.path(parent_dir, folder, base::paste0(folder, classification_filename))
  found_genes_this <- get_tusco_fsm_ism_found(class_file, tusco_genes) # Changed from bugsi_genes
  fn_this <- base::setdiff(tusco_genes, found_genes_this) # Changed from bugsi_genes
  fn_list[[folder]] <- fn_this
  
  base::message(folder, " => Found (TUSCO): ", base::length(found_genes_this), # Changed from BUGSI
                ", FN (TUSCO): ", base::length(fn_this)) # Changed from BUGSI
}

fn_list[["short_read"]] <- short_tusco_fn # Changed from short_bugsi_fn

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7) Improve labels for UpSet
#    Remove ES_ prefix, underscores, fix ONT->ONT, pacbio->PacBio, etc.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pretty_label <- function(x) {
  out <- stringr::str_remove(x, "^ES_") # Remove ES_ prefix
  out <- stringr::str_replace_all(out, "_", " ")
  out <- stringr::str_replace_all(out, "(?i)ont", "ONT")
  out <- stringr::str_replace_all(out, "(?i)pacbio", "PacBio")
  out <- stringr::str_replace(out, "short read", "Short Read") # Capitalize Short Read
  return(out)
}

original_names <- base::names(fn_list)
nice_names <- base::sapply(original_names, pretty_label, USE.NAMES = FALSE)

base::names(fn_list) <- nice_names

if (base::length(fn_list) == 0) {
  base::warning("fn_list is empty. Possibly no classification folders or short_read. Check your inputs.")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8) Create an UpSet plot using ComplexUpset & improve the look
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert list to binary incidence data.frame (no UpSetR dependency)
list_to_incidence_df <- function(x) {
  all_items <- unique(unlist(x, use.names = FALSE))
  if (length(all_items) == 0) return(data.frame())
  df <- data.frame(item = all_items, check.names = FALSE)
  for (nm in names(x)) {
    df[[nm]] <- df$item %in% x[[nm]]
  }
  df
}
upset_data <- list_to_incidence_df(fn_list)

sets_order <- base::rev(base::names(fn_list))

if (has_ComplexUpset && nrow(upset_data) > 0) {
  p <- ComplexUpset::upset(
    upset_data,
    width_ratio = 0.2,
    height_ratio = 0.7, # Matched human script for consistency
    intersect = sets_order,
    keep_empty_groups = TRUE,
    wrap = TRUE,
    set_sizes = FALSE,
    base_annotations = list(
      'Intersect\nFN TUSCO\nGenes' = ComplexUpset::intersection_size(counts = TRUE, text_size = 7, point_size = 2)
    )
  ) +
    ggplot2::labs(
      title = "False Negative TUSCO Genes UpSet Plot",
      subtitle = NULL,
      caption = NULL,
      y = "TUSCO\n(Mouse)"
    ) +
    ggplot2::theme_classic(base_family = "Helvetica", base_size = 7) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 7),
      axis.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(size = 7, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 7),
      plot.caption = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_text(size = 7, angle = 90, vjust = 0.5),
      strip.text.x = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      legend.position = "none"
    )
} else {
  p <- NULL
  base::message("ComplexUpset not available or empty data; skipping UpSet plot.")
}

## No separate output_dir; use plot_dir/tsv_dir defined above

if (FALSE) ggplot2::ggsave(
  filename = base::file.path(plot_dir, "figure3d-mouse.pdf"), # Disabled per request
  plot = p,
  width = 4.5,  # Matched human script
  height = 3.5, # Matched human script
  units = "in",
  device = grDevices::cairo_pdf, # Explicitly use grDevices::cairo_pdf
  dpi = 300
)

base::message("UpSet plot step completed (plot may be skipped).")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9) Print out the intersection of all sets (common FNs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
common_fn <- base::Reduce(base::intersect, fn_list)
base::message("Number of common FN TUSCO genes across all sets: ", base::length(common_fn)) # Changed from BUGSI

# Do not write extra files outside the consolidated TSV export

base::message("Done.")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 10) Plot distribution of FNs by number of samples shared
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine all FN genes from all samples into a single vector
all_fn_genes <- base::unlist(fn_list)

# Count how many times each gene appears (i.e., in how many samples it's an FN)
gene_counts <- base::table(all_fn_genes)

# Count how many genes appear in exactly 1 sample, 2 samples, etc.
sample_counts <- base::table(gene_counts)

# Convert to a data frame for plotting
sample_counts_df <- base::as.data.frame(sample_counts)
base::colnames(sample_counts_df) <- c("num_samples_shared", "num_genes")

# Ensure num_samples_shared is numeric (factor by default from table)
sample_counts_df$num_samples_shared <- base::as.numeric(base::as.character(sample_counts_df$num_samples_shared))

# Ensure all sample counts from 1 to N are present, even if count is 0
num_total_samples <- base::length(fn_list)
all_possible_counts <- base::data.frame(num_samples_shared = 1:num_total_samples)

plot_df <- base::merge(all_possible_counts, sample_counts_df, by = "num_samples_shared", all.x = TRUE)
plot_df$num_genes[base::is.na(plot_df$num_genes)] <- 0 # Replace NA with 0 for missing counts


# Create the bar plot
p_dist <- ggplot2::ggplot(plot_df, ggplot2::aes(x = base::factor(num_samples_shared), y = num_genes)) +
  ggplot2::geom_bar(stat = "identity", fill = "#1C9E77") + # Use specified color (same as original mouse script)
  ggplot2::labs(
    title = NULL,
    x = "#Mouse Sample", 
    y = "FN TUSCO Genes" # Changed from BUGSI
  ) +
  # Set fixed y-axis scale and breaks
  ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15)) + # Matched human script
  ggplot2::theme_minimal(base_size = 7) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 7),
    axis.text = ggplot2::element_text(size = 7),
    axis.title = ggplot2::element_text(size = 7),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    plot.margin = ggplot2::margin(5, 5, 5, 5),
    panel.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  # Add text labels on top of bars
  ggplot2::geom_text(ggplot2::aes(label = num_genes), vjust = -0.5, size = 2)

# Save the plot
dist_plot_file <- base::file.path(plot_dir, "figure3d-mouse.pdf") # Output only under ./plot
ggplot2::ggsave(
  filename = dist_plot_file,
  plot = p_dist,
  width = 2,    
  height = 1.2, 
  units = "in",
  device = grDevices::cairo_pdf, # Explicitly use grDevices::cairo_pdf
  dpi = 300
)

base::message("FN TUSCO distribution bar plot saved to: ", dist_plot_file)

# Export TSV with underlying data and minimal metadata
tsv_path <- base::file.path(tsv_dir, "figure3d-mouse.tsv")
plot_df_out <- plot_df %>% dplyr::mutate(figure_id = "fig-3", panel_id = "3d-mouse", record_type = "distribution_counts")
readr::write_tsv(plot_df_out, tsv_path)
