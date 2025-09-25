#!/usr/bin/env Rscript
#####################################
#####       TUSCO Report       ######
#####################################
### Author: Tianyuan Liu
### Last Modified: 11/07/2024
###

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do not change working directory; restrict outputs to local ./plot and ./tsv
# Mitigate OpenMP SHM issues in restricted environments
Sys.setenv(OMP_NUM_THREADS = "1", OMP_PROC_BIND = "FALSE", OMP_WAIT_POLICY = "PASSIVE", KMP_INIT_AT_FORK = "0")

# Helper to resolve preferred paths: prefer new data layout, fallback to legacy
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
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
# transcript_gtf_file <- args[3]
# utilities_path <- args[4]

# Data roots (new structure under data/, fallback to legacy figs/data)
parent_dir <- resolve_path(c(
  file.path("..", "..", "..", "data", "raw", "lrgasp", "human"),
  file.path("..", "data", "lrgasp", "human")
), is_dir = TRUE)

classification_folders <- c(
  "WTC11_cdna_ont",
  "WTC11_cdna_ont_ls",
  "WTC11_cdna_pacbio",
  "WTC11_cdna_pacbio_ls",
  "WTC11_drna_ont",
  "WTC11_drna_ont_ls"
)

classification_filename <- "_classification.txt"

# TUSCO TSV file (processed output)
tusco_file <- resolve_path(c(
  file.path("..", "..", "..", "data", "processed", "tusco", "hsa", "tusco_human.tsv"),
  file.path("..", "data", "tusco", "tusco_human.tsv")
))

# Short-read quant
short_read_quant <- resolve_path(c(
  file.path("..", "..", "..", "data", "raw", "lrgasp", "short_read_quant", "human_quant.genes.sf"),
  file.path("..", "data", "lrgasp", "short_read_quant", "human_quant.genes.sf")
))

# Output folders (project-relative under figs/figure-03)
plot_dir <- file.path("..", "plots")
tsv_dir  <- file.path("..", "tables")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(tsv_dir))  dir.create(tsv_dir,  recursive = TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Helper: read TUSCO TSV, get TUSCO gene set (no version)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
read_tusco_genes <- function(tsv_path) {
  # Robust read of header‑less Tusco annotation (falls back to whitespace if tabs not found)
  # Adapted from barplot_check_human.R
  tusco_df <- read_delim(
    tsv_path,
    delim = "\\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = cols(.default = "c"),
    trim_ws = TRUE
  )
  if (ncol(tusco_df) == 1) { # Fallback if not tab-delimited
    message("Falling back to space-delimited reading for TUSCO TSV.")
    tusco_df <- read_table2(
      tsv_path,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = cols(.default = "c")
    )
  }

  if (!"ensembl" %in% colnames(tusco_df)) {
    stop("Tusco TSV must contain 'ensembl' column.")
  }
  
  tusco_genes_raw <- tusco_df %>%
    mutate(ensembl_id = str_remove(ensembl, "\\..+$")) %>% # remove version from ensembl ID
    distinct(ensembl_id) %>%
    pull(ensembl_id)
  
  # Filter out invalid gene IDs like "#", empty strings, or NA
  tusco_genes_filtered <- tusco_genes_raw[!tusco_genes_raw %in% c("#", "") & !is.na(tusco_genes_raw)]
  
  message(paste("DEBUG: Raw TUSCO gene IDs count (pre-filtering):", length(tusco_genes_raw)))
  message(paste("DEBUG: Filtered TUSCO gene IDs count (post-filtering '#', empty, NA):", length(tusco_genes_filtered)))
  
  return(tusco_genes_filtered)
}

tusco_genes <- read_tusco_genes(tusco_file) # Changed from bugsi_genes
message("Total TUSCO genes: ", length(tusco_genes))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Define short-read FNs
#    A TUSCO gene is found if TPM > 0 in short_read_quant
#    => FN = TUSCO set - found set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
short_df <- read_tsv(short_read_quant, comment="#") %>%
  mutate(gene_id = str_remove(Name, "\\..+$"))  # remove version

short_found <- short_df %>%
  filter(TPM > 0) %>%
  distinct(gene_id) %>%
  pull(gene_id)

short_tusco_found <- dplyr::intersect(tusco_genes, short_found) # Changed from short_bugsi_found
short_tusco_fn <- dplyr::setdiff(tusco_genes, short_tusco_found) # Changed from short_bugsi_fn
message("Short-read FN (TUSCO): ", length(short_tusco_fn))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Define a function to read a classification file,
#    then identify which TUSCO genes are "found" by FSM/ISM.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_tusco_fsm_ism_found <- function(class_file, tusco_set) { # Changed from get_bugsi_fsm_ism_found
  # Read classification file
  if (!file.exists(class_file)) {
    stop("Classification file not found: ", class_file)
  }
  
  # 1) Read file
  df <- read_tsv(class_file, col_types = cols(.default = "c"))
  # Some versions of SQANTI classification have these columns in a certain order:
  # isoform, chrom, strand, length, exons, structural_category, associated_gene, ...
  # We'll rename at least up to "associated_gene":
  colnames(df)[1:7] <- c("isoform","chrom","strand","length","exons",
                         "structural_category","associated_gene")
  
  # 2) Split out fusions if structural_category == "fusion"
  #    Actually, in many pipelines, the "associated_gene" might have underscores
  #    if it's a fusion. You can separate them, but let's do similarly as you do:
  df_clean <- df %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      df %>%
        filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
    )
  
  # 3) Remove version from associated_gene
  df_clean <- df_clean %>%
    mutate(associated_gene = str_remove(associated_gene, "\\.\\d+$"))
  
  # 4) "Found" = gene that appears with structural_category in {FSM, ISM}
  #    full-splice_match, incomplete-splice_match
  found_genes <- df_clean %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    distinct(associated_gene) %>%
    pull(associated_gene)
  
  # 5) Intersect with TUSCO set => these are the TUSCO genes that are "found"
  tusco_found <- dplyr::intersect(tusco_set, found_genes) # Changed from bugsi_found
  return(tusco_found)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6) For each classification folder, define FN = TUSCO - found
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fn_list <- list()

for (folder in classification_folders) {
  class_file <- file.path(parent_dir, folder, paste0(folder, classification_filename))
  found_genes_this <- get_tusco_fsm_ism_found(class_file, tusco_genes) # Changed from get_bugsi_fsm_ism_found and bugsi_genes
  fn_this <- dplyr::setdiff(tusco_genes, found_genes_this) # Changed from bugsi_genes
  fn_list[[folder]] <- fn_this
  
  message(folder, " => Found (TUSCO): ", length(found_genes_this),
          ", FN (TUSCO): ", length(fn_this))
}

# Finally add short_read to the list
fn_list[["short_read"]] <- short_tusco_fn # Changed from short_bugsi_fn


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7) Improve labels for UpSet
#    Remove underscores, fix ONT->ONT, pacbio->PacBio, etc.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pretty_label <- function(x) {
  # 0) remove WTC11 prefix if present
  out <- str_remove(x, "^WTC11_")
  # 1) remove underscores
  out <- str_replace_all(out, "_", " ")
  # 2) "ont" -> "ONT", "pacbio" -> "PacBio"
  out <- str_replace_all(out, "(?i)ont", "ONT")       # case-insensitive
  out <- str_replace_all(out, "(?i)pacbio", "PacBio")
  # 3) Capitalize Short read
  out <- str_replace(out, "short read", "Short Read")
  return(out)
}

original_names <- names(fn_list)
nice_names <- sapply(original_names, pretty_label, USE.NAMES=FALSE)

# rename in the list, then create the UpSet data
names(fn_list) <- nice_names


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8) Create an UpSet plot using ComplexUpset & improve the look
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Convert the list of false negatives to a binary incidence data.frame
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

# Define the order of sets (reverse order so "Short read" is last)
sets_order <- rev(names(fn_list))

# Create the UpSet plot using ComplexUpset & improve the look
if (has_ComplexUpset && nrow(upset_data) > 0) {
  p <- ComplexUpset::upset(
    upset_data,
    width_ratio = 0.2,
    height_ratio = 0.7,
    intersect = sets_order,
    keep_empty_groups = TRUE,
    wrap = TRUE,
    set_sizes = FALSE,
    base_annotations = list(
      'Intersect\nFN TUSCO\nGenes' = ComplexUpset::intersection_size(counts = TRUE)
    )
  ) +
    labs(title = NULL, subtitle = NULL, caption = NULL, y = "TUSCO\n(Human)") +
    theme_classic(base_family = "Helvetica", base_size = 7) +
    theme(
      text = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.title.y = element_text(size = 7, angle = 90, vjust = 0.5),
      strip.text.x = element_text(size = 7),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "none"
    )
} else {
  p <- NULL
  message("ComplexUpset not available or empty data; skipping UpSet plot.")
}

# Save the plot to PDF
if (FALSE) ggsave(
  filename = file.path(output_dir, "figure3d-human.pdf"), # Disabled per request
  plot = p,
  width = 4.5, # Adjusted width slightly
  height = 3.5, # Adjusted height slightly
  units = "in",
  device = cairo_pdf, # Use cairo_pdf for better font embedding
  dpi = 300
)

message("UpSet plot step completed (plot may be skipped).")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9) Print out the intersection of all sets (common FNs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
common_fn <- Reduce(dplyr::intersect, fn_list)
message("Number of common FN TUSCO genes across all sets: ", length(common_fn))

# ---- START DEBUG ----
message("DEBUG: Content of common_fn below:")
if (length(common_fn) > 0) {
  print(common_fn)
} else {
  message("DEBUG: common_fn is empty.")
}
message("DEBUG: Structure of common_fn below:")
print(str(common_fn))
# ---- END DEBUG ----

# Do not write extra files outside the consolidated TSV export

message("Done.")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 10) Plot distribution of FNs by number of samples shared
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine all FN genes from all samples into a single vector
all_fn_genes <- unlist(fn_list)

# Count how many times each gene appears (i.e., in how many samples it's an FN)
gene_counts <- table(all_fn_genes)

# Count how many genes appear in exactly 1 sample, 2 samples, etc.
sample_counts <- table(gene_counts)

# Convert to a data frame for plotting
sample_counts_df <- as.data.frame(sample_counts)
colnames(sample_counts_df) <- c("num_samples_shared", "num_genes")

# Ensure num_samples_shared is numeric (factor by default from table)
sample_counts_df$num_samples_shared <- as.numeric(as.character(sample_counts_df$num_samples_shared))

# Ensure all sample counts from 1 to N are present, even if count is 0
num_total_samples <- length(fn_list)
all_possible_counts <- data.frame(num_samples_shared = 1:num_total_samples)

plot_df <- merge(all_possible_counts, sample_counts_df, by = "num_samples_shared", all.x = TRUE)
plot_df$num_genes[is.na(plot_df$num_genes)] <- 0 # Replace NA with 0 for missing counts


# Create the bar plot
p_dist <- ggplot(plot_df, aes(x = factor(num_samples_shared), y = num_genes)) +
  geom_bar(stat = "identity", fill = "#A8D5A0") + # Changed fill color to user specified hex code
  labs(
    title = NULL, # Removed title
    x = "#Sample", # Updated x-axis label
    y = "FN TUSCO Genes" # Updated y-axis label, Changed from BUGSI
  ) +
  # Add padding to the y-axis to ensure labels fit
  # ylim(0, max(plot_df$num_genes, na.rm = TRUE) * 1.15) + # Removed dynamic ylim
  # Set fixed y-axis scale and breaks
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15)) +
  theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5), # Matched margins
    # Remove background elements
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Ensure axis lines are still present if desired (optional, theme_minimal keeps them)
    # axis.line = element_line(colour = "black")
  ) +
  # Add text labels on top of bars
  geom_text(aes(label = num_genes), vjust = -0.5, size = 2) # Adjusted text size and position

# Save the plot
dist_plot_file <- file.path(plot_dir, "figure3d-human.pdf") # Output only under ./plot
ggsave(
  filename = dist_plot_file,
  plot = p_dist,
  width = 2, # Matched width
  height = 1.2, # Matched height
  units = "in",
  device = cairo_pdf, # Matched device
  dpi = 300
)

message("FN TUSCO distribution bar plot saved to: ", dist_plot_file)

# Export TSV with underlying data and minimal metadata
tsv_path <- file.path(tsv_dir, "figure3d-human.tsv")
plot_df_out <- plot_df %>% mutate(figure_id = "fig-3", panel_id = "3d-human", record_type = "distribution_counts")
readr::write_tsv(plot_df_out, tsv_path)
