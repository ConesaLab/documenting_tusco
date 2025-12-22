#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# --------------------------------------------------------------------------------------
# Resolve paths robustly based on this script's location
# --------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0) {
  # Fallback when run interactively; use current working directory
  script_dir <- normalizePath(getwd())
} else {
  script_dir <- normalizePath(dirname(script_path))
}
fig_dir  <- normalizePath(file.path(script_dir, ".."))
repo_root <- normalizePath(file.path(fig_dir, "../.."))

# Inputs and outputs
csv_file <- file.path(repo_root, "data","raw","lrgasp","FN_correlation_plot.csv")
human_parent_dir <- file.path(repo_root, "data","raw","lrgasp","human")
mouse_parent_dir <- file.path(repo_root, "data","raw","lrgasp","mouse")
tusco_human_file <- file.path(repo_root, "data","processed","tusco","hsa","tusco_human.tsv")
tusco_mouse_file <- file.path(repo_root, "data","processed","tusco","mmu","tusco_mouse.tsv")

plot_dir <- file.path(fig_dir, "plots")
tsv_dir  <- file.path(fig_dir, "tables")
run_log  <- file.path(fig_dir, "run.log")

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(..., collapse = ""), "\n")
  cat(msg, file = run_log, append = TRUE)
  message(paste(..., collapse = ""))
}

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tsv_dir, showWarnings = FALSE, recursive = TRUE)

log_message("[start] Figure S3 run in ", fig_dir)

#############################################
# Load CSV (metrics) and compute FN ourselves
#############################################
if (!file.exists(csv_file)) {
  log_message("[error] Missing input CSV: ", csv_file, ". Skipping plot generation.")
  quit(status = 0)
}

log_message("[info] Reading CSV: ", csv_file)
metrics_df <- tryCatch({
  read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8") %>%
    mutate(
      Sample = trimws(Sample),
      `Log(total Reads Number) * Median Read Length` = suppressWarnings(as.numeric(as.character(`Log(total Reads Number) * Median Read Length`))),
      `Total Reads Number` = suppressWarnings(as.numeric(as.character(`Total Reads Number`))),
      `Median Read Length` = suppressWarnings(as.numeric(as.character(`Median Read Length`)))
    )
}, error = function(e) {
  log_message("[error] Failed to read CSV: ", conditionMessage(e))
  NULL
})

if (is.null(metrics_df)) quit(status = 0)

required_cols <- c("Sample", "Log(total Reads Number) * Median Read Length")
missing_cols <- setdiff(required_cols, colnames(metrics_df))
if (length(missing_cols) > 0) {
  log_message("[error] CSV missing required columns: ", paste(missing_cols, collapse = ", "), ". Skipping plot generation.")
  quit(status = 0)
}

# --------------------------------------------------------------------------------------
# Helpers to read TUSCO and classification files
# --------------------------------------------------------------------------------------
read_tusco_genes <- function(tsv_path) {
  lines <- readr::read_lines(tsv_path, progress = FALSE)
  if (length(lines) == 0) return(character())
  first_tokens <- stringr::str_trim(stringr::str_extract(lines, "^[^\\s#]+"))
  ids <- first_tokens[!is.na(first_tokens) & first_tokens != ""]
  ids <- stringr::str_remove(ids, "\\..+$")
  unique(ids)
}

get_tusco_fsm_ism_found <- function(class_file, tusco_set) {
  if (!file.exists(class_file)) {
    stop("Classification file not found: ", class_file)
  }
  df <- read_tsv(class_file, col_types = cols(.default = "c"))
  colnames(df)[1:7] <- c("isoform","chrom","strand","length","exons","structural_category","associated_gene")
  df_clean <- df %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      df %>%
        filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(associated_gene = str_remove(associated_gene, "\\.\\d+$"))
  found_genes <- df_clean %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") |
           subcategory == "mono-exon_by_intron_retention") %>%
    distinct(associated_gene) %>%
    pull(associated_gene)
  intersect(tusco_set, found_genes)
}

map_sample_to_classification <- function(sample_name) {
  parts <- str_split(sample_name, "-", simplify = TRUE)
  if (ncol(parts) < 3) {
    stop("Unexpected Sample format: ", sample_name)
  }
  cell <- parts[1]
  method <- tolower(parts[2])  # cdna or drna
  platform <- tolower(parts[3]) # ont or pacbio
  is_ls <- grepl("-LS$", sample_name)

  folder_prefix <- paste0(cell, "_", method, "_")
  folder <- paste0(folder_prefix, platform)
  if (is_ls) folder <- paste0(folder, "_ls")

  parent_dir <- if (cell == "WTC11") human_parent_dir else if (cell == "ES") mouse_parent_dir else stop("Unknown cell prefix: ", cell)

  standard_path <- file.path(parent_dir, folder, paste0(folder, "_classification.txt"))
  sqanti_path <- file.path(parent_dir, folder, paste0(folder, "_sqanti_classification.txt"))

  class_file <- if (file.exists(standard_path)) standard_path else if (file.exists(sqanti_path)) sqanti_path else NA_character_

  list(cell = cell, folder = folder, class_file = class_file, parent_dir = parent_dir)
}

# Pre-load TUSCO gene sets
if (!file.exists(tusco_human_file)) {
  log_message("[error] Missing TUSCO human TSV: ", tusco_human_file)
  quit(status = 0)
}
if (!file.exists(tusco_mouse_file)) {
  log_message("[error] Missing TUSCO mouse TSV: ", tusco_mouse_file)
  quit(status = 0)
}

tusco_genes_human <- read_tusco_genes(tusco_human_file)
tusco_genes_mouse <- read_tusco_genes(tusco_mouse_file)

# For each sample, compute FN percentage
compute_fn_percent_for_sample <- function(sample_name) {
  map_info <- map_sample_to_classification(sample_name)
  class_file <- map_info$class_file
  if (is.na(class_file)) {
    log_message("[warn] Classification file not found for sample ", sample_name, "; returning NA")
    return(NA_real_)
  }
  tusco_set <- if (map_info$cell == "WTC11") tusco_genes_human else tusco_genes_mouse
  found <- get_tusco_fsm_ism_found(class_file, tusco_set)
  total_tusco <- length(tusco_set)
  if (total_tusco == 0) return(NA_real_)
  fn_count <- length(setdiff(tusco_set, found))
  (fn_count / total_tusco) * 100
}

metrics_df <- metrics_df %>%
  rowwise() %>%
  mutate(`False Negatives (%)` = compute_fn_percent_for_sample(Sample)) %>%
  ungroup() %>%
  filter(!is.na(`False Negatives (%)`) & !is.na(`Log(total Reads Number) * Median Read Length`))

if (nrow(metrics_df) < 2) {
  log_message("[error] Not enough data points to plot (n=", nrow(metrics_df), "). Skipping plot generation.")
  quit(status = 0)
}

# --------------------------------------------------------------------------------------
# Correlation and plotting
# --------------------------------------------------------------------------------------
plot_df <- metrics_df

cor_test <- suppressWarnings(tryCatch({
  cor.test(
    plot_df$`Log(total Reads Number) * Median Read Length`,
    plot_df$`False Negatives (%)`,
    method = "pearson"
  )
}, error = function(e) NULL))

r_value <- if (!is.null(cor_test)) as.numeric(cor_test$estimate) else NA_real_
p_value <- if (!is.null(cor_test)) cor_test$p.value else NA_real_
n_points <- nrow(plot_df)

custom_theme <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.text  = element_text(size = rel(1), color = "black"),
    axis.title = element_text(size = rel(1.1)),
    axis.line  = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.title = element_text(size=rel(1.0), face="bold"),
    legend.text = element_text(size=rel(0.9))
  )

# Derive grouping aesthetics
plot_df <- plot_df %>%
  mutate(
    cell_line = str_extract(Sample, "^(ES|WTC11)"),
    method_simple = case_when(
      grepl("PacBio", Sample, fixed = TRUE) ~ "PacBio",
      grepl("ONT", Sample, fixed = TRUE) ~ "ONT",
      TRUE ~ "Other"
    ),
    method = case_when(
      method_simple == "PacBio" & grepl("cDNA", Sample, fixed=TRUE) ~ "PacBio-cDNA",
      method_simple == "ONT" & grepl("cDNA", Sample, fixed=TRUE) ~ "ONT-cDNA",
      method_simple == "ONT" & grepl("dRNA", Sample, fixed=TRUE) ~ "ONT-dRNA",
      TRUE ~ method_simple
    ),
    is_ls = grepl("-LS$", Sample),
    point_size = if_else(is_ls, 4.5, 3)
  )

# Select a minimal informative subset of labels to avoid stacking
label_select_n <- 3
x_order <- order(plot_df$`Log(total Reads Number) * Median Read Length`, decreasing = TRUE)
y_order <- order(plot_df$`False Negatives (%)`, decreasing = TRUE)
label_candidates <- unique(c(
  plot_df$Sample[head(x_order, label_select_n)],
  plot_df$Sample[tail(x_order, label_select_n)],
  plot_df$Sample[head(y_order, label_select_n)],
  plot_df$Sample[tail(y_order, label_select_n)],
  plot_df$Sample[plot_df$is_ls]
))
# Always label these key samples if present
always_label_samples <- c("ES-cDNA-ONT")
label_candidates <- unique(c(label_candidates, plot_df$Sample[plot_df$Sample %in% always_label_samples]))
plot_df <- plot_df %>% mutate(label_me = Sample %in% label_candidates)

subtitle_text <- if (!is.na(r_value) && !is.na(p_value)) {
  paste0("Pearson r = ", round(r_value, 2), ", p = ", format(p_value, scientific = TRUE, digits = 2))
} else {
  "Pearson correlation not available"
}

p <- ggplot(plot_df, aes(x = `Log(total Reads Number) * Median Read Length`, y = `False Negatives (%)`)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgray", fill = "gray", linetype = "solid", alpha = 0.2, linewidth = 0.4) +
  geom_point(aes(color = cell_line, shape = method, size = point_size), alpha = 0.9, stroke = 0.5) +
  {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(
        data = subset(plot_df, label_me),
        aes(label = Sample),
        color = "black",
        size = 2.5,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.6,
        force = 3,
        direction = "both",
        seed = 123,
        max.time = 2,
        segment.color = 'grey50',
        segment.size = 0.2,
        segment.alpha = 0.6,
        min.segment.length = 0
      )
    } else {
      geom_text(
        data = subset(plot_df, label_me),
        aes(label = Sample),
        size = 2.5,
        vjust = -0.6,
        color = "black",
        check_overlap = TRUE
      )
    }
  } +
  scale_color_manual(name = "Cell Line", values = c("WTC11" = "#a6dba0", "ES" = "#1C9E77")) +
  scale_shape_manual(name = "Method", values = c("PacBio-cDNA" = 16, "ONT-cDNA" = 17, "ONT-dRNA" = 18)) +
  scale_size_identity(guide = "legend", name = "Short-read Support", breaks = c(3, 4.5), labels = c("No", "Yes")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  guides(color = guide_legend(title="Cell Line", override.aes = list(size=4)),
         shape = guide_legend(title="Method", override.aes = list(size=4)),
         size = guide_legend(title="Short-read Support")) +
  labs(x = "Log(Total Reads) * Median Read Length",
       y = "False Negatives (%)",
       subtitle = subtitle_text) +
  custom_theme +
  coord_cartesian(clip = "off") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(10, 10, 10, 10))

outfile_pdf <- file.path(plot_dir, "fig-s3.pdf")
ggsave(filename = outfile_pdf, plot = p, width = 5, height = 3.5, device = "pdf")
log_message("[ok] Wrote PDF: ", outfile_pdf)

# --------------------------------------------------------------------------------------
# Write TSV with underlying data and metadata
# --------------------------------------------------------------------------------------
outfile_tsv <- file.path(tsv_dir, "fig-s3.tsv")
tsv_df <- plot_df %>%
  transmute(
    figure_id = "fig-s3",
    panel_id = NA_character_,
    sample_id = Sample,
    x_name = "Log(total Reads Number) * Median Read Length",
    y_name = "False Negatives (%)",
    x_value = `Log(total Reads Number) * Median Read Length`,
    y_value = `False Negatives (%)`,
    cell_line = cell_line,
    method = method,
    short_read_support = if_else(is_ls, "Yes", "No"),
    n = n_points,
    pearson_r = r_value,
    p_value = p_value
  )
readr::write_tsv(tsv_df, outfile_tsv)
log_message("[ok] Wrote TSV: ", outfile_tsv)

log_message("[done] Figure S3 completed. n=", n_points, ", r=", round(r_value, 4), ", p=", format(p_value, scientific = TRUE, digits = 3))
