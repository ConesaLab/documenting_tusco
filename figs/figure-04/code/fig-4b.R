###############################################################
# Script for TUSCO Evaluation (Figure 4)
# Generates a single overall 8x5 radar grid PDF: fig-4b.pdf
###############################################################

message("Loading necessary libraries...")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(fmsb)
  library(cowplot)
  library(grid)
  library(ggplotify)
})
# Helper to draw a group header with a bracket line below the title
make_group_header <- function(title_text, title_size = 9) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.75, label = title_text, family = "Helvetica", fontface = "bold", size = title_size/3) +
    annotate("segment", x = 0.05, xend = 0.95, y = 0.12, yend = 0.12, linewidth = 0.35) +
    annotate("segment", x = 0.05, xend = 0.05, y = 0.12, yend = 0.02, linewidth = 0.35) +
    annotate("segment", x = 0.95, xend = 0.95, y = 0.12, yend = 0.02, linewidth = 0.35) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
}

# ------------------------------------------------------------------
# Paths (updated to NIH data and repo plots directory)
# ------------------------------------------------------------------
# Localize all paths to this fig-4 folder only
argv <- commandArgs(trailingOnly = FALSE)
script_path <- tryCatch({
  sub("^--file=", "", argv[grep("^--file=", argv)][1])
}, error = function(e) NA_character_)
if (is.na(script_path) || script_path == "") {
  # Fallback: assume working directory contains figure-04/code
  script_path <- file.path(getwd(), "figs", "figure-04", "code", "fig-4b.R")
}
fig_dir  <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
figs_dir <- normalizePath(file.path(fig_dir, ".."), winslash = "/", mustWork = FALSE) # repository figs/
repo_dir <- normalizePath(file.path(figs_dir, ".."), winslash = "/", mustWork = FALSE)

first_existing <- function(paths) {
  for (p in paths) {
    if (!is.na(p) && !is.null(p) && file.exists(p)) return(p)
  }
  # return first candidate even if missing (to keep consistent path intent)
  if (length(paths) > 0) return(paths[[1]]) else return(NA_character_)
}

first_existing_dir <- function(paths) {
  for (p in paths) {
    if (!is.na(p) && !is.null(p) && dir.exists(p)) return(p)
  }
  if (length(paths) > 0) return(paths[[1]]) else return(NA_character_)
}

# TUSCO references: prefer local figure-04/data, fallback to processed data under repo
tusco_human_file <- first_existing(list(
  file.path(fig_dir,  "data", "tusco_human_multi_exon.tsv"),
  file.path(repo_dir, "data", "processed", "tusco", "hsa", "tusco_human_multi_exon.tsv")
))
tusco_mouse_file <- first_existing(list(
  file.path(fig_dir,  "data", "tusco_mouse_multi_exon.tsv"),
  file.path(repo_dir, "data", "processed", "tusco", "mmu", "tusco_mouse_multi_exon.tsv")
))

# Base data directory: prefer local figure-04/data, fallback to repo data/raw/lrgasp/tusco_novel_evl
base_data_dir <- first_existing_dir(list(
  file.path(fig_dir,  "data", "lrgasp", "tusco_novel_evl"),
  file.path(repo_dir, "data", "raw", "lrgasp", "tusco_novel_evl")
))

# Output directories (restricted to fig-4)
output_dir <- file.path(fig_dir, "plots")
tsv_dir    <- file.path(fig_dir, "tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(tsv_dir)) dir.create(tsv_dir, recursive = TRUE)

# ------------------------------------------------------------------
# Helpers retained (minimal set required to render the overall grid)
# ------------------------------------------------------------------

# Read TSV safely
read_tsv_safe <- function(file_path, col_names = TRUE, ...) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  tryCatch(
    {
      data <- read_tsv(file_path, col_names = col_names, show_col_types = FALSE, ...)
      return(data)
    },
    error = function(e) {
      warning("Error reading file ", file_path, ": ", e$message)
      return(NULL)
    }
  )
}

# Compute TUSCO metrics (Sn, nrPre, 1/red, 1-FDR, PDR, rPre) in percent
calculate_tusco_metrics <- function(classification_data, tusco_annotation_df, rTUSCO_total_ref_transcripts) {
  if (is.null(classification_data) || nrow(classification_data) == 0) {
    return(rep(0, 6))
  }

  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)[0-9]{11}(\\.[0-9]+)?$",
    refseq    = "^(NM_|NR_|NP_)[0-9]{6,}(\\.[0-9]+)?$",
    gene_name = "^[A-Za-z0-9][A-Za-z0-9.-]*[A-Za-z0-9]$"
  )

  classification_data_cleaned <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.[0-9]+$"),
      id_type = case_when(
        str_detect(associated_gene, patterns$ensembl) ~ "ensembl",
        str_detect(associated_gene, patterns$refseq)  ~ "refseq",
        associated_gene %in% tusco_annotation_df$gene_name ~ "gene_name",
        TRUE ~ "unknown"
      )
    ) %>%
    filter(id_type != "unknown") %>%
    distinct(isoform, associated_gene, .keep_all = TRUE)

  if (nrow(classification_data_cleaned) == 0) {
    return(rep(0, 6))
  }

  id_summary <- classification_data_cleaned %>% count(id_type, sort = TRUE)
  # Choose the majority (most frequent) non-unknown id_type; fallback to gene_name
  top_id_type <- id_summary %>%
    dplyr::filter(id_type != "unknown") %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::pull(id_type)
  if (length(top_id_type) == 0 || is.na(top_id_type)) top_id_type <- "gene_name"

  tusco_transcripts <- classification_data_cleaned %>%
    filter(associated_gene %in% tusco_annotation_df[[top_id_type]])

  TUSCO_RM <- tusco_transcripts %>%
    filter(
      # TP rule (uniform): RM or FSM mono-exon with both ends within 50bp
      subcategory == "reference_match" |
      (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
         !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
         abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50) |
      # Rule 3: For ref_length > 3000, FSM with TSS and TTS within 100bp is TP
      (!is.na(ref_length) & ref_length > 3000 & 
       structural_category == "full-splice_match" &
       !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
       abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100)
    )
  TP_tusco <- TUSCO_RM
  TP_TUSCO_unique_genes <- unique(TP_tusco$associated_gene)

  PTP_tusco <- tusco_transcripts %>%
    filter((structural_category %in% c("full-splice_match", "incomplete-splice_match") |
            subcategory == "mono-exon_by_intron_retention") &
           !isoform %in% TP_tusco$isoform)

  detected_tusco_genes <- unique(tusco_transcripts$associated_gene)
  FN_tusco <- tusco_annotation_df %>%
    filter(!(!!sym(top_id_type) %in% detected_tusco_genes))

  FP_tusco <- tusco_transcripts %>%
    filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog", "genic", "fusion", "antisense", "intergenic", "genic_intron") &
           subcategory != "mono-exon_by_intron_retention")

  fsm_ism_count_tusco <- tusco_transcripts %>%
    filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    nrow()

  non_redundant_sensitivity_tusco <- if (rTUSCO_total_ref_transcripts > 0) length(TP_TUSCO_unique_genes) / rTUSCO_total_ref_transcripts else 0
  # Make nrPre transcript-based: TP transcripts / total detected transcripts
  non_redundant_precision_tusco  <- if (nrow(tusco_transcripts) > 0) nrow(TP_tusco) / nrow(tusco_transcripts) else 0

  unique_tp_ptp_tusco_genes <- length(unique(c(TP_TUSCO_unique_genes, PTP_tusco$associated_gene)))
  redundancy_tusco <- if (unique_tp_ptp_tusco_genes > 0) fsm_ism_count_tusco / unique_tp_ptp_tusco_genes else 0
  inv_redundancy_metric <- if (redundancy_tusco > 0) (1 / redundancy_tusco) else 0

  false_discovery_rate_tusco <- if (nrow(tusco_transcripts) > 0) (nrow(tusco_transcripts) - nrow(TUSCO_RM)) / nrow(tusco_transcripts) else 0
  one_minus_fdr_metric <- 1 - false_discovery_rate_tusco

  positive_detection_rate_tusco <- if (rTUSCO_total_ref_transcripts > 0) unique_tp_ptp_tusco_genes / rTUSCO_total_ref_transcripts else 0
  redundant_precision_tusco <- if (nrow(tusco_transcripts) > 0) (nrow(TP_tusco) + nrow(PTP_tusco)) / nrow(tusco_transcripts) else 0

  metrics_values <- c(
    Sn    = non_redundant_sensitivity_tusco * 100,
    nrPre = non_redundant_precision_tusco * 100,
    `1/red` = inv_redundancy_metric * 100,
    `1-FDR` = one_minus_fdr_metric * 100,
    PDR   = positive_detection_rate_tusco * 100,
    rPre  = redundant_precision_tusco * 100
  )
  metrics_values[is.na(metrics_values)] <- 0
  return(metrics_values)
}

# Render a single radar chart grob for a given pipeline/sample/species
metrics_records <- list()
area_ratio_records <- list()

# helper: append a metrics record safely
append_record <- function(rec_df) {
  idx <- length(metrics_records) + 1L
  metrics_records[[idx]] <<- rec_df
}

append_area_record <- function(rec_df) {
  idx <- length(area_ratio_records) + 1L
  area_ratio_records[[idx]] <<- rec_df
}

compute_radar_area <- function(metric_vec, max_value = 100) {
  if (is.null(metric_vec)) return(NA_real_)
  values <- as.numeric(metric_vec)
  if (length(values) < 3) return(NA_real_)
  if (anyNA(values)) return(NA_real_)
  scaled <- values / max_value
  angles <- seq(0, 2 * pi, length.out = length(values) + 1)
  angles <- angles[-length(angles)]
  x <- scaled * cos(angles)
  y <- scaled * sin(angles)
  area <- 0.5 * abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])))
  area
}

compute_area_performance_drop <- function(vec_ref, vec_novel, max_value = 100) {
  area_ref <- compute_radar_area(vec_ref, max_value)
  area_novel <- compute_radar_area(vec_novel, max_value)
  drop_pct <- if (!is.na(area_ref) && area_ref > 0 && !is.na(area_novel)) {
    (area_ref - area_novel) / area_ref * 100
  } else {
    NA_real_
  }
  list(
    area_reference = area_ref,
    area_novel = area_novel,
    performance_drop_pct = drop_pct
  )
}

sanitize_names <- function(x) {
  x <- gsub("/", "_", x, fixed = TRUE)
  x <- gsub("-", "_", x, fixed = TRUE)
  x
}

process_and_plot_pipeline_pair <- function(pipeline_prefix, species, tusco_ref_filepath, data_dir_base,
                                           is_ml_processing = FALSE, omit_title = TRUE,
                                           pipeline_display = NULL, pipeline_identifier = NULL) {
  main_pipeline_type <- basename(data_dir_base)
  pipeline_label <- if (!is.null(pipeline_display)) pipeline_display else main_pipeline_type
  pipeline_id <- if (!is.null(pipeline_identifier)) pipeline_identifier else main_pipeline_type

  # Choose SQANTI directory and filename patterns per pipeline
  if (main_pipeline_type == "ujc_sq3") {
    sqanti_dir_name <- "sqanti3"
    ref_class_filename <- paste0(pipeline_prefix, "_sqanti3_ref_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_sqanti3_ref_classification.txt")
  } else if (main_pipeline_type == "stringtie_sq3") {
    sqanti_dir_name <- "sqanti3_out"
    ref_class_filename <- paste0(pipeline_prefix, "_stringtie_ref_sqanti_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_stringtie_tusco_sqanti_classification.txt")
  } else if (main_pipeline_type == "all-ref_sq3") {
    sqanti_dir_name <- "sqanti3"
    ref_class_filename <- paste0(pipeline_prefix, "_sqanti3_all_ref_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_sqanti3_all_ref_classification.txt")
  } else if (main_pipeline_type == "isoseq_sq3") {
    sqanti_dir_name <- "sqanti3_out"
    ref_class_filename <- "isoforms_classification.txt"
    novel_class_filename <- "isoforms_classification.txt"
  } else {
    sqanti_dir_name <- "sqanti3_out"
    ref_class_filename <- paste0(pipeline_prefix, "_sqanti_classification.txt")
    novel_class_filename <- paste0(pipeline_prefix, "_sqanti_classification.txt")
  }

  ref_evl_path <- file.path(data_dir_base, "ref_evl", pipeline_prefix, sqanti_dir_name)
  novel_evl_path <- file.path(data_dir_base, "novel_evl", pipeline_prefix, sqanti_dir_name)

  ref_class_file <- if (dir.exists(ref_evl_path)) file.path(ref_evl_path, ref_class_filename) else NULL
  novel_class_file <- if (dir.exists(novel_evl_path)) file.path(novel_evl_path, novel_class_filename) else NULL

  # Fallback: dynamically search for isoforms_classification.txt within the sample folder
  if (is.null(ref_class_file) || !file.exists(ref_class_file)) {
    sample_ref_dir <- file.path(data_dir_base, "ref_evl", pipeline_prefix)
    if (dir.exists(sample_ref_dir)) {
      hits <- list.files(sample_ref_dir, pattern = "^isoforms_classification\\.txt$", recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) ref_class_file <- hits[[1]]
    }
  }
  if (is.null(novel_class_file) || !file.exists(novel_class_file)) {
    sample_nov_dir <- file.path(data_dir_base, "novel_evl", pipeline_prefix)
    if (dir.exists(sample_nov_dir)) {
      hits <- list.files(sample_nov_dir, pattern = "^isoforms_classification\\.txt$", recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) novel_class_file <- hits[[1]]
    }
  }

  tusco_cols <- c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq")
  tusco_df <- NULL
  if (file.exists(tusco_ref_filepath)) {
    tusco_df <- tryCatch({
      read_delim(tusco_ref_filepath, delim = "\t", col_names = tusco_cols, col_types = cols(.default = "c"), trim_ws = TRUE, comment = "#", show_col_types = FALSE, guess_max = 2000)
    }, error = function(e) NULL)
    if (!is.null(tusco_df) && ncol(tusco_df) == 1) {
      tusco_df <- tryCatch({
        read_table2(tusco_ref_filepath, col_names = tusco_cols, col_types = cols(.default = "c"), comment = "#", show_col_types = FALSE, guess_max = 2000)
      }, error = function(e) NULL)
    }
  }
  if (is.null(tusco_df) || nrow(tusco_df) == 0) {
    append_record(data.frame(
      figure_id = "fig-4b",
      pipeline  = main_pipeline_type,
      pipeline_id = pipeline_id,
      pipeline_label = pipeline_label,
      sample    = pipeline_prefix,
      species   = species,
      eval_type = NA_character_,
      status    = "no_tusco_ref",
      stringsAsFactors = FALSE
    ))
    return(ggdraw())
  }

  tusco_annotation_for_calc <- tusco_df %>% select(any_of(c("ensembl", "refseq", "gene_name"))) %>% distinct()
  rTUSCO <- nrow(tusco_annotation_for_calc)
  if (rTUSCO == 0) {
    return(ggdraw() + draw_label("Empty TUSCO ref", fontface = "bold"))
  }

  ref_classification_data <- if (!is.null(ref_class_file) && file.exists(ref_class_file)) read_tsv_safe(ref_class_file) else NULL
  novel_classification_data <- if (!is.null(novel_class_file) && file.exists(novel_class_file)) read_tsv_safe(novel_class_file) else NULL

  metrics_ref <- if (!is.null(ref_classification_data)) calculate_tusco_metrics(ref_classification_data, tusco_annotation_for_calc, rTUSCO) else NULL
  metrics_novel <- if (!is.null(novel_classification_data)) calculate_tusco_metrics(novel_classification_data, tusco_annotation_for_calc, rTUSCO) else NULL

  if (is.null(metrics_ref) && is.null(metrics_novel)) {
    return(ggdraw())
  }

  if (!is.null(metrics_ref) && !is.null(metrics_novel)) {
    radar_data_df <- rbind(rep(100, 6), rep(0, 6), metrics_ref, metrics_novel)
    plot_colors <- rep(if (species == "human") "#a8d5a0" else "#1b9e77", 2)
    plot_linetypes <- c(1, 2)
    # log metrics
    mref <- as.list(metrics_ref)
    names(mref) <- sanitize_names(names(mref))
    mnov <- as.list(metrics_novel)
    names(mnov) <- sanitize_names(names(mnov))
    append_record(as.data.frame(c(
      list(figure_id = "fig-4b", pipeline = main_pipeline_type, pipeline_id = pipeline_id, pipeline_label = pipeline_label,
           sample = pipeline_prefix, species = species, eval_type = "ref", status = "ok"),
      mref
    ), check.names = FALSE, stringsAsFactors = FALSE))
    append_record(as.data.frame(c(
      list(figure_id = "fig-4b", pipeline = main_pipeline_type, pipeline_id = pipeline_id, pipeline_label = pipeline_label,
           sample = pipeline_prefix, species = species, eval_type = "novel", status = "ok"),
      mnov
    ), check.names = FALSE, stringsAsFactors = FALSE))
    area_metrics <- compute_area_performance_drop(metrics_ref, metrics_novel)
    append_area_record(tibble(
      figure_id = "fig-4b",
      pipeline = main_pipeline_type,
      pipeline_id = pipeline_id,
      pipeline_label = pipeline_label,
      sample = pipeline_prefix,
      species = species,
      area_reference = area_metrics$area_reference,
      area_novel = area_metrics$area_novel,
      performance_drop_pct = area_metrics$performance_drop_pct
    ))
  } else if (!is.null(metrics_ref)) {
    radar_data_df <- rbind(rep(100, 6), rep(0, 6), metrics_ref)
    plot_colors <- if (species == "human") "#a8d5a0" else "#1b9e77"
    plot_linetypes <- 1
    mref <- as.list(metrics_ref)
    names(mref) <- sanitize_names(names(mref))
    append_record(as.data.frame(c(
      list(figure_id = "fig-4b", pipeline = main_pipeline_type, pipeline_id = pipeline_id, pipeline_label = pipeline_label,
           sample = pipeline_prefix, species = species, eval_type = "ref", status = "ok"),
      mref
    ), check.names = FALSE, stringsAsFactors = FALSE))
  } else {
    radar_data_df <- rbind(rep(100, 6), rep(0, 6), metrics_novel)
    plot_colors <- if (species == "human") "#a8d5a0" else "#1b9e77"
    plot_linetypes <- 2
    mnov <- as.list(metrics_novel)
    names(mnov) <- sanitize_names(names(mnov))
    append_record(as.data.frame(c(
      list(figure_id = "fig-4b", pipeline = main_pipeline_type, pipeline_id = pipeline_id, pipeline_label = pipeline_label,
           sample = pipeline_prefix, species = species, eval_type = "novel", status = "ok"),
      mnov
    ), check.names = FALSE, stringsAsFactors = FALSE))
  }

  colnames(radar_data_df) <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
  radar_data_df <- as.data.frame(radar_data_df)

  # Capture radarchart as a grob without writing temp files
  plot_obj <- NULL
  tryCatch({
    plot_obj <- ggplotify::as.ggplot(function() {
      op <- par(mar = c(0.2, 0.2, if (omit_title) 0.2 else 1.0, 0.2), family = "Helvetica", font = 2, cex = 0.5)
      on.exit(par(op), add = TRUE)
      radarchart(radar_data_df,
                 axistype = 0,
                 pcol = plot_colors,
                 pfcol = NA,
                 plwd = 1.2,
                 plty = plot_linetypes,
                 pty = 16,
                 cglcol = "grey", cglty = 1, cglwd = 0.35,
                 axislabcol = "black",
                 vlabels = rep("", ncol(radar_data_df)),
                 vlcex = 0,
                 caxislabels = NULL,
                 centerzero = FALSE,
                 title = if (omit_title) "" else pipeline_prefix)
    }) +
      coord_fixed(ratio = 1) +
      theme_void() +
      theme(plot.margin = unit(rep(0.5, 4), "mm"))
  }, error = function(e) {
    plot_obj <<- ggdraw() + draw_label("Plot Error", fontface = "bold")
  })
  return(plot_obj)
}

# ------------------------------------------------------------------
# Build legend used underneath the grid
# ------------------------------------------------------------------
legend_line_df <- data.frame(
  Annotation = factor(c("Gencode reference annotation", "TUSCO-sim annotation"),
                      levels = c("Gencode reference annotation", "TUSCO-sim annotation")),
  x = 0, xend = 1, y = c(1, 2), yend = c(1, 2)
)
p_legend_lines <- ggplot(legend_line_df) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend, linetype = Annotation),
               linewidth = 1.2, color = "black") +
  scale_linetype_manual(values = c("Gencode reference annotation" = "solid",
                                   "TUSCO-sim annotation" = "dashed")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 7, family = "Helvetica"))
legend_linetype <- cowplot::get_legend(p_legend_lines)

legend_species_df <- data.frame(
  Species = factor(c("Human", "Mouse"), levels = c("Human", "Mouse")),
  x = 1, y = c(1, 2)
)
p_legend_species <- ggplot(legend_species_df, aes(x = x, y = y, color = Species)) +
  geom_point(shape = 15, size = 5) +
  scale_color_manual(values = c("Human" = "#a8d5a0", "Mouse" = "#1b9e77")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 7, family = "Helvetica"))
legend_species <- cowplot::get_legend(p_legend_species)

radar_plot_legend <- plot_grid(legend_linetype, legend_species, ncol = 2, rel_widths = c(2, 1))

# ------------------------------------------------------------------
# Generate the overall 8x5 radar grid and save as figure4.pdf
# ------------------------------------------------------------------
message("Generating overall 8x5 radar grid figure...")

# Ordered sample prefixes (columns)
overall_sample_prefixes <- c(
  "wtc11_captrap_pacbio",
  "wtc11_cdna_pacbio",
  "wtc11_cdna_ont",
  "wtc11_captrap_ont",
  "es_captrap_pacbio",
  "es_cdna_pacbio",
  "es_cdna_ont",
  "es_captrap_ont"
)

# Pipeline rows (updated labels, removed All-ref)
pipeline_row_specs <- list(
  list(row_label = "Bambu",            pipeline_dir = "bambu_sq3",    pipeline_id = "bambu",              is_ml = FALSE),
  list(row_label = "StringTie2",       pipeline_dir = "stringtie_sq3", pipeline_id = "stringtie2",         is_ml = FALSE),
  list(row_label = "FLAIR",            pipeline_dir = "flair_sq3",    pipeline_id = "flair",              is_ml = FALSE),
  list(row_label = "Iso-Seq + SQ3 ML", pipeline_dir = "isoseq_sq3",    pipeline_id = "isoseq_sqanti3ml",   is_ml = FALSE)
)

pipeline_order <- vapply(pipeline_row_specs, function(x) x$pipeline_id, character(1))

# Collect grobs row-wise
overall_radar_grobs <- list()
for (row_spec in pipeline_row_specs) {
  pipeline_path <- file.path(base_data_dir, row_spec$pipeline_dir)
  for (sample_prefix in overall_sample_prefixes) {
    species <- if (grepl("^wtc11", sample_prefix, ignore.case = TRUE)) "human" else "mouse"
    tusco_file <- if (species == "human") tusco_human_file else tusco_mouse_file
    grob_cell <- tryCatch({
      process_and_plot_pipeline_pair(
        pipeline_prefix    = sample_prefix,
        species            = species,
        tusco_ref_filepath = tusco_file,
        data_dir_base      = pipeline_path,
        is_ml_processing   = row_spec$is_ml,
        omit_title         = TRUE,
        pipeline_display   = row_spec$row_label,
        pipeline_identifier = row_spec$pipeline_id
      )
    }, error = function(e) ggdraw() + draw_label("Error", fontface = "bold"))
    overall_radar_grobs <- append(overall_radar_grobs, list(grob_cell))
  }
}

# Bottom column labels (generic, repeated for both species)
bottom_column_labels <- c("CapTrap\nPacBio", "cDNA\nPacBio", "cDNA\nONT", "CapTrap\nONT",
                          "CapTrap\nPacBio", "cDNA\nPacBio", "cDNA\nONT", "CapTrap\nONT")
bottom_label_grobs <- lapply(bottom_column_labels, function(txt) {
  ggdraw() + draw_label(txt, fontface = "plain", size = 7, fontfamily = "Helvetica")
})

# Build rows with row labels (wider left label and larger text)
row_labels <- vapply(pipeline_row_specs, function(x) x$row_label, character(1))
# Wrap long row label to three lines: "Iso-Seq", "+", "SQ3 ML"
row_labels[row_labels == "Iso-Seq + SQ3 ML"] <- "Iso-Seq\n+\nSQ3 ML"
row_panels <- list()
for (i in seq_along(row_labels)) {
  row_label_grob <- ggdraw() +
    draw_label(row_labels[i], fontfamily = "Helvetica", fontface = "bold", size = 7, angle = 0, hjust = 1) +
    theme(plot.margin = unit(c(0, 1, 0, 2), "mm"))
  row_grobs <- overall_radar_grobs[((i-1)*8 + 1):(i*8)]
  row_panel <- plot_grid(plotlist = c(list(row_label_grob), row_grobs), ncol = 9, rel_widths = c(0.9, rep(1,8)))
  row_panels[[i]] <- row_panel
}

# Group headers at top: Human (first 4 cols) and Mouse (last 4 cols)
blank_corner <- ggdraw()
human_header <- make_group_header("Human", title_size = 9)
mouse_header <- make_group_header("Mouse", title_size = 9)
group_header_row <- plot_grid(blank_corner, human_header, mouse_header, ncol = 3, rel_widths = c(0.9, 4, 4))

# Bottom row with column labels (force center alignment under each radar cell)
bottom_labels_row <- plot_grid(
  plotlist = c(list(blank_corner), bottom_label_grobs),
  ncol = 9,
  rel_widths = c(0.9, rep(1, 8)),
  align = "h",
  axis = "tb"
)

# Combine header + rows + bottom labels, add legend
overall_grid <- plot_grid(plotlist = c(list(group_header_row), row_panels, list(bottom_labels_row)), ncol = 1, rel_heights = c(0.15, rep(1, length(row_panels)), 0.5))
overall_grid_padded <- overall_grid + theme(plot.margin = unit(c(10, 2, 2, 6), "mm"))
overall_grid_with_legend <- plot_grid(overall_grid_padded, radar_plot_legend, ncol = 1, rel_heights = c(1, 0.18))

# Save as fig-4b.pdf (width 180 mm; height derived from 4 rows + headers + legend)
cell_width_mm <- 180 / 8
overall_height_mm <- cell_width_mm * 5.6 + 12
  figure_path <- file.path(output_dir, "fig-4b.pdf")
# Convert dimensions from mm to inches
width_inches <- 180/25.4
height_inches <- overall_height_mm/25.4
ggsave(figure_path, overall_grid_with_legend, width = width_inches, height = height_inches, units = "in", device = "pdf", limitsize = FALSE)
message("Saved plot: ", figure_path)

# Write TSV of underlying metrics/metadata
tsv_path <- file.path(tsv_dir, "fig-4b.tsv")
if (length(metrics_records) > 0) {
  metrics_df <- tryCatch({
    bind_rows(metrics_records)
  }, error = function(e) NULL)
  if (!is.null(metrics_df)) {
    # Ensure standard column order
    std_cols <- c("figure_id", "pipeline", "sample", "species", "eval_type", "status")
    other_cols <- setdiff(colnames(metrics_df), std_cols)
    metrics_df <- metrics_df[, c(std_cols, other_cols)]
    write_tsv(metrics_df, tsv_path)
    message("Saved data: ", tsv_path)
  } else {
    # Fallback minimal TSV when no records (shouldn't happen)
    write_tsv(tibble(figure_id = "fig-4b", note = "no_metrics_records"), tsv_path)
    message("Saved data: ", tsv_path)
  }
} else {
  # No cells processed or no data everywhere
  write_tsv(tibble(figure_id = "fig-4b", status = "no_data"), tsv_path)
  message("Saved data: ", tsv_path)
}

area_ratio_tsv_path <- file.path(tsv_dir, "fig-4b_area_ratio.tsv")
if (length(area_ratio_records) > 0) {
  area_ratio_df <- tryCatch({
    bind_rows(area_ratio_records)
  }, error = function(e) NULL)
  if (!is.null(area_ratio_df)) {
    area_ratio_summary <- area_ratio_df %>%
      mutate(
        area_reference = as.numeric(area_reference),
        area_novel = as.numeric(area_novel),
        performance_drop_pct = as.numeric(performance_drop_pct)
      ) %>%
      mutate(species = factor(species, levels = c("human", "mouse"))) %>%
      group_by(pipeline_id, pipeline_label, species) %>%
      summarise(
        n_pairs = sum(!is.na(performance_drop_pct)),
        mean_performance_drop_pct = mean(performance_drop_pct, na.rm = TRUE),
        median_performance_drop_pct = median(performance_drop_pct, na.rm = TRUE),
        median_area_reference = median(area_reference, na.rm = TRUE),
        median_area_novel = median(area_novel, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        mean_performance_drop_pct = ifelse(is.nan(mean_performance_drop_pct), NA_real_, mean_performance_drop_pct),
        median_performance_drop_pct = ifelse(is.nan(median_performance_drop_pct), NA_real_, median_performance_drop_pct),
        median_area_reference = ifelse(is.nan(median_area_reference), NA_real_, median_area_reference),
        median_area_novel = ifelse(is.nan(median_area_novel), NA_real_, median_area_novel)
      ) %>%
      mutate(pipeline_id = factor(pipeline_id, levels = pipeline_order)) %>%
      arrange(pipeline_id, species) %>%
      mutate(
        tool_id = as.character(pipeline_id),
        tool = pipeline_label
      ) %>%
      select(tool_id, tool, species, n_pairs, mean_performance_drop_pct, median_performance_drop_pct,
             median_area_reference, median_area_novel)
    write_tsv(area_ratio_summary, area_ratio_tsv_path)
    message("Saved area ratio summary: ", area_ratio_tsv_path)
  } else {
    write_tsv(tibble(figure_id = "fig-4b", status = "invalid_area_records"), area_ratio_tsv_path)
    message("Saved area ratio summary: ", area_ratio_tsv_path)
  }
} else {
  write_tsv(tibble(figure_id = "fig-4b", status = "no_area_records"), area_ratio_tsv_path)
  message("Saved area ratio summary: ", area_ratio_tsv_path)
}

message("Done.")
