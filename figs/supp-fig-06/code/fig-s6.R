#!/usr/bin/env Rscript

###############################################################
# Figure S6 - SIRV metrics across replicate combinations
# Usage: Rscript fig-s6.R [output_dir] [width] [height]
###############################################################

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

ctx <- figure_context(defaults = list(
  out_dir = "..",
  width = 7.09,
  height = 3.55
))
params <- ctx$params

local_only <- tolower(Sys.getenv("LOCAL_ONLY", "0")) %in% c("1", "true", "t", "yes", "y")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(cowplot)
  library(scales)
  library(grid)
})

if (!local_only) {
  suppressPackageStartupMessages(library(rtracklayer))
}

plot_dir <- ctx$plot_dir
tsv_dir <- ctx$table_dir
out_pdf <- file.path(plot_dir, "fig-s6.pdf")

message("Output directories:")
message("  Plots: ", plot_dir)
message("  Tables: ", tsv_dir)

fig_dir <- ctx$figure_dir
repo_root <- ctx$repo_root

resolve_existing <- function(paths, type = c("file", "directory"), required = FALSE) {
  type <- match.arg(type)
  for (p in paths) {
    resolved <- ctx$resolve_input(candidates = p, type = if (type == "file") "file" else "directory", required = FALSE)
    if (type == "file") {
      if (!is.null(resolved) && file.exists(resolved)) return(resolved)
    } else {
      if (!is.null(resolved) && dir.exists(resolved)) return(resolved)
    }
  }
  if (required) {
    stop("Unable to locate required path. Checked: ", paste(paths, collapse = ", "))
  }
  NA_character_
}

intersection_all_comb_root <- ctx$resolve_input("nih", "intersection_all_comb", type = "directory", required = FALSE)
isoseq_ar_root <- ctx$resolve_input("nih", "single_sample", type = "directory", required = FALSE)
if (!isTRUE(dir.exists(intersection_all_comb_root))) intersection_all_comb_root <- NA_character_
if (!isTRUE(dir.exists(isoseq_ar_root))) isoseq_ar_root <- NA_character_

detect_sirv_gtf <- function() {
  resolve_existing(c(
    file.path("nih", "SIRVs.gtf"),
    file.path("spike-ins", "lrgasp_sirvs.gtf"),
    file.path("spike-ins", "lrgasp_sirvs4.gtf"),
    file.path("reference", "SIRVs.gtf")
  ), type = "file")
}

sirv_file <- if (!local_only) detect_sirv_gtf() else NA_character_
if (!local_only && (is.null(sirv_file) || is.na(sirv_file) || !file.exists(sirv_file))) {
  stop("SIRV GTF not found. Checked spike-ins and reference locations under: ",
       paste(unique(ctx$data_roots), collapse = ", "))
}

read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  tryCatch(read_tsv(file_path, show_col_types = FALSE, ...),
           error = function(e) stop("Unable to read ", file_path, ": ", e$message))
}

find_classification_file <- function(dir_path) {
  files <- list.files(dir_path, full.names = TRUE)
  f <- files[grepl("(?<!_union)_classification\\.txt$", files, perl = TRUE)]
  if (length(f) > 0) return(f[1])
  f2 <- files[grepl("_union_classification\\.txt$", files)]
  if (length(f2) > 0) return(f2[1])
  f3 <- files[grepl("classification\\.txt$", files)]
  if (length(f3) > 0) return(f3[1])
  stop("No classification file found in ", dir_path)
}

compute_sirv_metrics_from_classification <- function(classification_data) {
  sirv_gtf_df <- as.data.frame(rtracklayer::import(sirv_file))
  sirv_exons <- sirv_gtf_df[sirv_gtf_df$type == "exon", ]
  rSIRV <- length(unique(sirv_exons$transcript_id))

  classification_data_cleaned_sirv <- classification_data %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(classification_data %>%
                dplyr::filter(structural_category == "fusion") %>%
                tidyr::separate_rows(associated_transcript, sep = "_")) %>%
    dplyr::mutate(
      associated_gene = stringr::str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = stringr::str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    dplyr::distinct(isoform, associated_transcript, .keep_all = TRUE)

  sirv_chromosomes <- unique(sirv_gtf_df$seqnames)
  classification_data_cleaned_sirv <- classification_data_cleaned_sirv %>%
    dplyr::filter(chrom %in% sirv_chromosomes)
  SIRV_transcripts <- classification_data_cleaned_sirv %>% dplyr::filter(grepl("SIRV", chrom))
  SIRV_RM <- SIRV_transcripts %>%
    dplyr::filter(
      subcategory == "reference_match" |
      (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
         !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
         abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50) |
      (!is.na(ref_length) & ref_length > 3000 &
       structural_category == "full-splice_match" &
       !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
       abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100)
    )
  TP_sirv <- SIRV_RM
  PTP_sirv <- SIRV_transcripts %>%
    dplyr::filter((structural_category %in% c("full-splice_match", "incomplete-splice_match") |
                    subcategory == "mono-exon_by_intron_retention") &
             !associated_transcript %in% TP_sirv$associated_transcript)
  fsm_ism_count_sirv <- SIRV_transcripts %>%
    dplyr::filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    nrow()

  non_redundant_sensitivity_sirv <- length(unique(TP_sirv$associated_transcript)) / rSIRV
  non_redundant_precision_sirv  <- if (nrow(SIRV_transcripts) > 0) nrow(TP_sirv) / nrow(SIRV_transcripts) else NA
  redundancy_sirv              <- if (length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) > 0)
                                     fsm_ism_count_sirv / length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) else NA
  false_discovery_rate_sirv    <- if (nrow(SIRV_transcripts) > 0) (nrow(SIRV_transcripts) - nrow(SIRV_RM)) / nrow(SIRV_transcripts) else NA
  positive_detection_rate_sirv <- if (rSIRV > 0) length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) / rSIRV else NA
  f1_sirv <- if (!is.na(non_redundant_sensitivity_sirv) && !is.na(non_redundant_precision_sirv) &&
                  (non_redundant_sensitivity_sirv + non_redundant_precision_sirv) > 0) {
    2 * non_redundant_sensitivity_sirv * non_redundant_precision_sirv /
      (non_redundant_sensitivity_sirv + non_redundant_precision_sirv)
  } else { NA_real_ }

  c(
    Sensitivity = non_redundant_sensitivity_sirv * 100,
    `Non-redundant Precision` = non_redundant_precision_sirv * 100,
    `Inv. Redundancy` = if (!is.na(redundancy_sirv) && redundancy_sirv != 0) (1/redundancy_sirv) * 100 else 0,
    `1 - FDR` = if (!is.na(false_discovery_rate_sirv)) (100 - (false_discovery_rate_sirv * 100)) else 0,
    `PDR` = if (!is.na(positive_detection_rate_sirv)) positive_detection_rate_sirv * 100 else 0,
    `F1` = if (!is.na(f1_sirv)) f1_sirv * 100 else NA_real_
  )
}

collect_points_for_sirv <- function() {
  single_points <- list()
  if (!is.na(isoseq_ar_root) && dir.exists(isoseq_ar_root)) {
    single_sample_dirs <- list.dirs(isoseq_ar_root, full.names = TRUE, recursive = FALSE)
    single_sample_dirs <- single_sample_dirs[grepl("/(B|K)[0-9]+\\.isoforms$", single_sample_dirs)]
    for (sd in single_sample_dirs) {
      sample_name <- basename(sd)
      class_file  <- file.path(sd, paste0(sample_name, "_classification.txt"))
      if (!file.exists(class_file)) next
      tissue <- if (startsWith(sample_name, "B")) "Brain" else if (startsWith(sample_name, "K")) "Kidney" else "Unknown"
      classification <- read_tsv_safe(class_file)
      mv <- compute_sirv_metrics_from_classification(classification)
      df <- tibble(
        Tissue = tissue,
        Samples = 1L,
        Combo = sample_name,
        Metric = names(mv),
        Value = as.numeric(mv)
      )
      single_points[[length(single_points) + 1]] <- df
    }
  }
  single_points_df <- if (length(single_points) > 0) {
    dplyr::bind_rows(single_points)
  } else {
    tibble::tibble(
      Tissue = character(),
      Samples = integer(),
      Combo = character(),
      Metric = character(),
      Value = double()
    )
  }

  combo_points <- list()
  if (!is.na(intersection_all_comb_root) && dir.exists(intersection_all_comb_root)) {
    combo_dirs <- list.dirs(intersection_all_comb_root, full.names = TRUE, recursive = FALSE)
    combo_dirs <- combo_dirs[grepl("/(B|K)[0-9]+(_(B|K)[0-9]+)*$", combo_dirs)]
    is_pure_tissue <- function(name) {
      parts <- strsplit(name, "_")[[1]]
      starts_with_b <- all(startsWith(parts, "B"))
      starts_with_k <- all(startsWith(parts, "K"))
      starts_with_b || starts_with_k
    }
    combo_dirs <- combo_dirs[vapply(basename(combo_dirs), is_pure_tissue, logical(1))]

    for (cd in combo_dirs) {
      combo_name <- basename(cd)
      parts <- strsplit(combo_name, "_")[[1]]
      tissue <- if (all(startsWith(parts, "B"))) "Brain" else if (all(startsWith(parts, "K"))) "Kidney" else "Mixed"
      if (tissue == "Mixed") next
      samples_n <- length(parts)
      cf <- find_classification_file(cd)
      classification <- read_tsv_safe(cf)
      mv <- compute_sirv_metrics_from_classification(classification)
      df <- tibble(
        Tissue = tissue,
        Samples = samples_n,
        Combo = combo_name,
        Metric = names(mv),
        Value = as.numeric(mv)
      )
      combo_points[[length(combo_points) + 1]] <- df
    }
  }

  combo_points_df <- if (length(combo_points) > 0) {
    dplyr::bind_rows(combo_points)
  } else {
    tibble::tibble(
      Tissue = character(),
      Samples = integer(),
      Combo = character(),
      Metric = character(),
      Value = double()
    )
  }

  bind_rows(single_points_df, combo_points_df)
}

summarize_points_df <- function(points_df) {
  points_df %>%
    dplyr::group_by(Tissue, Samples, Metric) %>%
    dplyr::summarise(
      N = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SE = ifelse(N > 0, SD / sqrt(pmax(N, 1)), NA_real_),
      df_t = ifelse(!is.na(N) & N > 1, pmax(N - 1, 1), 1),
      use_ci = !is.na(N) & N > 1,
      CI_t_raw = qt(0.975, df = df_t),
      CI_t = ifelse(use_ci, CI_t_raw, NA_real_),
      CI_lower = ifelse(use_ci, Mean - CI_t_raw * SE, NA_real_),
      CI_upper = ifelse(use_ci, Mean + CI_t_raw * SE, NA_real_),
      SamplesLabel = as.character(Samples)
    ) %>%
    dplyr::select(-df_t, -use_ci, -CI_t_raw)
}

metric_order_sirv <- c("Sensitivity", "Non-redundant Precision", "Inv. Redundancy", "1 - FDR", "PDR", "F1")
points_s6_path <- file.path(plot_dir, "fig-s6_points.tsv")
bars_s6_path   <- file.path(plot_dir, "fig-s6_bars.tsv")
points_s6_tsv  <- file.path(tsv_dir,  "fig-s6_points.tsv")
bars_s6_tsv    <- file.path(tsv_dir,  "fig-s6_bars.tsv")

if (!local_only) {
  all_points_sirv_df <- collect_points_for_sirv()
  if (!all(c("Tissue","Samples","Combo","Metric","Value") %in% names(all_points_sirv_df))) {
    warning("No points found or missing columns; proceeding with empty dataset for fig-s6.")
    all_points_sirv_df <- tibble::tibble(Tissue = character(), Samples = integer(), Combo = character(), Metric = character(), Value = double())
  }
  all_points_sirv_df <- all_points_sirv_df %>%
    filter(Tissue %in% c("Brain", "Kidney")) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Brain", "Kidney")),
      Metric = factor(Metric, levels = metric_order_sirv),
      SamplesLabel = as.character(Samples)
    )
  summary_sirv_df <- summarize_points_df(all_points_sirv_df)

  readr::write_tsv(all_points_sirv_df, points_s6_path)
  readr::write_tsv(summary_sirv_df, bars_s6_path)
  readr::write_tsv(all_points_sirv_df, points_s6_tsv)
  readr::write_tsv(summary_sirv_df, bars_s6_tsv)
} else {
  if (!file.exists(points_s6_path) || !file.exists(bars_s6_path)) {
    stop("LOCAL_ONLY is set, but required SIRV TSVs are missing in ", plot_dir, ". Expected: ",
         basename(points_s6_path), ", ", basename(bars_s6_path))
  }
  all_points_sirv_df <- readr::read_tsv(points_s6_path, show_col_types = FALSE) %>%
    mutate(
      Tissue = factor(Tissue, levels = c("Brain", "Kidney")),
      Metric = factor(Metric, levels = metric_order_sirv),
      SamplesLabel = as.character(Samples)
    )
  summary_sirv_df <- readr::read_tsv(bars_s6_path, show_col_types = FALSE) %>%
    mutate(SamplesLabel = as.character(Samples))
}

create_metric_barplot_all <- function(metric_name, points_df, summary_df, color_value = "#2ca02c") {
  metric_points <- points_df %>% dplyr::filter(Metric == metric_name)
  metric_summary <- summary_df %>% dplyr::filter(Metric == metric_name)

  x_levels <- as.character(sort(unique(points_df$Samples)))
  x_levels <- intersect(as.character(1:5), x_levels)
  metric_points$SamplesLabel <- factor(metric_points$SamplesLabel, levels = as.character(1:5))
  metric_summary$SamplesLabel <- factor(metric_summary$SamplesLabel, levels = as.character(1:5))

  improvement_data <- metric_summary %>%
    dplyr::filter(Samples %in% c(1, 3), Tissue %in% c("Brain", "Kidney")) %>%
    dplyr::arrange(Tissue, Samples) %>%
    dplyr::group_by(Tissue) %>%
    dplyr::filter(dplyr::n() == 2) %>%
    dplyr::summarise(
      improvement = Mean[Samples == 3] - Mean[Samples == 1],
      sample1_value = Mean[Samples == 1],
      sample3_value = Mean[Samples == 3],
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      improvement_text = ifelse(improvement >= 0, paste0("+", round(improvement, 1), "pp"), paste0(round(improvement, 1), "pp")),
      annotation_y = pmax(sample1_value, sample3_value) + 18,
      bracket_y = pmax(sample1_value, sample3_value) + 8,
      x_start = factor("1", levels = as.character(1:5)),
      x_end   = factor("3", levels = as.character(1:5)),
      x_mid   = factor("2", levels = as.character(1:5))
    )

  facet_levels <- c("Brain", "Kidney")
  base_facet_df <- data.frame(Tissue = factor(facet_levels, levels = facet_levels))

  ggplot(base_facet_df) +
    geom_bar(
      data = metric_summary,
      aes(x = SamplesLabel, y = Mean, fill = Tissue),
      stat = "identity", position = position_dodge2(width = 0.75, preserve = "single"),
      width = 0.7, color = "white", linewidth = 0.3
    ) +
    geom_errorbar(
      data = metric_summary,
      aes(x = SamplesLabel, ymin = pmax(0, CI_lower), ymax = pmin(130, CI_upper), group = Tissue),
      position = position_dodge2(width = 0.75, preserve = "single"),
      width = 0.5, linewidth = 0.5, color = "black"
    ) +
    geom_point(
      data = metric_points,
      aes(x = SamplesLabel, y = Value, group = Tissue),
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      size = 0.35, alpha = 0.45, color = "#1f1f1f"
    ) +
    geom_text(
      data = improvement_data,
      aes(x = x_mid, y = annotation_y, label = improvement_text),
      inherit.aes = FALSE,
      hjust = 0.5, vjust = 0.5, size = 7/.pt, fontface = "bold", color = "#d62728"
    ) +
    geom_segment(
      data = improvement_data,
      aes(x = x_start, xend = x_end, y = bracket_y + 2, yend = bracket_y + 2, group = Tissue),
      inherit.aes = FALSE, color = "#d62728", linewidth = 0.4
    ) +
    geom_segment(
      data = improvement_data,
      aes(x = x_start, xend = x_start, y = bracket_y, yend = bracket_y + 2, group = Tissue),
      inherit.aes = FALSE, color = "#d62728", linewidth = 0.4
    ) +
    geom_segment(
      data = improvement_data,
      aes(x = x_end, xend = x_end, y = bracket_y, yend = bracket_y + 2, group = Tissue),
      inherit.aes = FALSE, color = "#d62728", linewidth = 0.4
    ) +
    facet_grid(. ~ Tissue, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(Brain = alpha(color_value, 0.90), Kidney = alpha(color_value, 0.70))) +
    scale_y_continuous(limits = c(0, 130), labels = scales::percent_format(scale = 1), breaks = seq(0, 100, 20)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
    coord_cartesian(clip = "off") +
    labs(title = metric_name, x = NULL, y = "Performance (%)") +
    theme_classic(base_family = "Helvetica", base_size = 7) +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5, face = "plain", margin = margin(b = 6)),
      axis.title.x = element_blank(), axis.title.y = element_text(size = 7, margin = margin(r = 4)),
      axis.text = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.8),
      axis.line.y = element_blank(),
      strip.text = element_text(size = 7, face = "plain", margin = margin(t = 4, b = 4)),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.spacing = unit(0.5, "cm"), plot.margin = margin(2, 2, 2, 2), legend.position = "none"
    )
}

metric_levels_sirv <- metric_order_sirv
plots_s6 <- lapply(metric_levels_sirv, function(m) create_metric_barplot_all(m, all_points_sirv_df, summary_sirv_df, "#cab2d6"))
fig_s6_panel <- plot_grid(plotlist = plots_s6, ncol = 3, nrow = 2, align = "hv")

sirv_underlying_out <- file.path(tsv_dir, "fig-s6.tsv")
sirv_all_points_out <- all_points_sirv_df %>%
  arrange(Tissue, Samples, Combo, Metric) %>%
  mutate(figure_id = "fig-s6", panel_id = "s6_points")
sirv_bar_means_out <- summary_sirv_df %>%
  arrange(Tissue, Samples, Metric) %>%
  select(Tissue, Samples, Metric, N, Mean, SD, SE, CI_lower, CI_upper) %>%
  mutate(figure_id = "fig-s6", panel_id = "s6_bars")
sirv_export_df <- bind_rows(
  sirv_all_points_out %>% mutate(dataset = "s6_points"),
  sirv_bar_means_out %>% mutate(dataset = "s6_bars")
)
readr::write_tsv(sirv_export_df, sirv_underlying_out)
message("Saved data: ", sirv_underlying_out)

pdf(out_pdf, width = params$width, height = params$height)
print(fig_s6_panel)
dev.off()
message("Saved plot: ", out_pdf)

extraction_file <- file.path(repo_root, "manuscript", "extraction_values.tsv")
extraction_data_list <- list()

if (exists("summary_sirv_df") && nrow(summary_sirv_df) > 0) {
  sirv_brain <- summary_sirv_df %>% filter(Tissue == "Brain")
  sirv_kidney <- summary_sirv_df %>% filter(Tissue == "Kidney")

  brain_1rep_sens <- sirv_brain %>% filter(Samples == 1, Metric == "Sensitivity")
  brain_1rep_prec <- sirv_brain %>% filter(Samples == 1, Metric == "Non-redundant Precision")

  if (nrow(brain_1rep_sens) > 0) {
    extraction_data_list[["sirv_brain_1rep_sens"]] <- data.frame(
      figure_id = "fig-s6",
      metric_name = "brain_sirv_1rep_sensitivity",
      value = formatC(brain_1rep_sens$Mean, format="f", digits=1),
      ci_lower = if (!is.na(brain_1rep_sens$CI_lower)) formatC(brain_1rep_sens$CI_lower, format="f", digits=1) else NA,
      ci_upper = if (!is.na(brain_1rep_sens$CI_upper)) formatC(brain_1rep_sens$CI_upper, format="f", digits=1) else NA,
      notes = "SIRV brain 1-replicate sensitivity (95% CI)",
      stringsAsFactors = FALSE
    )
  }

  if (nrow(brain_1rep_prec) > 0) {
    extraction_data_list[["sirv_brain_1rep_prec"]] <- data.frame(
      figure_id = "fig-s6",
      metric_name = "brain_sirv_1rep_precision",
      value = formatC(brain_1rep_prec$Mean, format="f", digits=1),
      ci_lower = if (!is.na(brain_1rep_prec$CI_lower)) formatC(brain_1rep_prec$CI_lower, format="f", digits=1) else NA,
      ci_upper = if (!is.na(brain_1rep_prec$CI_upper)) formatC(brain_1rep_prec$CI_upper, format="f", digits=1) else NA,
      notes = "SIRV brain 1-replicate precision (95% CI)",
      stringsAsFactors = FALSE
    )
  }

  kidney_1rep_sens <- sirv_kidney %>% filter(Samples == 1, Metric == "Sensitivity")
  kidney_1rep_prec <- sirv_kidney %>% filter(Samples == 1, Metric == "Non-redundant Precision")

  if (nrow(kidney_1rep_sens) > 0) {
    extraction_data_list[["sirv_kidney_1rep_sens"]] <- data.frame(
      figure_id = "fig-s6",
      metric_name = "kidney_sirv_1rep_sensitivity",
      value = formatC(kidney_1rep_sens$Mean, format="f", digits=1),
      ci_lower = if (!is.na(kidney_1rep_sens$CI_lower)) formatC(kidney_1rep_sens$CI_lower, format="f", digits=1) else NA,
      ci_upper = if (!is.na(kidney_1rep_sens$CI_upper)) formatC(kidney_1rep_sens$CI_upper, format="f", digits=1) else NA,
      notes = "SIRV kidney 1-replicate sensitivity (95% CI)",
      stringsAsFactors = FALSE
    )
  }

  if (nrow(kidney_1rep_prec) > 0) {
    extraction_data_list[["sirv_kidney_1rep_prec"]] <- data.frame(
      figure_id = "fig-s6",
      metric_name = "kidney_sirv_1rep_precision",
      value = formatC(kidney_1rep_prec$Mean, format="f", digits=1),
      ci_lower = if (!is.na(kidney_1rep_prec$CI_lower)) formatC(kidney_1rep_prec$CI_lower, format="f", digits=1) else NA,
      ci_upper = if (!is.na(kidney_1rep_prec$CI_upper)) formatC(kidney_1rep_prec$CI_upper, format="f", digits=1) else NA,
      notes = "SIRV kidney 1-replicate precision (95% CI)",
      stringsAsFactors = FALSE
    )
  }
}

if (length(extraction_data_list) > 0) {
  extraction_data <- bind_rows(extraction_data_list)

  if (!file.exists(extraction_file)) {
    write.table(extraction_data, extraction_file, sep = "\t",
                row.names = FALSE, quote = FALSE)
    cat("Created extraction file:", extraction_file, "\n")
  } else {
    write.table(extraction_data, extraction_file, sep = "\t",
                row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)
  }

  cat("Extracted", nrow(extraction_data), "SIRV metrics to:", extraction_file, "\n")
} else {
  cat("No data available for extraction\n")
}
