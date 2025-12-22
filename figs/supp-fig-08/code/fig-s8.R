#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(cowplot)
  library(patchwork)
})

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

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- if (length(script_path) == 1 && nzchar(script_path)) dirname(normalizePath(script_path, mustWork = FALSE)) else getwd()
repo_root <- find_repo_root_for_utils(script_dir)
source(file.path(repo_root, "scripts", "figure_utils.R"))

ctx <- figure_context(defaults = list(out_dir = "..", width = NATURE_DIMS$double_column, height = 7.5))
reviewer3_root <- dirname(ctx$figure_dir)

pdf_device <- if ("cairo_pdf" %in% getNamespaceExports("grDevices")) grDevices::cairo_pdf else "pdf"

read_tsv_quiet <- function(path) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)

SERIES_LEVELS <- c("TUSCO (Human)", "TUSCO (Mouse)", "GENCODE (Human)", "GENCODE (Mouse)")
SERIES_COLORS <- c(
  "TUSCO (Human)" = TUSCO_COLORS$human_tusco,
  "TUSCO (Mouse)" = TUSCO_COLORS$mouse_tusco,
  "GENCODE (Human)" = TUSCO_COLORS$human_gencode,
  "GENCODE (Mouse)" = TUSCO_COLORS$mouse_gencode
)

make_tss_ecdf <- function(tss_dir, threshold = 300L, show_legend = TRUE) {
  hsa_within <- read_tsv_quiet(file.path(tss_dir, "tss_within_gene_distances_corrected_human.tsv"))
  hsa_between <- read_tsv_quiet(file.path(tss_dir, "tss_between_gene_distances_corrected_human.tsv"))
  mmu_within <- read_tsv_quiet(file.path(tss_dir, "tss_within_gene_distances_corrected_mouse.tsv"))
  mmu_between <- read_tsv_quiet(file.path(tss_dir, "tss_between_gene_distances_corrected_mouse.tsv"))

  ecdf_data <- bind_rows(
    transmute(hsa_within, distance = distance_to_reftss, Series = "TUSCO (Human)"),
    transmute(mmu_within, distance = distance_to_reftss, Series = "TUSCO (Mouse)"),
    transmute(hsa_between, distance = distance_to_other_gene, Series = "GENCODE (Human)"),
    transmute(mmu_between, distance = distance_to_other_gene, Series = "GENCODE (Mouse)")
  ) %>%
    mutate(
      distance_plot = pmax(as.numeric(distance), 1),
      Series = factor(Series, levels = SERIES_LEVELS)
    )

  p <- ggplot(ecdf_data, aes(x = distance_plot, color = Series)) +
    stat_ecdf(geom = "step", linewidth = 0.4, pad = FALSE) +
    geom_vline(xintercept = threshold, linetype = "dotted", color = "gray30", linewidth = 0.25) +
    annotate("text",
      x = threshold,
      y = 0.92,
      label = paste0(threshold, "bp"),
      hjust = 0.5,
      vjust = 0,
      size = 2.4,
      color = "gray30"
    ) +
    scale_x_log10(
      breaks = c(1, 10, 100, 300, 1000, 10000, 100000, 300000),
      labels = scales::label_number(scale_cut = scales::cut_short_scale(), accuracy = 1),
      expand = expansion(mult = c(0.03, 0.08))
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.08))
    ) +
    scale_color_manual(values = SERIES_COLORS, breaks = SERIES_LEVELS) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0.8))) +
    coord_cartesian(xlim = c(1, 4e5), ylim = c(0, 1.02), clip = "off") +
    labs(x = "Window size (bp)", y = "Cumulative fraction") +
    theme_tusco(base_size = 7) +
    theme(
      legend.title = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "horizontal",
      legend.key.width = unit(0.9, "lines"),
      legend.spacing.x = unit(0.35, "lines"),
      plot.margin = margin(2, 6, 2, 2)
    )
  p
}

make_tss_peaks <- function(tss_dir, threshold = 300L, show_legend = FALSE) {
  hsa_peaks <- read_tsv_quiet(file.path(tss_dir, "reftss_peak_counts_by_window_human.tsv"))
  mmu_peaks <- read_tsv_quiet(file.path(tss_dir, "reftss_peak_counts_by_window_mouse.tsv"))

  peak_data <- bind_rows(
    mutate(hsa_peaks, Species = "Human"),
    mutate(mmu_peaks, Species = "Mouse")
  ) %>%
    mutate(q75_plot = if_else(as.numeric(q75) > 0, as.numeric(q75), NA_real_))

  p <- ggplot(peak_data, aes(x = window_size, y = q75_plot, color = Species)) +
    geom_line(linewidth = 0.4, na.rm = TRUE) +
    geom_point(size = 1, na.rm = TRUE) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray30", linewidth = 0.25) +
    geom_vline(xintercept = threshold, linetype = "dotted", color = "gray30", linewidth = 0.25) +
    annotate("text",
      x = threshold,
      y = max(peak_data$q75_plot[peak_data$window_size <= 1500], na.rm = TRUE) * 0.98,
      label = paste0(threshold, "bp"),
      hjust = 0.5,
      vjust = -0.5,
      size = 2.4,
      color = "gray30"
    ) +
    scale_x_log10(
      breaks = c(1, 10, 100, 300, 1000, 10000, 100000),
      labels = scales::comma,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_log10(
      breaks = c(0.1, 0.5, 1, 2, 5, 10, 20),
      labels = scales::comma,
      expand = c(0.02, 0.02)
    ) +
    scale_color_manual(
      values = c("Human" = TUSCO_COLORS$human_tusco, "Mouse" = TUSCO_COLORS$mouse_tusco),
      guide = "none"
    ) +
    labs(x = "Window size (bp)", y = "Peaks per gene (Q75)") +
    theme_tusco(base_size = 7) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      plot.margin = margin(2, 2, 2, 2)
    )
  p
}

make_sirv_config <- function(sirv_dir) {
  df <- read_tsv_quiet(file.path(sirv_dir, "sirv_config_comparison.tsv"))

  df <- df %>%
    mutate(
      species_label = if_else(species == "hsa", "Human", "Mouse"),
      config_label = factor(config_label, levels = c("Strict", "Current", "Lenient", "More Lenient"))
    )

  # Use distinct gradients for human vs mouse (no legend; facets keep it readable)
  config_colors_human <- c(
    "Strict" = "#f7fcf5",
    "Current" = "#c7e9c0",
    "Lenient" = "#74c476",
    "More Lenient" = "#238b45"
  )
  config_colors_mouse <- c(
    "Strict" = "#3fb08f",
    "Current" = TUSCO_COLORS$mouse_tusco,
    "Lenient" = "#0f7a5b",
    "More Lenient" = "#055844"
  )
  df <- df %>%
    mutate(fill_key = paste0(species_label, "|", as.character(config_label)))
  config_colors <- c(
    setNames(config_colors_human, paste0("Human|", names(config_colors_human))),
    setNames(config_colors_mouse, paste0("Mouse|", names(config_colors_mouse)))
  )

  baselines <- df %>%
    group_by(species_label) %>%
    summarise(median_sirv = first(median_sirv_fn_pct), .groups = "drop")

  p <- ggplot(df, aes(x = config_label, y = median_tusco_fn_pct, fill = fill_key)) +
    geom_col(color = "black", width = 0.7, linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = q1_tusco_fn_pct, ymax = q3_tusco_fn_pct),
      width = 0.2,
      linewidth = 0.3
    ) +
    geom_hline(
      data = baselines,
      aes(yintercept = median_sirv),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = TUSCO_COLORS$sirvs,
      linewidth = 0.5
    ) +
    facet_wrap(~species_label, scales = "free_y") +
    scale_fill_manual(values = config_colors, guide = "none") +
    labs(x = "TSS threshold", y = "Median FN rate (%)") +
    theme_tusco(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95", linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2)
    )
  p
}

make_alphagenome_density <- function(alphagenome_dir) {
  plot_df <- read_tsv_quiet(file.path(alphagenome_dir, "tables", "alphagenome_expression_data.tsv"))

  plot_df <- plot_df %>%
    mutate(
      species = factor(species, levels = c("Human", "Mouse")),
      is_tusco = as.logical(is_tusco)
    )

  tusco_medians <- plot_df %>%
    filter(is_tusco) %>%
    group_by(species) %>%
    summarise(median_log10 = median(log10_expression, na.rm = TRUE), .groups = "drop")

  species_colors <- c("Human" = TUSCO_COLORS$human_tusco, "Mouse" = TUSCO_COLORS$mouse_tusco)

  p <- ggplot() +
    geom_density(
      data = filter(plot_df, !is_tusco),
      aes(x = log10_expression),
      fill = "#D1D3D4",
      color = "#999999",
      alpha = 0.7,
      linewidth = 0.3
    ) +
    geom_density(
      data = filter(plot_df, is_tusco),
      aes(x = log10_expression, fill = species, color = species),
      alpha = 0.7,
      linewidth = 0.4
    ) +
    geom_vline(
      data = tusco_medians,
      aes(xintercept = median_log10, color = species),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    facet_wrap(~species, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = species_colors, guide = "none") +
    scale_color_manual(values = species_colors, guide = "none") +
    scale_x_continuous(breaks = seq(-2, 2, by = 1), labels = scales::math_format(10^.x)) +
    labs(x = "Median RPKM across tissues", y = "Density") +
    theme_tusco(base_size = 7) +
    theme(
      strip.background = element_rect(fill = "grey95", linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2)
    )
  p
}

make_splice_ecdf <- function(splice_dir) {
  hsa <- read_tsv_quiet(file.path(splice_dir, "hsa_gene_junction_read_counts.tsv.gz"))
  mmu <- read_tsv_quiet(file.path(splice_dir, "mmu_gene_junction_read_counts.tsv.gz"))

  prep <- function(df, species) {
    df <- df %>%
      mutate(
        species = species,
        is_tusco = as.integer(is_tusco),
        annotated_intron_count = as.integer(annotated_intron_count),
        known_observed_junctions = as.integer(known_observed_junctions),
        novel_junctions = as.integer(novel_junctions),
        novel_max_over_known_mean = suppressWarnings(as.numeric(novel_max_over_known_mean))
      ) %>%
      filter(
        annotated_intron_count >= 1,
        (known_observed_junctions + novel_junctions) >= 1,
        is.finite(novel_max_over_known_mean),
        novel_max_over_known_mean >= 0
      ) %>%
      mutate(
        novel_max_over_known_mean_plot = pmax(novel_max_over_known_mean, 1e-4),
        Series = case_when(
          species == "Human" & is_tusco == 1 ~ "TUSCO (Human)",
          species == "Mouse" & is_tusco == 1 ~ "TUSCO (Mouse)",
          species == "Human" & is_tusco == 0 ~ "GENCODE (Human)",
          species == "Mouse" & is_tusco == 0 ~ "GENCODE (Mouse)",
          TRUE ~ NA_character_
        )
      )
    df
  }

  plot_df <- bind_rows(prep(hsa, "Human"), prep(mmu, "Mouse")) %>%
    mutate(species = factor(species, levels = c("Human", "Mouse")))

  alpha_df <- tibble::tibble(species = c("Human", "Mouse"), alpha = c(0.01, 0.05))

  plot_df <- plot_df %>%
    mutate(Series = factor(Series, levels = SERIES_LEVELS))

  p <- ggplot(plot_df, aes(x = novel_max_over_known_mean_plot, color = Series)) +
    stat_ecdf(geom = "step", linewidth = 0.4) +
    geom_vline(
      data = alpha_df,
      aes(xintercept = alpha),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.25,
      color = "black"
    ) +
    facet_wrap(~species, nrow = 1, scales = "fixed") +
    scale_x_log10() +
    coord_cartesian(xlim = c(1e-4, 1e1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = SERIES_COLORS, breaks = SERIES_LEVELS) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0.8))) +
    labs(x = "max(novel junction reads) / mean(annotated junction reads)", y = "ECDF across genes") +
    theme_tusco(base_size = 7) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95", linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2)
    )
  p
}

# Data directories - use absolute paths from repo root
reviewer3_data <- file.path(repo_root, "reviewer_response", "review_plots", "reviewer3")
tss_dir <- file.path(reviewer3_data, "TSS")
sirv_dir <- file.path(reviewer3_data, "expression", "lrgasp_evaluation")
splice_dir <- file.path(reviewer3_data, "splice")
alphagenome_dir <- file.path(reviewer3_data, "alphagenome")

p_tss_ecdf <- make_tss_ecdf(tss_dir, show_legend = TRUE)
p_tss_peaks <- make_tss_peaks(tss_dir, show_legend = FALSE)
p_sirv <- make_sirv_config(sirv_dir)
p_alpha <- make_alphagenome_density(alphagenome_dir)
p_splice <- make_splice_ecdf(splice_dir) + theme(legend.position = "none")

row1 <- (p_tss_ecdf + p_tss_peaks) + plot_layout(widths = c(1.15, 0.85))
row2 <- (p_sirv + p_alpha) + plot_layout(widths = c(1, 1))
combined <- row1 / row2 / p_splice +
  plot_layout(guides = "collect", heights = c(1, 1, 1.05)) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.95, "lines"),
    legend.spacing.x = unit(0.35, "lines"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
    legend.box.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
    legend.margin = margin(2, 4, 2, 4)
  )

combined <- combined + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

ctx$save_plot(combined, "fig-s8.pdf", width = ctx$params$width, height = ctx$params$height, units = "in", device = pdf_device)
ctx$save_plot(combined, "fig-s8.png", width = ctx$params$width, height = ctx$params$height, units = "in", dpi = 300)
