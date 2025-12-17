#!/usr/bin/env Rscript

###############################################################
# Figure 3A Human (Alternative) - Dumbbell metrics plot
# Reads existing fig-3a-human.tsv and renders a facet dumbbell plot
# Usage: Rscript fig-3a-dumbbell-human.R [output_dir] [width] [height]
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Prefer repo helper functions when available (keeps behavior aligned with other figure scripts)
if (file.exists("../../../scripts/figure_utils.R")) {
  source("../../../scripts/figure_utils.R")
} else {
  parse_figure_args <- function(defaults = list(out_dir = "..", width = 11, height = 4.6)) {
    args <- commandArgs(trailingOnly = TRUE)
    result <- defaults
    if (length(args) > 0) result$out_dir <- args[1]
    if (length(args) > 1) result$width <- as.numeric(args[2])
    if (length(args) > 2) result$height <- as.numeric(args[3])
    result
  }
}

params <- parse_figure_args(defaults = list(out_dir = "..", width = 24, height = 4))
plot_dir <- file.path(params$out_dir, "plots")
tsv_dir  <- file.path(params$out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_tsv <- file.path(tsv_dir, "fig-3a-human.tsv")
if (!file.exists(input_tsv)) {
  stop("Missing input TSV: ", input_tsv, "\nRun figs/figure-03/code/figure3a-human.R first.")
}

# Match ordering used by the radar plot scripts
pipelines <- c(
  "WTC11_drna_ont",
  "WTC11_cdna_ont",
  "WTC11_cdna_pacbio",
  "WTC11_drna_ont_ls",
  "WTC11_cdna_ont_ls",
  "WTC11_cdna_pacbio_ls"
)

format_pipeline_name <- function(prefix) {
  parts <- strsplit(prefix, "_")[[1]]
  name_map <- c(drna = "dRNA", cdna = "cDNA", ont = "ONT", pacbio = "PacBio")
  formatted_parts <- c("ES")
  for (p in parts[-1]) {
    if (p %in% names(name_map)) {
      formatted_parts <- c(formatted_parts, name_map[p])
    } else if (p == "ls") {
      formatted_parts <- c(formatted_parts, "LS")
    } else {
      formatted_parts <- c(formatted_parts, toupper(p))
    }
  }
  paste(formatted_parts, collapse = " ")
}

metrics_labels <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
metrics_order_for_y <- rev(metrics_labels)
metric_label_width <- 0.25

raw <- read_tsv(input_tsv, show_col_types = FALSE)
metrics_df <- raw %>%
  filter(record_type == "radar_metrics") %>%
  transmute(
    pipeline = pipeline,
    metric = metric,
    sirv = as.numeric(sirv),
    tusco = as.numeric(tusco)
  )

metrics_df <- metrics_df %>%
  mutate(
    pipeline_label = vapply(pipeline, format_pipeline_name, character(1)),
    pipeline_label = factor(pipeline_label, levels = vapply(pipelines, format_pipeline_name, character(1))),
    metric = factor(metric, levels = metrics_order_for_y),
    metric_idx = as.integer(metric)
  )

seg_df <- metrics_df %>%
  filter(!is.na(sirv), !is.na(tusco)) %>%
  transmute(pipeline_label, metric_idx, x = sirv, xend = tusco, y = metric_idx, yend = metric_idx)

point_df <- metrics_df %>%
  pivot_longer(cols = c(sirv, tusco), names_to = "Type", values_to = "value") %>%
  mutate(
    Type = recode(Type, sirv = "SIRVs", tusco = "TUSCO"),
    Type = factor(Type, levels = c("SIRVs", "TUSCO"))
  )

stripe_df <- tibble(metric_idx = seq_along(metrics_order_for_y)) %>%
  filter(metric_idx %% 2 == 0) %>%
  mutate(ymin = metric_idx - 0.5, ymax = metric_idx + 0.5)

type_colors <- c("SIRVs" = "#cab2d6", "TUSCO" = "#a8d5a0")

y_offset <- 0.12
stripe_fill <- "#FAFAFA"
connector_color <- "grey90"
connector_lwd <- 0.75
point_size_primary <- 4.2
point_size_secondary <- 3.8
base_font <- 12

make_metric_label_plot <- function() {
  ggplot() +
    geom_rect(
      data = stripe_df,
      aes(ymin = ymin, ymax = ymax),
      xmin = -Inf, xmax = Inf,
      fill = stripe_fill, color = NA
    ) +
    geom_blank(aes(x = 0, y = seq_along(metrics_order_for_y))) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq_along(metrics_order_for_y),
      labels = metrics_order_for_y,
      expand = c(0.06, 0.06)
    ) +
    theme_classic(base_family = "Helvetica", base_size = base_font) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_text(size = base_font, hjust = 1, margin = margin(r = 1.0, unit = "mm")),
      plot.margin = margin(0.2, 2.4, 0.2, 2.4, unit = "mm")
    )
}

make_pipeline_plot <- function(pipeline_id) {
  d <- metrics_df %>% filter(pipeline == pipeline_id)
  if (nrow(d) == 0) {
    return(
      ggplot() +
        theme_void() +
        theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "mm"))
    )
  }

  seg_df <- d %>%
    filter(!is.na(sirv), !is.na(tusco)) %>%
    transmute(
      x = sirv, xend = tusco,
      y = metric_idx,
      yend = metric_idx
    )

  point_df <- d %>%
    pivot_longer(cols = c(sirv, tusco), names_to = "Type", values_to = "value") %>%
    mutate(
      Type = recode(Type, sirv = "SIRVs", tusco = "TUSCO"),
      Type = factor(Type, levels = c("SIRVs", "TUSCO")),
      y = metric_idx
    )

  p <- ggplot() +
    geom_rect(
      data = stripe_df,
      aes(ymin = ymin, ymax = ymax),
      xmin = -Inf, xmax = Inf,
      fill = stripe_fill, color = NA
    ) +
    geom_segment(
      data = seg_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = connector_lwd,
      color = connector_color,
      lineend = "round"
    ) +
    geom_point(
      data = point_df %>% filter(Type == "TUSCO"),
      aes(x = value, y = y, fill = Type, shape = Type),
      size = point_size_primary,
      stroke = 0.35,
      color = "white",
      show.legend = FALSE
    ) +
    geom_point(
      data = point_df %>% filter(Type == "SIRVs"),
      aes(x = value, y = y, fill = Type, shape = Type),
      size = point_size_secondary,
      stroke = 0.35,
      color = "white",
      show.legend = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(
      breaks = seq_along(metrics_order_for_y),
      labels = metrics_order_for_y,
      expand = c(0.06, 0.06)
    ) +
    scale_fill_manual(
      values = type_colors,
      breaks = c("SIRVs", "TUSCO"),
      labels = c("SIRVs" = "SIRVs", "TUSCO" = "TUSCO evaluation")
    ) +
    scale_shape_manual(
      values = c("SIRVs" = 21, "TUSCO" = 24),
      breaks = c("SIRVs", "TUSCO"),
      labels = c("SIRVs" = "SIRVs", "TUSCO" = "TUSCO evaluation")
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_family = "Helvetica", base_size = base_font) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.25),
      axis.line.x = element_line(linewidth = 0.25),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "mm")
    )

  p
}

metric_label_plot <- make_metric_label_plot()

pipeline_titles <- lapply(pipelines, function(p) {
  ggdraw() +
    draw_label(format_pipeline_name(p), fontfamily = "Helvetica", fontface = "bold", size = 14) +
    theme(plot.margin = margin(0, 2.4, 0, 2.4, unit = "mm"))
})
blank_pad <- ggdraw() + theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))

gap_rel_width <- 0.08
add_spacers <- function(plot_list) {
  out <- list()
  for (i in seq_along(plot_list)) {
    out <- c(out, list(plot_list[[i]]))
    if (i < length(plot_list)) out <- c(out, list(ggplot() + theme_void()))
  }
  out
}

pipeline_titles_gapped <- add_spacers(pipeline_titles)
pipeline_plots <- lapply(pipelines, make_pipeline_plot)
pipeline_plots <- cowplot::align_plots(plotlist = pipeline_plots, align = "hv", axis = "tblr")
pipeline_plots_gapped <- add_spacers(pipeline_plots)

pipeline_rel_widths <- c()
for (i in seq_along(pipelines)) {
  pipeline_rel_widths <- c(pipeline_rel_widths, 1)
  if (i < length(pipelines)) pipeline_rel_widths <- c(pipeline_rel_widths, gap_rel_width)
}

top_labels_row <- plot_grid(
  plotlist = c(list(blank_pad), pipeline_titles_gapped),
  ncol = length(pipeline_titles_gapped) + 1,
  rel_widths = c(metric_label_width, pipeline_rel_widths)
)

main_row <- plot_grid(
  plotlist = c(list(metric_label_plot), pipeline_plots_gapped),
  ncol = length(pipeline_plots_gapped) + 1,
  rel_widths = c(metric_label_width, pipeline_rel_widths)
)

legend_df <- data.frame(Type = c("SIRVs", "TUSCO"), x = c(1, 2), y = c(1, 2))
p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = Type, shape = Type)) +
  geom_point(size = 4.0, stroke = 0.35, color = "white") +
  scale_fill_manual(values = type_colors,
                    breaks = c("SIRVs", "TUSCO"),
                    labels = c("SIRVs" = "SIRVs", "TUSCO" = "TUSCO evaluation")) +
  scale_shape_manual(values = c("SIRVs" = 21, "TUSCO" = 24),
                     breaks = c("SIRVs", "TUSCO"),
                     labels = c("SIRVs" = "SIRVs", "TUSCO" = "TUSCO evaluation")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, family = "Helvetica"),
    legend.direction = "horizontal",
    legend.key.size = unit(0.9, "lines"),
    legend.spacing.x = unit(0.25, "cm")
  )
legend_grob <- cowplot::get_legend(p_legend)
legend_plot <- cowplot::ggdraw() + cowplot::draw_grob(legend_grob, 0, 0, 1, 1)
legend_row <- plot_grid(
  blank_pad,
  legend_plot,
  ncol = 2,
  rel_widths = c(metric_label_width, length(pipelines))
)

final_plot <- plot_grid(
  top_labels_row,
  main_row,
  legend_row,
  ncol = 1,
  rel_heights = c(0.07, 1, 0.06)
)

out_pdf <- file.path(plot_dir, "fig-3a-human-dumbbell.pdf")
ggsave(out_pdf, final_plot, width = params$width, height = params$height, units = "in", device = "pdf", limitsize = FALSE)
message("Saved plot: ", out_pdf)

# Also write a high-resolution PNG for quick preview/embedding
out_png <- file.path(plot_dir, "fig-3a-human-dumbbell.png")
png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
ggsave(
  out_png,
  final_plot,
  width = params$width,
  height = params$height,
  units = "in",
  dpi = 300,
  device = png_device,
  bg = "white",
  limitsize = FALSE
)
message("Saved plot: ", out_png)
