#!/usr/bin/env Rscript

###############################################################
# Figure 3A Combined - Human + Mouse Dumbbell metrics plot
# Reads existing fig-3a-human.tsv and fig-3a-mouse.tsv and renders
# a combined facet dumbbell plot with both species stacked vertically
# Usage: Rscript fig-3a-dumbbell.R [output_dir] [width] [height]
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

# 180mm = 7.09 inches width; 45mm = 1.77 inches height
params <- parse_figure_args(defaults = list(out_dir = "..", width = 7.09, height = 1.77))
plot_dir <- file.path(params$out_dir, "plots")
tsv_dir  <- file.path(params$out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Check for input TSV files
input_tsv_human <- file.path(tsv_dir, "fig-3a-human.tsv")
input_tsv_mouse <- file.path(tsv_dir, "fig-3a-mouse.tsv")

if (!file.exists(input_tsv_human)) {
  stop("Missing input TSV: ", input_tsv_human, "\nRun figs/figure-03/code/figure3a-human.R first.")
}
if (!file.exists(input_tsv_mouse)) {
  stop("Missing input TSV: ", input_tsv_mouse, "\nRun figs/figure-03/code/figure3a-mouse.R first.")
}

# Pipeline definitions
pipelines_human <- c(
  "WTC11_drna_ont", "WTC11_cdna_ont", "WTC11_cdna_pacbio",
  "WTC11_drna_ont_ls", "WTC11_cdna_ont_ls", "WTC11_cdna_pacbio_ls"
)
pipelines_mouse <- c(
  "es_drna_ont", "es_cdna_ont", "es_cdna_pacbio",
  "es_drna_ont_ls", "es_cdna_ont_ls", "es_cdna_pacbio_ls"
)

# Format pipeline names for display
# Human uses WTC11 cell line, Mouse uses ES cell line
format_pipeline_name_human <- function(prefix) {
  parts <- strsplit(prefix, "_")[[1]]
  name_map <- c(drna = "dRNA", cdna = "cDNA", ont = "ONT", pacbio = "PacBio")
  formatted_parts <- c()  # No cell line prefix - just show method

  for (p in parts[-1]) {
    if (p %in% names(name_map)) {
      formatted_parts <- c(formatted_parts, name_map[p])
    } else if (p == "ls") {
      formatted_parts <- c(formatted_parts, "LS")
    }
  }
  paste(formatted_parts, collapse = " ")
}

format_pipeline_name_mouse <- function(prefix) {
  parts <- strsplit(prefix, "_")[[1]]
  name_map <- c(drna = "dRNA", cdna = "cDNA", ont = "ONT", pacbio = "PacBio")
  formatted_parts <- c()  # No cell line prefix - just show method

  for (p in parts[-1]) {
    if (p %in% names(name_map)) {
      formatted_parts <- c(formatted_parts, name_map[p])
    } else if (p == "ls") {
      formatted_parts <- c(formatted_parts, "LS")
    }
  }
  paste(formatted_parts, collapse = " ")
}

# Metrics setup
metrics_labels <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
metrics_order_for_y <- rev(metrics_labels)
metric_label_width <- 0.26
species_label_width <- 0.32

# Styling parameters (scaled for 7pt base font to match theme_tusco)
base_font <- 7
stripe_fill <- "#FAFAFA"
connector_color <- "grey90"
connector_lwd <- 0.5
point_size_primary <- 2.5
point_size_secondary <- 2.2

# Colors: SIRVs purple, Human TUSCO light green, Mouse TUSCO dark green
type_colors_human <- c("SIRVs" = "#cab2d6", "TUSCO" = "#a8d5a0")
type_colors_mouse <- c("SIRVs" = "#cab2d6", "TUSCO" = "#1b9e77")

# Read and process data
raw_human <- read_tsv(input_tsv_human, show_col_types = FALSE)
raw_mouse <- read_tsv(input_tsv_mouse, show_col_types = FALSE)

process_metrics <- function(raw, pipelines, format_fn) {
  raw %>%
    filter(record_type == "radar_metrics") %>%
    transmute(
      pipeline = pipeline,
      metric = metric,
      sirv = as.numeric(sirv),
      tusco = as.numeric(tusco)
    ) %>%
    mutate(
      pipeline_label = vapply(pipeline, format_fn, character(1)),
      pipeline_label = factor(pipeline_label, levels = vapply(pipelines, format_fn, character(1))),
      metric = factor(metric, levels = metrics_order_for_y),
      metric_idx = as.integer(metric)
    )
}

metrics_human <- process_metrics(raw_human, pipelines_human, format_pipeline_name_human)
metrics_mouse <- process_metrics(raw_mouse, pipelines_mouse, format_pipeline_name_mouse)

# Stripe data for alternating backgrounds
stripe_df <- tibble(metric_idx = seq_along(metrics_order_for_y)) %>%
  filter(metric_idx %% 2 == 0) %>%
  mutate(ymin = metric_idx - 0.5, ymax = metric_idx + 0.5)

# Function to create metric label plot (y-axis labels)
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
      axis.text.y = element_text(size = base_font, hjust = 1, margin = margin(r = 2, unit = "mm")),
      plot.margin = margin(0.5, 2, 0.5, 0.5, unit = "mm")
    ) +
    coord_cartesian(clip = "off")
}

# Function to create a single pipeline dumbbell plot
make_pipeline_plot <- function(metrics_df, pipeline_id, type_colors) {
  d <- metrics_df %>% filter(pipeline == pipeline_id)
  if (nrow(d) == 0) {
    return(ggplot() + theme_void() + theme(plot.margin = margin(0.5, 1.5, 0.5, 1.5, unit = "mm")))
  }

  seg_df <- d %>%
    filter(!is.na(sirv), !is.na(tusco)) %>%
    transmute(x = sirv, xend = tusco, y = metric_idx, yend = metric_idx)

  point_df <- d %>%
    pivot_longer(cols = c(sirv, tusco), names_to = "Type", values_to = "value") %>%
    mutate(
      Type = recode(Type, sirv = "SIRVs", tusco = "TUSCO"),
      Type = factor(Type, levels = c("SIRVs", "TUSCO")),
      y = metric_idx
    )

  ggplot() +
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
      stroke = 0.25,
      color = "white",
      show.legend = FALSE
    ) +
    geom_point(
      data = point_df %>% filter(Type == "SIRVs"),
      aes(x = value, y = y, fill = Type, shape = Type),
      size = point_size_secondary,
      stroke = 0.25,
      color = "white",
      show.legend = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 50, 100),
      expand = c(0.08, 0.08)
    ) +
    scale_y_continuous(
      breaks = seq_along(metrics_order_for_y),
      labels = metrics_order_for_y,
      expand = c(0.06, 0.06)
    ) +
    scale_fill_manual(values = type_colors) +
    scale_shape_manual(values = c("SIRVs" = 21, "TUSCO" = 24)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_family = "Helvetica", base_size = base_font) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.2),
      axis.line.x = element_line(linewidth = 0.2),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0.5, 1.5, 0.5, 1.5, unit = "mm")
    )
}

# Helper to add spacers between plots
gap_rel_width <- 0.06
add_spacers <- function(plot_list) {
  out <- list()
  for (i in seq_along(plot_list)) {
    out <- c(out, list(plot_list[[i]]))
    if (i < length(plot_list)) out <- c(out, list(ggplot() + theme_void()))
  }
  out
}

# Build relative widths for pipeline columns
pipeline_rel_widths <- c()
for (i in seq_along(pipelines_human)) {
  pipeline_rel_widths <- c(pipeline_rel_widths, 1)
  if (i < length(pipelines_human)) pipeline_rel_widths <- c(pipeline_rel_widths, gap_rel_width)
}

# Create pipeline title labels (shared across species)
pipeline_display_labels <- vapply(pipelines_human, format_pipeline_name_human, character(1))
pipeline_titles <- lapply(pipeline_display_labels, function(lab) {
  ggdraw() +
    draw_label(lab, fontfamily = "Helvetica", fontface = "bold", size = base_font) +
    theme(plot.margin = margin(0, 1.2, 0, 1.2, unit = "mm"))
})
pipeline_titles_gapped <- add_spacers(pipeline_titles)

blank_pad <- ggdraw() + theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))

# Top row: pipeline titles
top_labels_row <- plot_grid(
  plotlist = c(list(blank_pad), list(blank_pad), pipeline_titles_gapped),
  ncol = length(pipeline_titles_gapped) + 2,
  rel_widths = c(species_label_width, metric_label_width, pipeline_rel_widths)
)

# Human row
metric_label_plot <- make_metric_label_plot()
human_plots <- lapply(pipelines_human, function(p) make_pipeline_plot(metrics_human, p, type_colors_human))
human_plots <- cowplot::align_plots(plotlist = human_plots, align = "hv", axis = "tblr")
human_plots_gapped <- add_spacers(human_plots)

human_species_label <- ggdraw() +
  draw_label("WTC11", fontfamily = "Helvetica", fontface = "bold", size = base_font, angle = 90) +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))

human_row <- plot_grid(
  plotlist = c(list(human_species_label), list(metric_label_plot), human_plots_gapped),
  ncol = length(human_plots_gapped) + 2,
  rel_widths = c(species_label_width, metric_label_width, pipeline_rel_widths)
)

# Mouse row
mouse_plots <- lapply(pipelines_mouse, function(p) make_pipeline_plot(metrics_mouse, p, type_colors_mouse))
mouse_plots <- cowplot::align_plots(plotlist = mouse_plots, align = "hv", axis = "tblr")
mouse_plots_gapped <- add_spacers(mouse_plots)

mouse_species_label <- ggdraw() +
  draw_label("ES", fontfamily = "Helvetica", fontface = "bold", size = base_font, angle = 90) +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))

mouse_row <- plot_grid(
  plotlist = c(list(mouse_species_label), list(metric_label_plot), mouse_plots_gapped),
  ncol = length(mouse_plots_gapped) + 2,
  rel_widths = c(species_label_width, metric_label_width, pipeline_rel_widths)
)

# Legend (shared)
legend_df <- data.frame(
  Type = c("SIRVs", "TUSCO Human", "TUSCO Mouse"),
  x = c(1, 2, 3),
  y = c(1, 2, 3)
)
p_legend <- ggplot(legend_df, aes(x = x, y = y, fill = Type, shape = Type)) +
  geom_point(size = 2.5, stroke = 0.25, color = "white") +
  scale_fill_manual(
    values = c("SIRVs" = "#cab2d6", "TUSCO Human" = "#a8d5a0", "TUSCO Mouse" = "#1b9e77"),
    breaks = c("SIRVs", "TUSCO Human", "TUSCO Mouse"),
    labels = c("SIRVs" = "SIRVs", "TUSCO Human" = "TUSCO (Human)", "TUSCO Mouse" = "TUSCO (Mouse)")
  ) +
  scale_shape_manual(
    values = c("SIRVs" = 21, "TUSCO Human" = 24, "TUSCO Mouse" = 24),
    breaks = c("SIRVs", "TUSCO Human", "TUSCO Mouse"),
    labels = c("SIRVs" = "SIRVs", "TUSCO Human" = "TUSCO (Human)", "TUSCO Mouse" = "TUSCO (Mouse)")
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = base_font, family = "Helvetica"),
    legend.direction = "horizontal",
    legend.key.size = unit(0.6, "lines"),
    legend.spacing.x = unit(0.15, "cm")
  )
legend_grob <- cowplot::get_legend(p_legend)
legend_plot <- cowplot::ggdraw() + cowplot::draw_grob(legend_grob, 0, 0, 1, 1)

legend_row <- plot_grid(
  blank_pad,
  blank_pad,
  legend_plot,
  ncol = 3,
  rel_widths = c(species_label_width, metric_label_width, sum(pipeline_rel_widths))
)

# Final combined plot
final_plot <- plot_grid(
  top_labels_row,
  human_row,
  mouse_row,
  legend_row,
  ncol = 1,
  rel_heights = c(0.08, 1, 1, 0.1)
)

# Save outputs
out_pdf <- file.path(plot_dir, "fig-3a-dumbbell.pdf")
ggsave(out_pdf, final_plot, width = params$width, height = params$height, units = "in", device = "pdf", limitsize = FALSE)
message("Saved plot: ", out_pdf)

out_png <- file.path(plot_dir, "fig-3a-dumbbell.png")
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
