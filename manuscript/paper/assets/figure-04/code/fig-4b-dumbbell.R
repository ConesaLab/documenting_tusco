#!/usr/bin/env Rscript

###############################################################
# Figure 4B (Alternative) - Dumbbell metrics plot
# Reads existing fig-4b.tsv and renders a sample-resolved grid (tool x sample),
# preserving the original Fig 4B ordering/layout.
# Usage: Rscript fig-4b-dumbbell.R [width_in] [height_in]
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(cowplot)
  library(grid)
})

# Localize paths to this fig-4 folder (same approach as fig-4b.R)
argv <- commandArgs(trailingOnly = FALSE)
script_path <- tryCatch({
  sub("^--file=", "", argv[grep("^--file=", argv)][1])
}, error = function(e) NA_character_)
if (is.na(script_path) || script_path == "") {
  script_path <- file.path(getwd(), "figs", "figure-04", "code", "fig-4b-dumbbell.R")
}
fig_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(fig_dir, "plots")
tsv_dir <- file.path(fig_dir, "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
plot_width <- if (length(args) >= 1) as.numeric(args[[1]]) else NA_real_
plot_height <- if (length(args) >= 2) as.numeric(args[[2]]) else NA_real_

input_tsv <- file.path(tsv_dir, "fig-4b.tsv")
if (!file.exists(input_tsv)) {
  stop("Missing input TSV: ", input_tsv, "\nRun figs/figure-04/code/fig-4b.R first.")
}

# Helper to draw a group header with a bracket line below the title
make_group_header <- function(title_text, title_size = 9) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.75, label = title_text, family = "Helvetica", fontface = "bold", size = title_size/3) +
    annotate("segment", x = 0.05, xend = 0.95, y = 0.12, yend = 0.12, linewidth = 0.35) +
    annotate("segment", x = 0.05, xend = 0.05, y = 0.12, yend = 0.02, linewidth = 0.35) +
    annotate("segment", x = 0.95, xend = 0.95, y = 0.12, yend = 0.02, linewidth = 0.35) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
}

metrics_labels <- c("Sn", "nrPre", "1/red", "1-FDR", "PDR", "rPre")
metrics_order_for_y <- rev(metrics_labels)

raw <- read_tsv(input_tsv, show_col_types = FALSE)

# Ordered sample prefixes (columns) – match fig-4b.R
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

# Bottom column labels (generic, repeated for both species) – match fig-4b.R
bottom_column_labels <- c(
  "CapTrap\nPacBio", "cDNA\nPacBio", "cDNA\nONT", "CapTrap\nONT",
  "CapTrap\nPacBio", "cDNA\nPacBio", "cDNA\nONT", "CapTrap\nONT"
)

pipeline_row_specs <- list(
  list(row_label = "Bambu",            pipeline_id = "bambu"),
  list(row_label = "StringTie2",       pipeline_id = "stringtie2"),
  list(row_label = "FLAIR",            pipeline_id = "flair"),
  list(row_label = "Iso-Seq + SQ3 ML", pipeline_id = "isoseq_sqanti3ml")
)
pipeline_order <- vapply(pipeline_row_specs, function(x) x$row_label, character(1))
pipeline_id_order <- vapply(pipeline_row_specs, function(x) x$pipeline_id, character(1))

ok <- raw %>%
  filter(status == "ok") %>%
  mutate(
    pipeline_id = factor(pipeline_id, levels = pipeline_id_order),
    pipeline_label = factor(pipeline_label, levels = pipeline_order),
    sample = factor(sample, levels = overall_sample_prefixes),
    species = factor(species, levels = c("human", "mouse")),
    eval_type = factor(eval_type, levels = c("ref", "novel"))
  ) %>%
  select(pipeline_id, pipeline_label, sample, species, eval_type, Sn, nrPre, `1_red`, `1_FDR`, PDR, rPre)

long <- ok %>%
  pivot_longer(cols = c(Sn, nrPre, `1_red`, `1_FDR`, PDR, rPre), names_to = "metric_raw", values_to = "value") %>%
  mutate(
    metric = case_when(
      metric_raw == "1_red" ~ "1/red",
      metric_raw == "1_FDR" ~ "1-FDR",
      TRUE ~ metric_raw
    ),
    metric = factor(metric, levels = metrics_order_for_y),
    metric_idx = as.integer(metric),
    value = as.numeric(value)
  )

stripe_df <- tibble(metric_idx = seq_along(metrics_order_for_y)) %>%
  filter(metric_idx %% 2 == 0) %>%
  mutate(ymin = metric_idx - 0.5, ymax = metric_idx + 0.5)

shape_map <- c("ref" = 21, "novel" = 24)
eval_labels <- c("ref" = "Gencode reference annotation", "novel" = "TUSCO-novel")
species_colors <- c("human" = "#a8d5a0", "mouse" = "#1b9e77")

# Visual tuning (keep original color style, reduce background clutter)
stripe_fill <- "#FAFAFA"
connector_color <- "grey82"
connector_lwd <- 0.75
cell_margin_mm <- 0.8

# Layout widths (relative to a single sample column)
row_label_width <- 0.78
metric_label_width <- 0.55

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
    theme_classic(base_family = "Helvetica", base_size = 6.5) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_text(size = 6.5, hjust = 1, margin = margin(r = 2.4, unit = "mm")),
      plot.margin = margin(0.6, 0.6, 0.6, 0.6, unit = "mm")
    )
}

metric_label_plot <- make_metric_label_plot()

make_cell_plot <- function(pipeline_id_value, sample_prefix_value, species_value) {
  d <- long %>%
    filter(pipeline_id == pipeline_id_value, sample == sample_prefix_value, species == species_value)

  if (nrow(d) == 0) {
    return(
      ggplot() +
        geom_blank(
          data = data.frame(x = c(0, 100), y = c(1, length(metrics_order_for_y))),
          aes(x = x, y = y)
        ) +
        scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0.02, 0.02)) +
        scale_y_continuous(breaks = seq_along(metrics_order_for_y), expand = c(0.06, 0.06)) +
        theme_void() +
        theme(plot.margin = margin(cell_margin_mm, cell_margin_mm, cell_margin_mm, cell_margin_mm, unit = "mm"))
    )
  }

  d_seg <- d %>%
    select(metric_idx, eval_type, value) %>%
    tidyr::pivot_wider(names_from = eval_type, values_from = value) %>%
    filter(!is.na(ref) & !is.na(novel)) %>%
    transmute(
      x = ref, xend = novel,
      y = metric_idx, yend = metric_idx
    )

  d_points <- d

  species_color <- species_colors[[as.character(species_value)]]

  p <- ggplot() +
    geom_rect(
      data = stripe_df,
      aes(ymin = ymin, ymax = ymax),
      xmin = -Inf, xmax = Inf,
      fill = stripe_fill, color = NA
    ) +
    geom_segment(
      data = d_seg,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = connector_lwd,
      color = connector_color,
      lineend = "round"
    ) +
    geom_vline(
      xintercept = c(25, 50, 75),
      linewidth = 0.2,
      color = "grey96"
    ) +
    geom_point(
      data = d_points %>% filter(eval_type == "novel"),
      aes(x = value, y = metric_idx, shape = eval_type),
      size = 2.05,
      stroke = 0.35,
      color = "white",
      fill = species_color,
      show.legend = FALSE
    ) +
    geom_point(
      data = d_points %>% filter(eval_type == "ref"),
      aes(x = value, y = metric_idx, shape = eval_type),
      size = 1.8,
      stroke = 0.35,
      color = "white",
      fill = species_color,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = shape_map) +
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
    coord_cartesian(clip = "off") +
    theme_classic(base_family = "Helvetica", base_size = 6.5) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      axis.title = element_blank(),
      plot.margin = margin(cell_margin_mm, cell_margin_mm, cell_margin_mm, cell_margin_mm, unit = "mm")
    ) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

# ------------------------------------------------------------------
# Build a manual legend (shape = annotation; color = species) with equal spacing
# ------------------------------------------------------------------
make_legend_item <- function(shape, fill, color, label, size = 3.6, text_size = 6.2) {
  ggplot() +
    geom_point(aes(x = 0.18, y = 0.5), shape = shape, fill = fill, color = color, size = size, stroke = 0.7) +
    geom_text(
      aes(x = 0.30, y = 0.5, label = label),
      hjust = 0,
      family = "Helvetica",
      size = text_size / 2.845,
      lineheight = 0.95
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void()
}

legend_plot <- plot_grid(
  make_legend_item(21, "white", "black", "Gencode reference annotation", size = 3.6),
  make_legend_item(24, "white", "black", "TUSCO-novel annotation", size = 3.6),
  make_legend_item(22, species_colors[["human"]], "grey35", "Human", size = 3.9),
  make_legend_item(22, species_colors[["mouse"]], "grey35", "Mouse", size = 3.9),
  ncol = 4,
  rel_widths = c(1.15, 1.15, 0.6, 0.6)
)

legend_box <- ggdraw() +
  draw_grob(
    grid::rectGrob(
      x = 0.5, y = 0.5, width = 0.965, height = 0.52,
      gp = grid::gpar(fill = "white", col = "grey60", lwd = 0.7)
    )
  ) +
  draw_plot(legend_plot, x = 0.04, y = 0.26, width = 0.92, height = 0.46)

# ------------------------------------------------------------------
# Assemble sample-resolved grid (same layout cues as fig-4b.R)
# ------------------------------------------------------------------
blank_row_label <- ggdraw()
blank_metric_label <- ggdraw()

make_column_axis_plot <- function() {
  ggplot() +
    geom_blank(aes(x = 0, y = 0)) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)) +
    theme_classic(base_family = "Helvetica", base_size = 6.5) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(0, cell_margin_mm, cell_margin_mm, cell_margin_mm, unit = "mm")
    )
}

axis_bottom_grobs <- replicate(length(overall_sample_prefixes), make_column_axis_plot(), simplify = FALSE)
axis_bottom_row <- plot_grid(
  plotlist = c(list(blank_row_label), list(blank_metric_label), axis_bottom_grobs),
  ncol = 10,
  rel_widths = c(row_label_width, metric_label_width, rep(1, 8)),
  align = "h",
  axis = "tb"
)

left_placeholder <- ggdraw()
human_header <- make_group_header("Human", title_size = 9)
mouse_header <- make_group_header("Mouse", title_size = 9)
group_header_row <- plot_grid(
  left_placeholder, human_header, mouse_header,
  ncol = 3,
  rel_widths = c(row_label_width + metric_label_width, 4, 4)
)

bottom_label_grobs <- lapply(bottom_column_labels, function(txt) {
  ggdraw() + draw_label(txt, fontface = "plain", size = 7, fontfamily = "Helvetica")
})
bottom_labels_row <- plot_grid(
  plotlist = c(list(blank_row_label), list(blank_metric_label), bottom_label_grobs),
  ncol = 10,
  rel_widths = c(row_label_width, metric_label_width, rep(1, 8)),
  align = "h",
  axis = "tb"
)

row_panels <- list()
for (i in seq_along(pipeline_row_specs)) {
  row_label <- pipeline_row_specs[[i]]$row_label
  if (row_label == "Iso-Seq + SQ3 ML") {
    row_label <- "Iso-Seq\n+\nSQ3 ML"
  }
  row_label_grob <- ggdraw() +
    draw_label(
      row_label,
      x = 0.02,
      hjust = 0,
      fontfamily = "Helvetica",
      fontface = "bold",
      size = 7,
      angle = 0
    ) +
    theme(plot.margin = unit(c(0, 1, 0, 2), "mm"))

  pipeline_id <- pipeline_row_specs[[i]]$pipeline_id

  cell_grobs <- lapply(seq_along(overall_sample_prefixes), function(j) {
    sample_prefix <- overall_sample_prefixes[[j]]
    species <- if (grepl("^wtc11", sample_prefix, ignore.case = TRUE)) "human" else "mouse"
    make_cell_plot(pipeline_id_value = pipeline_id, sample_prefix_value = sample_prefix, species_value = species)
  })

  row_panels[[i]] <- plot_grid(
    plotlist = c(list(row_label_grob), list(metric_label_plot), cell_grobs),
    ncol = 10,
    rel_widths = c(row_label_width, metric_label_width, rep(1, 8))
  )
}

overall_grid <- plot_grid(
  plotlist = c(list(group_header_row), row_panels, list(axis_bottom_row), list(bottom_labels_row)),
  ncol = 1,
  rel_heights = c(0.15, rep(1, length(row_panels)), 0.22, 0.5)
)
overall_grid_padded <- overall_grid + theme(plot.margin = unit(c(6, 2, 2, 4), "mm"))
overall_grid_with_legend <- plot_grid(overall_grid_padded, legend_box, ncol = 1, rel_heights = c(1, 0.2))

# Default size: match fig-4b.pdf width (180mm) and a slightly taller height (includes axis row)
cell_width_mm <- 180 / 8
overall_height_mm <- cell_width_mm * 5.6 + 12
default_width_in <- 180 / 25.4
default_height_in <- overall_height_mm / 25.4
if (is.na(plot_width)) plot_width <- default_width_in
if (is.na(plot_height)) plot_height <- default_height_in

out_pdf <- file.path(output_dir, "fig-4b-dumbbell.pdf")
ggsave(out_pdf, overall_grid_with_legend, width = plot_width, height = plot_height, units = "in", device = "pdf", limitsize = FALSE)
message("Saved plot: ", out_pdf)

# Also write a high-resolution PNG for quick preview/embedding
out_png <- file.path(output_dir, "fig-4b-dumbbell.png")
png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
ggsave(
  out_png,
  overall_grid_with_legend,
  width = plot_width,
  height = plot_height,
  units = "in",
  dpi = 600,
  device = png_device,
  bg = "white",
  limitsize = FALSE
)
message("Saved plot: ", out_png)
