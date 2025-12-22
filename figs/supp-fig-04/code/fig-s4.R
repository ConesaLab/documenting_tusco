#!/usr/bin/env Rscript
#
# fig-s4.R
# Supplementary Figure 4: Raincloud plot showing FN rates by platform and method
# Addresses Reviewer 4 comment: shows platform is primary driver of FN differences
#

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# Note: ggsignif is no longer used; we use manual annotations with precomputed
# Holm-adjusted p-values for consistency between logged and displayed values

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
fig_dir   <- normalizePath(file.path(script_dir, ".."))
repo_root <- normalizePath(file.path(fig_dir, "../.."))

# Inputs and outputs
csv_file <- file.path(repo_root, "data", "raw", "lrgasp", "FN_correlation_plot.csv")
tusco_human_file <- file.path(repo_root, "data", "processed", "tusco", "hsa", "tusco_human.tsv")
tusco_mouse_file <- file.path(repo_root, "data", "processed", "tusco", "mmu", "tusco_mouse.tsv")

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

log_message("[start] Figure S3 Raincloud run in ", fig_dir)

# --------------------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------------------
if (!file.exists(csv_file)) {
  log_message("[error] Missing input CSV: ", csv_file, ". Skipping plot generation.")
  quit(status = 0)
}

log_message("[info] Reading CSV: ", csv_file)
metrics_df <- tryCatch({
  read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8") %>%
    mutate(Sample = trimws(Sample))
}, error = function(e) {
  log_message("[error] Failed to read CSV: ", conditionMessage(e))
  NULL
})

if (is.null(metrics_df)) quit(status = 0)

# Count TUSCO genes for each species
count_tusco_genes <- function(tsv_path) {
  if (!file.exists(tsv_path)) return(NA_integer_)
  lines <- readr::read_lines(tsv_path, progress = FALSE)
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  length(lines)
}

n_tusco_human <- count_tusco_genes(tusco_human_file)
n_tusco_mouse <- count_tusco_genes(tusco_mouse_file)

log_message("[info] TUSCO gene counts: human=", n_tusco_human, ", mouse=", n_tusco_mouse)

# --------------------------------------------------------------------------------------
# Data preparation
# --------------------------------------------------------------------------------------
df <- metrics_df %>%
  mutate(
    # Parse sample name components
    cell_line = str_extract(Sample, "^(ES|WTC11)"),
    platform = case_when(
      grepl("PacBio", Sample, fixed = TRUE) ~ "PacBio",
      grepl("ONT", Sample, fixed = TRUE) ~ "ONT",
      TRUE ~ "Unknown"
    ),
    method = case_when(
      grepl("cDNA", Sample, fixed = TRUE) ~ "cDNA",
      grepl("dRNA", Sample, fixed = TRUE) ~ "dRNA",
      TRUE ~ "Unknown"
    ),
    short_read = ifelse(grepl("-LS$", Sample), "Yes", "No"),

    # Create grouping variable
    platform_method = paste0(platform, "-", method),

    # Calculate FN percentage
    # Human (WTC11): n_tusco_human genes, Mouse (ES): n_tusco_mouse genes
    total_tusco = ifelse(cell_line == "WTC11", n_tusco_human, n_tusco_mouse),
    FN_pct = (`False Negatives` / total_tusco) * 100
  ) %>%
  filter(!is.na(FN_pct))

log_message("[info] Processed ", nrow(df), " samples")

# Order platform_method factor for consistent display
df$platform_method <- factor(df$platform_method,
                             levels = c("PacBio-cDNA", "ONT-dRNA", "ONT-cDNA"))
df$cell_line <- factor(df$cell_line)

# --------------------------------------------------------------------------------------
# Statistical tests
# --------------------------------------------------------------------------------------

# --- Pooled analysis (all samples) ---
log_message("[stat] === POOLED ANALYSIS (all samples) ===")
kw_test <- kruskal.test(FN_pct ~ platform_method, data = df)
kw_p <- kw_test$p.value
log_message("[stat] Kruskal-Wallis test: chi-squared = ", round(kw_test$statistic, 3),
            ", df = ", kw_test$parameter, ", p = ", format(kw_p, digits = 3))

# Pairwise Wilcoxon tests with Holm correction
pairwise_tests <- pairwise.wilcox.test(df$FN_pct, df$platform_method,
                                        p.adjust.method = "holm", exact = FALSE)
log_message("[stat] Pairwise Wilcoxon tests (Holm-adjusted):")
pmat <- pairwise_tests$p.value
log_message("[stat]   P-value matrix:")
for (i in seq_len(nrow(pmat))) {
  for (j in seq_len(ncol(pmat))) {
    if (!is.na(pmat[i, j])) {
      log_message("[stat]     ", rownames(pmat)[i], " vs ", colnames(pmat)[j],
                  ": p = ", format(pmat[i, j], digits = 3))
    }
  }
}

# --- Covariate-adjusted permutation ANCOVA (Freedman-Lane) ---
log_message("[stat] === PERMUTATION ANCOVA (Freedman-Lane; cell_line) ===")
set.seed(42)
perm_n <- 10000L

permute_within <- function(x, groups) {
  out <- x
  for (g in unique(groups)) {
    idx <- which(groups == g)
    out[idx] <- sample(x[idx], length(idx), replace = FALSE)
  }
  out
}

full_fit <- lm(FN_pct ~ platform_method + cell_line, data = df)
full_drop <- drop1(full_fit, test = "F")
obs_f <- if ("platform_method" %in% rownames(full_drop)) {
  full_drop["platform_method", "F value"]
} else {
  NA_real_
}

reduced_fit <- lm(FN_pct ~ cell_line, data = df)
reduced_fitted <- fitted(reduced_fit)
reduced_resid <- residuals(reduced_fit)

if (!is.na(obs_f)) {
  perm_f <- numeric(perm_n)
  for (i in seq_len(perm_n)) {
    perm_resid <- permute_within(reduced_resid, df$cell_line)
    y_perm <- reduced_fitted + perm_resid
    fit_perm <- lm(y_perm ~ platform_method + cell_line, data = df)
    drop_perm <- drop1(fit_perm, test = "F")
    perm_f[i] <- drop_perm["platform_method", "F value"]
  }
  perm_p_platform <- (sum(perm_f >= obs_f, na.rm = TRUE) + 1) / (perm_n + 1)
  log_message("[stat]   platform_method: F = ", round(obs_f, 3),
              ", perm p = ", format(perm_p_platform, digits = 3),
              " (n_perm = ", perm_n, ")")
} else {
  perm_p_platform <- NA_real_
  log_message("[stat]   platform_method: not available")
}

cell_drop <- drop1(reduced_fit, test = "F")
cell_line_p_lm <- if ("cell_line" %in% rownames(cell_drop)) {
  cell_drop["cell_line", "Pr(>F)"]
} else {
  NA_real_
}

if (!is.na(cell_line_p_lm)) {
  log_message("[stat]   cell_line (parametric): p = ", format(cell_line_p_lm, digits = 3))
}

# --- Stratified analysis by cell line ---
log_message("[stat] === STRATIFIED ANALYSIS BY CELL LINE ===")
cell_lines <- unique(df$cell_line)
stratified_results <- list()

for (cl in cell_lines) {
  df_cl <- df %>% filter(cell_line == cl)
  log_message("[stat] --- Cell line: ", cl, " (n=", nrow(df_cl), ") ---")

  # Skip if insufficient data for analysis
  n_groups <- length(unique(df_cl$platform_method))
  if (n_groups < 2) {
    log_message("[stat]   Skipping: fewer than 2 platform-method groups")
    next
  }

  # Kruskal-Wallis for this cell line
  kw_cl <- kruskal.test(FN_pct ~ platform_method, data = df_cl)
  log_message("[stat]   Kruskal-Wallis: chi-squared = ", round(kw_cl$statistic, 3),
              ", df = ", kw_cl$parameter, ", p = ", format(kw_cl$p.value, digits = 3))

  # Pairwise Wilcoxon with Holm correction for this cell line
  if (n_groups >= 2) {
    pw_cl <- pairwise.wilcox.test(df_cl$FN_pct, df_cl$platform_method,
                                   p.adjust.method = "holm", exact = FALSE)
    pmat_cl <- pw_cl$p.value
    log_message("[stat]   Pairwise Wilcoxon (Holm-adjusted):")
    for (i in seq_len(nrow(pmat_cl))) {
      for (j in seq_len(ncol(pmat_cl))) {
        if (!is.na(pmat_cl[i, j])) {
          log_message("[stat]     ", rownames(pmat_cl)[i], " vs ", colnames(pmat_cl)[j],
                      ": p = ", format(pmat_cl[i, j], digits = 3))
        }
      }
    }
    stratified_results[[cl]] <- list(kw = kw_cl, pairwise = pw_cl, pmat = pmat_cl)
  }
}

# Build subtitle using covariate-adjusted p-value when available
perm_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

perm_label <- if (!is.na(perm_p_platform)) {
  paste0(
    "Perm ANCOVA (platform|cell) p=",
    format(perm_p_platform, digits = 2, scientific = TRUE),
    " ",
    perm_to_stars(perm_p_platform)
  )
} else {
  paste0("Kruskal-Wallis p = ", format(kw_p, digits = 2, scientific = TRUE))
}

# --------------------------------------------------------------------------------------
# Custom theme matching TUSCO style
# --------------------------------------------------------------------------------------
custom_theme <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.text  = element_text(size = rel(1), color = "black"),
    axis.title = element_text(size = rel(1.1)),
    axis.line  = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.title = element_text(size = rel(1.0), face = "bold"),
    legend.text = element_text(size = rel(0.9)),
    legend.position = "right",
    plot.subtitle = element_text(size = rel(0.9), color = "gray40")
  )

# --------------------------------------------------------------------------------------
# Raincloud plot (using base ggplot2 geoms)
# --------------------------------------------------------------------------------------
# Color palette - Official brand colors
platform_colors <- c(
  "PacBio-cDNA" = "#DF1894",   # PacBio magenta
  "ONT-dRNA"    = "#4BA3B5",   # ONT teal (lighter variant for dRNA)
  "ONT-cDNA"    = "#00789A"    # ONT teal
)

# Create numeric x positions for manual offset
df$x_numeric <- as.numeric(df$platform_method)

# Build plot with manual positioning for raincloud effect
p <- ggplot(df, aes(x = platform_method, y = FN_pct)) +

  # Half-violin (distribution) - using geom_violin with position nudge
  geom_violin(
    aes(fill = platform_method),
    trim = TRUE,
    scale = "width",
    width = 0.5,
    alpha = 0.6,
    color = NA,
    position = position_nudge(x = 0.15)
  ) +

  # Box plot (summary stats) - narrow, slightly offset
  geom_boxplot(
    aes(fill = platform_method),
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.8,
    color = "black",
    linewidth = 0.3,
    position = position_nudge(x = -0.05)
  ) +

  # Jittered points (raw data) - positioned to the left
  geom_point(
    aes(x = x_numeric - 0.2, shape = cell_line, color = short_read),
    position = position_jitter(width = 0.04, height = 0, seed = 42),
    size = 2.2,
    alpha = 0.9
  ) +

  annotate(
    "text",
    x = 2,
    y = 52.5,
    label = perm_label,
    size = 2.4,
    color = "gray30"
  ) +

  # Add significance annotations using precomputed Holm-adjusted p-values
  # This ensures consistency between logged p-values and plot annotations
  {
    # Helper function to format p-value as significance stars
    p_to_stars <- function(p) {
      if (is.na(p)) return("")
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      return("")
    }

    # Extract Holm-adjusted p-values from pmat
    # pmat structure: rows are "compared to", columns are "reference"
    p_pacbio_ont_cdna <- pmat["ONT-cDNA", "PacBio-cDNA"]
    p_pacbio_ont_drna <- pmat["ONT-dRNA", "PacBio-cDNA"]

    # Create annotation labels with Holm-adjusted p-values
    label1 <- paste0(p_to_stars(p_pacbio_ont_cdna),
                     " p=", format(p_pacbio_ont_cdna, digits = 2))
    label2 <- paste0(p_to_stars(p_pacbio_ont_drna),
                     " p=", format(p_pacbio_ont_drna, digits = 2))

    # Only show annotations if p-values are not NA
    annotations <- list()
    if (!is.na(p_pacbio_ont_cdna) && p_pacbio_ont_cdna < 0.05) {
      annotations <- c(annotations, list(
        annotate("segment", x = 1, xend = 3, y = 48, yend = 48,
                 color = "gray30", linewidth = 0.3),
        annotate("segment", x = 1, xend = 1, y = 47, yend = 48,
                 color = "gray30", linewidth = 0.3),
        annotate("segment", x = 3, xend = 3, y = 47, yend = 48,
                 color = "gray30", linewidth = 0.3),
        annotate("text", x = 2, y = 49.5, label = label1,
                 size = 2.2, color = "gray30")
      ))
    }
    if (!is.na(p_pacbio_ont_drna) && p_pacbio_ont_drna < 0.05) {
      annotations <- c(annotations, list(
        annotate("segment", x = 1, xend = 2, y = 42, yend = 42,
                 color = "gray30", linewidth = 0.3),
        annotate("segment", x = 1, xend = 1, y = 41, yend = 42,
                 color = "gray30", linewidth = 0.3),
        annotate("segment", x = 2, xend = 2, y = 41, yend = 42,
                 color = "gray30", linewidth = 0.3),
        annotate("text", x = 1.5, y = 43.5, label = label2,
                 size = 2.2, color = "gray30")
      ))
    }
    annotations
  } +

  # Styling
  scale_fill_manual(values = platform_colors, guide = "none") +
  scale_shape_manual(
    name = "Cell Line",
    values = c("ES" = 16, "WTC11" = 17)
  ) +
  scale_color_manual(
    name = "Short-read\nSupport",
    values = c("No" = "gray30", "Yes" = "#FF7F00")
  ) +

  # Labels
  labs(
    x = "Sequencing Platform & Method",
    y = "False Negative Rate (%)"
  ) +

  # Axis adjustments
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, 55), expand = expansion(mult = c(0, 0.05))) +

  custom_theme +

  # Fine-tune legend
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 3)),
    color = guide_legend(order = 2, override.aes = list(size = 3))
  )

# --------------------------------------------------------------------------------------
# Save outputs
# --------------------------------------------------------------------------------------
outfile_pdf <- file.path(plot_dir, "fig-s4.pdf")
ggsave(filename = outfile_pdf, plot = p, width = 4.5, height = 3.5, device = "pdf")
log_message("[ok] Wrote PDF: ", outfile_pdf)

# Save underlying data as TSV
outfile_tsv <- file.path(tsv_dir, "fig-s4.tsv")
tsv_df <- df %>%
  transmute(
    figure_id = "fig-s4",
    sample_id = Sample,
    cell_line = cell_line,
    platform = platform,
    method = method,
    platform_method = as.character(platform_method),
    short_read_support = short_read,
    total_tusco_genes = total_tusco,
    FN_count = `False Negatives`,
    FN_pct = round(FN_pct, 2),
    kruskal_wallis_p = kw_p,
    perm_ancova_platform_p = perm_p_platform,
    cell_line_p_lm = cell_line_p_lm
  )
readr::write_tsv(tsv_df, outfile_tsv)
log_message("[ok] Wrote TSV: ", outfile_tsv)

# Summary statistics
summary_stats <- df %>%
  group_by(platform_method) %>%
  summarise(
    n = n(),
    mean_FN_pct = round(mean(FN_pct), 1),
    sd_FN_pct = round(sd(FN_pct), 1),
    median_FN_pct = round(median(FN_pct), 1),
    .groups = "drop"
  )

log_message("[info] Summary by platform-method:")
for (i in seq_len(nrow(summary_stats))) {
  row <- summary_stats[i, ]
  log_message("[info]   ", row$platform_method, ": n=", row$n,
              ", mean=", row$mean_FN_pct, "% (SD=", row$sd_FN_pct,
              "%), median=", row$median_FN_pct, "%")
}

log_message("[done] Figure S3 Raincloud completed")
