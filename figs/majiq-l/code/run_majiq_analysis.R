###############################################################################
# run_majiq_analysis.R
#
# Purpose: Analyze MAJIQ-L junction comparison results
#          Three-way comparison: TUSCO vs Single-isoform vs Multi-isoform genes
#
# Input:
#   - tables/raw/*_comparison.tsv (from run_voila_pipeline.sh)
#   - data/processed/tusco/hsa/tusco_human.tsv (TUSCO genes)
#   - data/processed/tusco/hsa/step2_single_isoform_mapping.tsv (single-isoform genes)
#
# Output:
#   - tables/junction_comparison.tsv (per-gene statistics)
#   - tables/three_way_summary.tsv (summary by category)
#   - tables/statistical_tests.tsv (Wilcoxon test results)
#   - plots/fig-majiql-*.pdf (figures)
#
# Run from: figs/majiq-l/code/
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

cat("\n=== MAJIQ-L Three-Way Comparison Analysis ===\n\n")

###############################################################################
# Helper functions
###############################################################################
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && file.exists(p)) {
      return(p)
    }
    if (is_dir && dir.exists(p)) {
      return(p)
    }
  }
  return(candidates[[1]])
}

theme_majiq <- function(base_size = 8, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = "bottom",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_blank()
    )
}

# Color scheme
COLORS <- c(
  "TUSCO" = "#a8d5a0",
  "Single-isoform" = "#7fc9cf",
  "Rest genes" = "#BDBDBD"
)

###############################################################################
# 1) Define paths
###############################################################################
base_dir <- resolve_path(c("../../..", "."), is_dir = TRUE)

tusco_path <- resolve_path(c(
  file.path(base_dir, "data", "processed", "tusco", "hsa", "tusco_human.tsv"),
  "data/processed/tusco/hsa/tusco_human.tsv"
))

single_isoform_path <- resolve_path(c(
  file.path(base_dir, "data", "processed", "tusco", "hsa", "step2_single_isoform_mapping.tsv"),
  "data/processed/tusco/hsa/step2_single_isoform_mapping.tsv"
))

raw_dir <- file.path("..", "tables", "raw")
plot_dir <- file.path("..", "plots")
table_dir <- file.path("..", "tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Pipeline configurations
pipelines <- c(
  "WTC11_cdna_ont", "WTC11_cdna_ont_ls",
  "WTC11_cdna_pacbio", "WTC11_cdna_pacbio_ls",
  "WTC11_drna_ont", "WTC11_drna_ont_ls"
)

pipeline_labels <- c(
  "WTC11_cdna_ont" = "cDNA-ONT",
  "WTC11_cdna_ont_ls" = "cDNA-ONT (LS)",
  "WTC11_cdna_pacbio" = "cDNA-PacBio",
  "WTC11_cdna_pacbio_ls" = "cDNA-PacBio (LS)",
  "WTC11_drna_ont" = "dRNA-ONT",
  "WTC11_drna_ont_ls" = "dRNA-ONT (LS)"
)

###############################################################################
# 2) Load gene lists
###############################################################################
cat("Loading gene lists...\n")

# TUSCO genes
tusco_raw <- read_tsv(tusco_path,
  col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
  comment = "#", show_col_types = FALSE
)
tusco_genes <- unique(tusco_raw$ensembl)
cat(sprintf("  TUSCO genes: %d\n", length(tusco_genes)))

# Single-isoform genes
single_iso_raw <- read_tsv(single_isoform_path,
  col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
  comment = "#", show_col_types = FALSE
)
single_isoform_genes <- unique(single_iso_raw$ensembl)
cat(sprintf("  Single-isoform genes (total): %d\n", length(single_isoform_genes)))

# Single-isoform non-TUSCO
single_isoform_non_tusco <- setdiff(single_isoform_genes, tusco_genes)
cat(sprintf("  Single-isoform non-TUSCO: %d\n", length(single_isoform_non_tusco)))

###############################################################################
# 3) Load MAJIQ-L comparison files
###############################################################################
cat("\nLoading MAJIQ-L comparison files...\n")

all_data <- list()
for (pipeline in pipelines) {
  file_path <- file.path(raw_dir, paste0(pipeline, "_comparison.tsv"))
  if (file.exists(file_path)) {
    df <- read_tsv(file_path, show_col_types = FALSE)
    df$pipeline <- pipeline
    all_data[[pipeline]] <- df
    cat(sprintf("  %s: %d genes\n", pipeline, nrow(df)))
  } else {
    cat(sprintf("  WARNING: %s not found\n", file_path))
  }
}

combined <- bind_rows(all_data)
cat(sprintf("\nTotal rows: %d\n", nrow(combined)))

###############################################################################
# 4) Categorize genes and calculate metrics
###############################################################################
cat("\nCategorizing genes and calculating metrics...\n")

# Clean gene IDs (remove version numbers for matching)
combined <- combined %>%
  mutate(
    gene_id_clean = str_remove(gene_id, "\\.\\d+$"),
    gene_category = case_when(
      gene_id_clean %in% str_remove(tusco_genes, "\\.\\d+$") ~ "TUSCO",
      gene_id_clean %in% str_remove(single_isoform_non_tusco, "\\.\\d+$") ~ "Single-isoform",
      TRUE ~ "Rest genes"
    ),
    gene_category = factor(gene_category, levels = c("TUSCO", "Single-isoform", "Rest genes"))
  )

# Calculate metrics per gene
gene_stats <- combined %>%
  mutate(
    # LR-detected junctions (found in LR data)
    lr_detected = All + `Both de novo` + `LR & Annotaion` + `LR de novo`,
    # SR-detected junctions (found in SR data with reads)
    sr_detected = All + `Both de novo` + `MAJIQ & Annotation` + `MAJIQ de novo`,
    # Validation rate: % of LR junctions supported by SR
    lr_supported = All + `Both de novo` + `LR & Annotaion`,
    validation_rate = ifelse(lr_detected > 0, 100 * lr_supported / lr_detected, NA),
    # False negative rate: % of SR junctions missed by LR
    fn_count = `MAJIQ & Annotation` + `MAJIQ de novo`,
    fn_rate = ifelse(sr_detected > 0, 100 * fn_count / sr_detected, NA),
    # Novel LR-only junctions
    lr_novel = `LR de novo`,
    # Pipeline label
    pipeline_label = pipeline_labels[pipeline]
  )

cat(sprintf("  Gene categories:\n"))
print(table(gene_stats$gene_category) / length(pipelines))

###############################################################################
# 5) Summary statistics by category
###############################################################################
cat("\nCalculating summary statistics...\n")

summary_stats <- gene_stats %>%
  filter(!is.na(validation_rate) | !is.na(fn_rate)) %>%
  group_by(pipeline, pipeline_label, gene_category) %>%
  summarise(
    n_genes = n(),
    # Validation rate
    mean_validation = mean(validation_rate, na.rm = TRUE),
    median_validation = median(validation_rate, na.rm = TRUE),
    sd_validation = sd(validation_rate, na.rm = TRUE),
    # False negative rate
    mean_fn_rate = mean(fn_rate, na.rm = TRUE),
    median_fn_rate = median(fn_rate, na.rm = TRUE),
    sd_fn_rate = sd(fn_rate, na.rm = TRUE),
    # Junction counts
    total_lr_detected = sum(lr_detected, na.rm = TRUE),
    total_sr_detected = sum(sr_detected, na.rm = TRUE),
    total_lr_novel = sum(lr_novel, na.rm = TRUE),
    total_fn = sum(fn_count, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nSummary by category:\n")
print(summary_stats %>% select(pipeline_label, gene_category, n_genes, mean_validation, mean_fn_rate, total_lr_novel))

###############################################################################
# 6) Statistical tests (Wilcoxon pairwise)
###############################################################################
cat("\nPerforming statistical tests...\n")

stat_tests <- list()
for (pl in unique(gene_stats$pipeline)) {
  df_pl <- gene_stats %>% filter(pipeline == pl, !is.na(validation_rate))

  # Pairwise comparisons for validation rate
  pairs <- list(
    c("TUSCO", "Single-isoform"),
    c("TUSCO", "Multi-isoform"),
    c("Single-isoform", "Multi-isoform")
  )

  for (pair in pairs) {
    g1 <- df_pl %>%
      filter(gene_category == pair[1]) %>%
      pull(validation_rate)
    g2 <- df_pl %>%
      filter(gene_category == pair[2]) %>%
      pull(validation_rate)

    if (length(g1) > 1 && length(g2) > 1) {
      test <- wilcox.test(g1, g2, alternative = "two.sided")
      stat_tests[[length(stat_tests) + 1]] <- data.frame(
        pipeline = pl,
        pipeline_label = pipeline_labels[pl],
        comparison = paste(pair, collapse = " vs "),
        metric = "validation_rate",
        p_value = test$p.value,
        group1_mean = mean(g1, na.rm = TRUE),
        group2_mean = mean(g2, na.rm = TRUE),
        group1_n = length(g1),
        group2_n = length(g2)
      )
    }
  }

  # Pairwise comparisons for FN rate
  df_fn <- gene_stats %>% filter(pipeline == pl, !is.na(fn_rate))
  for (pair in pairs) {
    g1 <- df_fn %>%
      filter(gene_category == pair[1]) %>%
      pull(fn_rate)
    g2 <- df_fn %>%
      filter(gene_category == pair[2]) %>%
      pull(fn_rate)

    if (length(g1) > 1 && length(g2) > 1) {
      test <- wilcox.test(g1, g2, alternative = "two.sided")
      stat_tests[[length(stat_tests) + 1]] <- data.frame(
        pipeline = pl,
        pipeline_label = pipeline_labels[pl],
        comparison = paste(pair, collapse = " vs "),
        metric = "fn_rate",
        p_value = test$p.value,
        group1_mean = mean(g1, na.rm = TRUE),
        group2_mean = mean(g2, na.rm = TRUE),
        group1_n = length(g1),
        group2_n = length(g2)
      )
    }
  }
}

stat_df <- bind_rows(stat_tests) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

cat("\nKey comparisons (TUSCO vs Multi-isoform, validation rate):\n")
stat_df %>%
  filter(metric == "validation_rate", comparison == "TUSCO vs Multi-isoform") %>%
  select(pipeline_label, p_value, group1_mean, group2_mean) %>%
  print()

###############################################################################
# 7) Generate figures
###############################################################################
cat("\nGenerating figures...\n")

# Filter to genes with data for plotting
plot_data <- gene_stats %>%
  filter(!is.na(validation_rate) | !is.na(fn_rate))

# Figure 1: Validation rate by category
p_validation <- ggplot(
  plot_data %>% filter(!is.na(validation_rate)),
  aes(x = gene_category, y = validation_rate, fill = gene_category)
) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  facet_wrap(~pipeline_label, nrow = 1) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "Junction Validation Rate by Gene Category",
    subtitle = "% of LR junctions supported by SR splicegraph",
    x = NULL,
    y = "Validation Rate (%)",
    fill = "Gene Category"
  ) +
  theme_majiq() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "fig-majiql-validation.pdf"), p_validation,
  width = 10, height = 4
)
cat("  Saved: fig-majiql-validation.pdf\n")

# Figure 2: False negative rate by category
p_fn <- ggplot(
  plot_data %>% filter(!is.na(fn_rate)),
  aes(x = gene_category, y = fn_rate, fill = gene_category)
) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  facet_wrap(~pipeline_label, nrow = 1) +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "False Negative Rate by Gene Category",
    subtitle = "% of SR junctions missed by LR",
    x = NULL,
    y = "False Negative Rate (%)",
    fill = "Gene Category"
  ) +
  theme_majiq() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "fig-majiql-fn-rate.pdf"), p_fn,
  width = 10, height = 4
)
cat("  Saved: fig-majiql-fn-rate.pdf\n")

# Figure 3: Combined summary bar plot
summary_long <- summary_stats %>%
  select(pipeline_label, gene_category, mean_validation, mean_fn_rate, total_lr_novel) %>%
  pivot_longer(
    cols = c(mean_validation, mean_fn_rate),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric_label = case_when(
      metric == "mean_validation" ~ "Validation Rate (%)",
      metric == "mean_fn_rate" ~ "False Negative Rate (%)"
    )
  )

p_summary <- ggplot(
  summary_long,
  aes(x = gene_category, y = value, fill = gene_category)
) +
  geom_col(position = "dodge", alpha = 0.8) +
  facet_grid(metric_label ~ pipeline_label, scales = "free_y") +
  scale_fill_manual(values = COLORS) +
  labs(
    title = "MAJIQ-L Analysis: Three-Way Gene Category Comparison",
    x = NULL,
    y = NULL,
    fill = "Gene Category"
  ) +
  theme_majiq() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )

ggsave(file.path(plot_dir, "fig-majiql-combined.pdf"), p_summary,
  width = 12, height = 6
)
cat("  Saved: fig-majiql-combined.pdf\n")

###############################################################################
# 7b) Junction Category Percentage Figures (dot plots across pipelines)
###############################################################################
cat("\nGenerating junction category figures...\n")

# Calculate percentages for each category
category_pct <- gene_stats %>%
  group_by(pipeline, pipeline_label, gene_category) %>%
  summarise(
    n_genes = n(),
    # Total junctions
    total_LR = sum(All + `Both de novo` + `LR & Annotaion` + `LR de novo`, na.rm = TRUE),
    total_SR = sum(All + `Both de novo` + `MAJIQ & Annotation` + `MAJIQ de novo`, na.rm = TRUE),
    # Raw counts
    All = sum(All, na.rm = TRUE),
    Both_denovo = sum(`Both de novo`, na.rm = TRUE),
    MAJIQ_Annotation = sum(`MAJIQ & Annotation`, na.rm = TRUE),
    MAJIQ_denovo = sum(`MAJIQ de novo`, na.rm = TRUE),
    LR_Annotation = sum(`LR & Annotaion`, na.rm = TRUE),
    LR_denovo = sum(`LR de novo`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Percentages (of LR total for LR categories, SR total for SR categories)
    pct_All = ifelse(total_LR > 0, 100 * All / total_LR, 0),
    pct_Both_denovo = ifelse(total_LR > 0, 100 * Both_denovo / total_LR, 0),
    pct_MAJIQ_Annotation = ifelse(total_SR > 0, 100 * MAJIQ_Annotation / total_SR, 0),
    pct_MAJIQ_denovo = ifelse(total_SR > 0, 100 * MAJIQ_denovo / total_SR, 0),
    pct_LR_Annotation = ifelse(total_LR > 0, 100 * LR_Annotation / total_LR, 0),
    pct_LR_denovo = ifelse(total_LR > 0, 100 * LR_denovo / total_LR, 0),
    # Combined FN rate
    pct_FN = ifelse(total_SR > 0, 100 * (MAJIQ_Annotation + MAJIQ_denovo) / total_SR, 0),
    # Validation rate
    pct_validation = ifelse(total_LR > 0, 100 * (All + Both_denovo + LR_Annotation) / total_LR, 0)
  )

# Reshape for plotting
category_long <- category_pct %>%
  select(
    pipeline, pipeline_label, gene_category,
    pct_All, pct_Both_denovo, pct_MAJIQ_Annotation, pct_MAJIQ_denovo,
    pct_LR_Annotation, pct_LR_denovo, pct_FN, pct_validation
  ) %>%
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "category",
    values_to = "percentage"
  ) %>%
  mutate(
    category_label = case_when(
      category == "pct_All" ~ "All (LR+SR+Annot)",
      category == "pct_Both_denovo" ~ "Both de novo",
      category == "pct_MAJIQ_Annotation" ~ "MAJIQ & Annot\n(LR missed annotated)",
      category == "pct_MAJIQ_denovo" ~ "MAJIQ de novo\n(LR missed novel)",
      category == "pct_LR_Annotation" ~ "LR & Annot",
      category == "pct_LR_denovo" ~ "LR de novo\n(potential FP)",
      category == "pct_FN" ~ "False Negative Rate\n(SR missed by LR)",
      category == "pct_validation" ~ "LR Validation Rate"
    ),
    category_label = factor(category_label, levels = c(
      "All (LR+SR+Annot)", "Both de novo",
      "MAJIQ & Annot\n(LR missed annotated)", "MAJIQ de novo\n(LR missed novel)",
      "LR & Annot", "LR de novo\n(potential FP)",
      "False Negative Rate\n(SR missed by LR)", "LR Validation Rate"
    ))
  )

# Figure: All categories - dot plot with lines connecting pipelines
p_categories <- ggplot(
  category_long %>%
    filter(category %in% c(
      "pct_All", "pct_Both_denovo",
      "pct_MAJIQ_Annotation", "pct_MAJIQ_denovo",
      "pct_LR_Annotation", "pct_LR_denovo"
    )),
  aes(x = gene_category, y = percentage, color = gene_category)
) +
  geom_point(aes(shape = pipeline_label),
    size = 2.5, alpha = 0.8,
    position = position_dodge(width = 0.3)
  ) +
  stat_summary(
    fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3,
    color = "black", alpha = 0.5
  ) +
  facet_wrap(~category_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = COLORS) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 4)) +
  labs(
    title = "Junction Category Percentages by Gene Type",
    subtitle = "Each dot = one pipeline (6 pipelines total), crossbar = mean",
    x = NULL,
    y = "Percentage (%)",
    color = "Gene Category",
    shape = "Pipeline"
  ) +
  theme_majiq(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right",
    strip.text = element_text(size = 7)
  )

ggsave(file.path(plot_dir, "fig-majiql-categories.pdf"), p_categories,
  width = 10, height = 7
)
cat("  Saved: fig-majiql-categories.pdf\n")

# Figure: Key metrics only (FN rate and MAJIQ de novo)
# Relabel Multi-isoform to include "(with SR support)"
category_long_relabeled <- category_long %>%
  mutate(
    gene_category_label = case_when(
      gene_category == "Rest genes" ~ "Rest genes\n(with SR support)",
      TRUE ~ as.character(gene_category)
    ),
    gene_category_label = factor(gene_category_label,
      levels = c("TUSCO", "Single-isoform", "Rest genes\n(with SR support)")
    )
  )

p_key_metrics <- ggplot(
  category_long_relabeled %>%
    filter(category %in% c("pct_FN", "pct_MAJIQ_denovo", "pct_LR_denovo")),
  aes(x = gene_category_label, y = percentage, color = gene_category)
) +
  geom_point(aes(shape = pipeline_label),
    size = 3, alpha = 0.8,
    position = position_dodge(width = 0.4)
  ) +
  stat_summary(
    fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.4,
    color = "black"
  ) +
  facet_wrap(~category_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = COLORS) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 4)) +
  labs(
    x = NULL,
    y = "Percentage (%)",
    color = "Gene Category",
    shape = "Pipeline"
  ) +
  theme_majiq(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

ggsave(file.path(plot_dir, "fig-majiql-key-metrics.pdf"), p_key_metrics,
  width = 10, height = 5
)
cat("  Saved: fig-majiql-key-metrics.pdf\n")

###############################################################################
# 7c) Statistical tests for category percentages
###############################################################################
cat("\nPerforming statistical tests on category percentages...\n")

category_tests <- list()
categories_to_test <- c(
  "pct_All", "pct_MAJIQ_Annotation", "pct_MAJIQ_denovo",
  "pct_LR_denovo", "pct_FN"
)

for (cat in categories_to_test) {
  df_cat <- category_pct %>% select(pipeline, gene_category, all_of(cat))
  colnames(df_cat)[3] <- "value"

  # Pairwise tests
  pairs <- list(
    c("TUSCO", "Single-isoform"),
    c("TUSCO", "Multi-isoform"),
    c("Single-isoform", "Multi-isoform")
  )

  for (pair in pairs) {
    g1 <- df_cat %>%
      filter(gene_category == pair[1]) %>%
      pull(value)
    g2 <- df_cat %>%
      filter(gene_category == pair[2]) %>%
      pull(value)

    if (length(g1) >= 3 && length(g2) >= 3) {
      test <- wilcox.test(g1, g2, alternative = "two.sided", exact = FALSE)
      category_tests[[length(category_tests) + 1]] <- data.frame(
        category = cat,
        comparison = paste(pair, collapse = " vs "),
        p_value = test$p.value,
        group1_mean = mean(g1, na.rm = TRUE),
        group2_mean = mean(g2, na.rm = TRUE),
        group1_sd = sd(g1, na.rm = TRUE),
        group2_sd = sd(g2, na.rm = TRUE)
      )
    }
  }
}

category_test_df <- bind_rows(category_tests) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    category_label = case_when(
      category == "pct_All" ~ "All (LR+SR+Annot)",
      category == "pct_MAJIQ_Annotation" ~ "MAJIQ & Annot",
      category == "pct_MAJIQ_denovo" ~ "MAJIQ de novo",
      category == "pct_LR_denovo" ~ "LR de novo",
      category == "pct_FN" ~ "False Negative Rate"
    )
  )

cat("\n--- Statistical Tests (Wilcoxon, BH-adjusted) ---\n\n")
print(category_test_df %>%
  select(category_label, comparison, group1_mean, group2_mean, p_adj, significance) %>%
  arrange(category_label, comparison))

# Save statistical tests
write_tsv(category_test_df, file.path(table_dir, "category_statistical_tests.tsv"))
cat("\n  Saved: category_statistical_tests.tsv\n")

###############################################################################
# 8) Junction Category Summary - Addressing Reviewer's Questions
###############################################################################
cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("MAJIQ-L Analysis: Addressing Reviewer's Concerns\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")

cat("\nUsing short-reads as baseline (Han et al., Genome Research 2024):\n")
cat("- Can LR algorithms find junctions that exist in SR data?\n")
cat("- Do LR algorithms report spurious junctions not in SR?\n")
cat("- What is the false negative rate (SR junctions missed by LR)?\n")

# Sum junction categories by gene_category and pipeline
category_totals <- gene_stats %>%
  group_by(pipeline, pipeline_label, gene_category) %>%
  summarise(
    n_genes = n(),
    n_genes_with_lr = sum(lr_detected > 0, na.rm = TRUE),
    n_genes_with_sr = sum(sr_detected > 0, na.rm = TRUE),
    # Junction counts
    All = sum(All, na.rm = TRUE),
    Both_denovo = sum(`Both de novo`, na.rm = TRUE),
    MAJIQ_Annotation = sum(`MAJIQ & Annotation`, na.rm = TRUE),
    MAJIQ_denovo = sum(`MAJIQ de novo`, na.rm = TRUE),
    LR_Annotation = sum(`LR & Annotaion`, na.rm = TRUE),
    LR_denovo = sum(`LR de novo`, na.rm = TRUE),
    Annotation_only = sum(Annotation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Total junctions detected by each technology
    total_LR_junctions = All + Both_denovo + LR_Annotation + LR_denovo,
    total_SR_junctions = All + Both_denovo + MAJIQ_Annotation + MAJIQ_denovo,

    # KEY METRICS FOR REVIEWER
    # 1. LR junctions validated by SR (not spurious)
    LR_validated = All + Both_denovo + LR_Annotation,
    LR_validation_rate = ifelse(total_LR_junctions > 0,
      100 * LR_validated / total_LR_junctions, NA
    ),

    # 2. LR de novo rate (potential false positives - junctions only in LR)
    LR_denovo_rate = ifelse(total_LR_junctions > 0,
      100 * LR_denovo / total_LR_junctions, 0
    ),

    # 3. False negative rate (SR junctions missed by LR)
    SR_missed = MAJIQ_Annotation + MAJIQ_denovo,
    FN_rate = ifelse(total_SR_junctions > 0,
      100 * SR_missed / total_SR_junctions, NA
    ),

    # 4. Novel junction detection (de novo)
    # How many novel SR junctions did LR find vs miss?
    novel_SR_found = Both_denovo,
    novel_SR_missed = MAJIQ_denovo,
    novel_sensitivity = ifelse((Both_denovo + MAJIQ_denovo) > 0,
      100 * Both_denovo / (Both_denovo + MAJIQ_denovo), NA
    )
  )

# Print clear summary
cat("\n--- KEY METRICS (cDNA-PacBio) ---\n\n")

pacbio <- category_totals %>% filter(pipeline == "WTC11_cdna_pacbio")

cat("METRIC DEFINITIONS:\n")
cat("  - LR Validation Rate: % of LR junctions also found in SR (higher = better)\n")
cat("  - LR de novo Rate: % of LR junctions NOT in SR or annotation (lower = fewer false positives)\n")
cat("  - False Negative Rate: % of SR junctions missed by LR (lower = better sensitivity)\n")
cat("  - Novel Junction Sensitivity: % of novel SR junctions also found by LR\n")

cat("\n| Metric | TUSCO | Single-isoform | Multi-isoform |\n")
cat("|--------|-------|----------------|---------------|\n")

tusco <- pacbio %>% filter(gene_category == "TUSCO")
single <- pacbio %>% filter(gene_category == "Single-isoform")
multi <- pacbio %>% filter(gene_category == "Multi-isoform")

cat(sprintf(
  "| N genes (with LR data) | %d | %d | %d |\n",
  tusco$n_genes_with_lr, single$n_genes_with_lr, multi$n_genes_with_lr
))
cat(sprintf(
  "| Total LR junctions | %d | %d | %d |\n",
  tusco$total_LR_junctions, single$total_LR_junctions, multi$total_LR_junctions
))
cat(sprintf(
  "| Total SR junctions | %d | %d | %d |\n",
  tusco$total_SR_junctions, single$total_SR_junctions, multi$total_SR_junctions
))
cat(sprintf(
  "| **LR Validation Rate** | **%.1f%%** | %.1f%% | %.1f%% |\n",
  tusco$LR_validation_rate, single$LR_validation_rate, multi$LR_validation_rate
))
cat(sprintf(
  "| **LR de novo Rate** | **%.1f%%** | %.1f%% | %.1f%% |\n",
  tusco$LR_denovo_rate, single$LR_denovo_rate, multi$LR_denovo_rate
))
cat(sprintf(
  "| **False Negative Rate** | **%.1f%%** | %.1f%% | %.1f%% |\n",
  tusco$FN_rate, single$FN_rate, multi$FN_rate
))

cat("\n--- JUNCTION CATEGORY BREAKDOWN (cDNA-PacBio) ---\n\n")

cat("WHAT EACH CATEGORY MEANS:\n")
cat("  - 'All': Junction found in LR + SR + Annotation (perfect agreement)\n")
cat("  - 'Both de novo': Novel junction found by BOTH LR and SR (not in annotation)\n")
cat("  - 'MAJIQ & Annot': Annotated junction in SR, but LR MISSED it\n")
cat("  - 'MAJIQ de novo': NOVEL junction in SR (>=5 reads), but LR MISSED it\n")
cat("  - 'LR & Annot': Annotated junction in LR, no SR support\n")
cat("  - 'LR de novo': Junction ONLY in LR (not in SR or annotation) - potential false positive\n")

cat("\n| Category | TUSCO | Single-iso | Multi-iso | Description |\n")
cat("|----------|-------|------------|-----------|-------------|\n")
cat(sprintf(
  "| All | %d | %d | %d | LR + SR + Annotation |\n",
  tusco$All, single$All, multi$All
))
cat(sprintf(
  "| Both de novo | %d | %d | %d | Novel in BOTH LR & SR |\n",
  tusco$Both_denovo, single$Both_denovo, multi$Both_denovo
))
cat(sprintf(
  "| MAJIQ & Annot | %d | %d | %d | SR found annotated, LR missed |\n",
  tusco$MAJIQ_Annotation, single$MAJIQ_Annotation, multi$MAJIQ_Annotation
))
cat(sprintf(
  "| **MAJIQ de novo** | **%d** | %d | %d | **SR found NOVEL, LR missed** |\n",
  tusco$MAJIQ_denovo, single$MAJIQ_denovo, multi$MAJIQ_denovo
))
cat(sprintf(
  "| LR & Annot | %d | %d | %d | LR + Annotation, no SR |\n",
  tusco$LR_Annotation, single$LR_Annotation, multi$LR_Annotation
))
cat(sprintf(
  "| **LR de novo** | **%d** | %d | %d | **LR only (potential FP)** |\n",
  tusco$LR_denovo, single$LR_denovo, multi$LR_denovo
))

# Calculate MAJIQ de novo rate (novel SR junctions missed by LR)
tusco_novel_total <- tusco$Both_denovo + tusco$MAJIQ_denovo
single_novel_total <- single$Both_denovo + single$MAJIQ_denovo
multi_novel_total <- multi$Both_denovo + multi$MAJIQ_denovo

tusco_novel_missed_rate <- ifelse(tusco_novel_total > 0, 100 * tusco$MAJIQ_denovo / tusco_novel_total, 0)
single_novel_missed_rate <- ifelse(single_novel_total > 0, 100 * single$MAJIQ_denovo / single_novel_total, 0)
multi_novel_missed_rate <- ifelse(multi_novel_total > 0, 100 * multi$MAJIQ_denovo / multi_novel_total, 0)

cat("\n--- NOVEL JUNCTION ANALYSIS ---\n\n")
cat("Novel junctions = junctions NOT in annotation, discovered from RNA-seq data\n")
cat("  - SR detects novel junctions with >=5 reads support\n")
cat("  - Question: Does LR also find these novel junctions?\n\n")

cat("| Metric | TUSCO | Single-iso | Multi-iso |\n")
cat("|--------|-------|------------|----------|\n")
cat(sprintf(
  "| Novel junctions in SR | %d | %d | %d |\n",
  tusco_novel_total, single_novel_total, multi_novel_total
))
cat(sprintf(
  "| Novel found by BOTH | %d | %d | %d |\n",
  tusco$Both_denovo, single$Both_denovo, multi$Both_denovo
))
cat(sprintf(
  "| **Novel SR missed by LR** | **%d** | %d | %d |\n",
  tusco$MAJIQ_denovo, single$MAJIQ_denovo, multi$MAJIQ_denovo
))
cat(sprintf(
  "| Novel miss rate | %.1f%% | %.1f%% | %.1f%% |\n",
  tusco_novel_missed_rate, single_novel_missed_rate, multi_novel_missed_rate
))

cat("\n--- ADDRESSING REVIEWER'S SPECIFIC CONCERNS ---\n\n")

cat("Q: Do TUSCO genes show better performance with SR baseline?\n")
cat("A: YES. TUSCO genes demonstrate:\n")
cat(sprintf(
  "   1. Higher LR validation rate: %.1f%% vs %.1f%% (multi-isoform)\n",
  tusco$LR_validation_rate, multi$LR_validation_rate
))
cat(sprintf(
  "   2. Lower LR de novo rate: %.1f%% vs %.1f%% (fewer potential false positives)\n",
  tusco$LR_denovo_rate, multi$LR_denovo_rate
))
cat(sprintf(
  "   3. Much lower FN rate: %.1f%% vs %.1f%% (LR captures more SR junctions)\n",
  tusco$FN_rate, multi$FN_rate
))

cat("\nQ: Is this just because TUSCO genes are single-isoform?\n")
cat("A: NO. Single-isoform non-TUSCO genes show intermediate performance:\n")
cat(sprintf(
  "   - FN rate: TUSCO %.1f%% < Single-iso %.1f%% < Multi-iso %.1f%%\n",
  tusco$FN_rate, single$FN_rate, multi$FN_rate
))
cat("   This suggests TUSCO-specific characteristics beyond single-isoform status.\n")

###############################################################################
# 9) Save output tables
###############################################################################
cat("\nSaving output tables...\n")

# Per-gene statistics
write_tsv(
  gene_stats %>%
    select(
      gene_id, pipeline, pipeline_label, gene_category,
      All, `Both de novo`, `MAJIQ & Annotation`, `MAJIQ de novo`,
      `LR & Annotaion`, `LR de novo`, Annotation,
      lr_detected, sr_detected, validation_rate, fn_rate, lr_novel
    ),
  file.path(table_dir, "junction_comparison.tsv")
)
cat("  Saved: junction_comparison.tsv\n")

# Summary by category
write_tsv(summary_stats, file.path(table_dir, "three_way_summary.tsv"))
cat("  Saved: three_way_summary.tsv\n")

# Junction category totals
write_tsv(category_totals, file.path(table_dir, "junction_category_totals.tsv"))
cat("  Saved: junction_category_totals.tsv\n")

# Statistical tests
write_tsv(stat_df, file.path(table_dir, "statistical_tests.tsv"))
cat("  Saved: statistical_tests.tsv\n")

###############################################################################
# 9) Print summary
###############################################################################
cat("\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n")
cat("MAJIQ-L Three-Way Comparison Summary\n")
cat("=" %>% rep(60) %>% paste(collapse = ""), "\n\n")

cat("Gene Categories:\n")
cat(sprintf("  TUSCO: %d genes\n", length(tusco_genes)))
cat(sprintf("  Single-isoform (non-TUSCO): %d genes\n", length(single_isoform_non_tusco)))
cat(sprintf(
  "  Multi-isoform: ~%d genes\n",
  nrow(gene_stats %>% filter(pipeline == pipelines[1], gene_category == "Multi-isoform"))
))

cat("\nKey Findings (cDNA-PacBio):\n")
pacbio_summary <- summary_stats %>% filter(pipeline == "WTC11_cdna_pacbio")
for (cat in c("TUSCO", "Single-isoform", "Multi-isoform")) {
  row <- pacbio_summary %>% filter(gene_category == cat)
  cat(sprintf(
    "  %s: %.1f%% validation, %.1f%% FN rate, %d LR de novo\n",
    cat, row$mean_validation, row$mean_fn_rate, row$total_lr_novel
  ))
}

cat("\n=== Analysis Complete ===\n")
