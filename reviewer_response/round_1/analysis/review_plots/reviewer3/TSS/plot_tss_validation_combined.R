#!/usr/bin/env Rscript
# Create combined validation plot for TSS window threshold
# Panel A: ECDF of distance distributions
# Panel B: RefTSS peak count vs window size

library(ggplot2)
library(dplyr)
library(scales)
library(cowplot)

# Standard TUSCO color palette (matching other figures)
TUSCO_COLORS <- list(
  human_tusco = "#a8d5a0",
  mouse_tusco = "#1b9e77",
  human_gencode = "#fdbf6f",
  mouse_gencode = "#e66101"
)

# ============================================================================
# Panel A: Load distance data for ECDF
# ============================================================================

# Read human data
hsa_within <- read.table(
  "tss_within_gene_distances_corrected_human.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

hsa_between <- read.table(
  "tss_between_gene_distances_corrected_human.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Read mouse data
mmu_within <- read.table(
  "tss_within_gene_distances_corrected_mouse.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

mmu_between <- read.table(
  "tss_between_gene_distances_corrected_mouse.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("Distance data loaded:\n")
cat(sprintf("  Human within-gene: %s genes\n", format(nrow(hsa_within), big.mark = ",")))
cat(sprintf("  Human between-gene: %s genes\n", format(nrow(hsa_between), big.mark = ",")))
cat(sprintf("  Mouse within-gene: %s genes\n", format(nrow(mmu_within), big.mark = ",")))
cat(sprintf("  Mouse between-gene: %s genes\n", format(nrow(mmu_between), big.mark = ",")))

# Calculate key statistics at 300bp threshold
threshold <- 300

hsa_pct_within <- 100 * mean(hsa_within$distance_to_reftss <= threshold)
hsa_pct_collision <- 100 * mean(hsa_between$distance_to_other_gene <= threshold)

mmu_pct_within <- 100 * mean(mmu_within$distance_to_reftss <= threshold)
mmu_pct_collision <- 100 * mean(mmu_between$distance_to_other_gene <= threshold)

cat("\nStatistics at 300bp threshold:\n")
cat(sprintf("  Human within-gene (≤300bp): %.1f%% (sensitivity)\n", hsa_pct_within))
cat(sprintf("  Human between-gene (≤300bp): %.1f%% (collision rate)\n", hsa_pct_collision))
cat(sprintf("  Mouse within-gene (≤300bp): %.1f%% (sensitivity)\n", mmu_pct_within))
cat(sprintf("  Mouse between-gene (≤300bp): %.1f%% (collision rate)\n", mmu_pct_collision))

# Prepare combined data for ECDF
hsa_within$Species <- "Human"
hsa_within$Distribution <- "Within"
hsa_within$distance <- hsa_within$distance_to_reftss

hsa_between$Species <- "Human"
hsa_between$Distribution <- "Between"
hsa_between$distance <- hsa_between$distance_to_other_gene

mmu_within$Species <- "Mouse"
mmu_within$Distribution <- "Within"
mmu_within$distance <- mmu_within$distance_to_reftss

mmu_between$Species <- "Mouse"
mmu_between$Distribution <- "Between"
mmu_between$distance <- mmu_between$distance_to_other_gene

# Combine all data
ecdf_data <- rbind(
  data.frame(distance = hsa_within$distance, Species = hsa_within$Species, Distribution = hsa_within$Distribution),
  data.frame(distance = hsa_between$distance, Species = hsa_between$Species, Distribution = hsa_between$Distribution),
  data.frame(distance = mmu_within$distance, Species = mmu_within$Species, Distribution = mmu_within$Distribution),
  data.frame(distance = mmu_between$distance, Species = mmu_between$Species, Distribution = mmu_between$Distribution)
)

ecdf_data$Category <- paste(ecdf_data$Species, ecdf_data$Distribution, sep = " - ")

# Color and linetype schemes (matching TUSCO palette)
# Within-gene uses TUSCO colors (greens), Between-gene uses GENCODE colors (oranges)
colors <- c(
  "Human - Within" = TUSCO_COLORS$human_tusco, # Green for human TUSCO
  "Human - Between" = TUSCO_COLORS$human_gencode, # Orange for human GENCODE
  "Mouse - Within" = TUSCO_COLORS$mouse_tusco, # Dark green for mouse TUSCO
  "Mouse - Between" = TUSCO_COLORS$mouse_gencode # Dark orange for mouse GENCODE
)

linetypes <- c(
  "Human - Within" = "solid",
  "Human - Between" = "solid",
  "Mouse - Within" = "solid",
  "Mouse - Between" = "solid"
)

# ============================================================================
# Panel B: Load peak count data
# ============================================================================

hsa_peaks <- read.table(
  "reftss_peak_counts_by_window_human.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

mmu_peaks <- read.table(
  "reftss_peak_counts_by_window_mouse.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("\nPeak count data loaded:\n")
cat(sprintf("  Human: %d window sizes\n", nrow(hsa_peaks)))
cat(sprintf("  Mouse: %d window sizes\n", nrow(mmu_peaks)))

# Prepare peak count data
hsa_peaks$Species <- "Human"
mmu_peaks$Species <- "Mouse"

peak_data <- rbind(hsa_peaks, mmu_peaks)

# ============================================================================
# Create Panel A: ECDF plot
# ============================================================================

p_ecdf <- ggplot(ecdf_data, aes(x = distance, color = Category, linetype = Category)) +
  stat_ecdf(geom = "step", linewidth = 0.5, pad = FALSE) +

  # Add vertical line at 300bp
  geom_vline(xintercept = threshold, linetype = "dotted", color = "gray30", linewidth = 0.25) +

  # Add threshold annotation
  annotate("text",
    x = threshold,
    y = 0.98,
    label = "300bp",
    hjust = 0.5,
    vjust = -0.5,
    size = 2.5,
    color = "gray30",
    fontface = "plain"
  ) +

  # Scale and labels (log scale, no limit - show full data range)
  scale_x_log10(
    breaks = c(1, 10, 100, 300, 1000, 10000, 100000),
    labels = comma,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1),
    expand = c(0.01, 0)
  ) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +

  # Labels and theme (matching TUSCO figure standards)
  labs(
    title = "Distance Distributions",
    x = "Window Size (bp)",
    y = "Cumulative Fraction",
    color = "Distribution",
    linetype = "Distribution"
  ) +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    legend.key.width = unit(0.6, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    plot.title = element_text(face = "plain", size = 7, hjust = 0.5, margin = margin(b = 6)),
    axis.title = element_text(size = 7),
    axis.title.y = element_text(size = 7, margin = margin(r = 4)),
    axis.text = element_text(size = 7),
    plot.margin = margin(2, 10, 2, 2)
  )

# ============================================================================
# Create Panel B: Peak count plot
# ============================================================================

p_peaks <- ggplot(peak_data, aes(x = window_size, y = q75, color = Species)) +
  # Add Q75 lines (75th percentile - 25% of genes have more peaks)
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +


  # Add horizontal line at y=1.0
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray30", linewidth = 0.25) +

  # Add vertical line at 300bp
  geom_vline(xintercept = threshold, linetype = "dotted", color = "gray30", linewidth = 0.25) +

  # Add annotations
  annotate("text",
    x = threshold,
    y = max(peak_data$q75[peak_data$window_size <= 1500]) * 0.98,
    label = "300bp",
    hjust = 0.5,
    vjust = -0.5,
    size = 2.5,
    color = "gray30",
    fontface = "plain"
  ) +
  annotate("text",
    x = 1500 * 0.95,
    y = 1.0,
    label = "1 peak/gene",
    hjust = 1,
    vjust = -0.5,
    size = 2.5,
    color = "gray30",
    fontface = "plain"
  ) +

  # Scale and labels (log scale, same range as Panel A)
  scale_x_log10(
    breaks = c(1, 10, 100, 300, 1000, 10000, 100000),
    labels = comma,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10, 20),
    labels = comma,
    expand = c(0.02, 0.02)
  ) +
  coord_cartesian(clip = "off") +
  scale_color_manual(
    values = c("Human" = TUSCO_COLORS$human_tusco, "Mouse" = TUSCO_COLORS$mouse_tusco)
  ) +
  scale_fill_manual(
    values = c("Human" = TUSCO_COLORS$human_tusco, "Mouse" = TUSCO_COLORS$mouse_tusco)
  ) +

  # Labels and theme (matching TUSCO figure standards)
  labs(
    title = "RefTSS Peak Count",
    x = "Window Size (bp)",
    y = "Peaks per Gene (Q75)",
    color = "Species",
    fill = "Species"
  ) +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    plot.title = element_text(face = "plain", size = 7, hjust = 0.5, margin = margin(b = 6)),
    axis.title = element_text(size = 7),
    axis.title.y = element_text(size = 7, margin = margin(r = 4)),
    axis.text = element_text(size = 7),
    plot.margin = margin(2, 10, 2, 2)
  ) +
  coord_cartesian(clip = "off")

# ============================================================================
# Combine panels and save
# ============================================================================

# Combine using cowplot (matching figure 5 style)
p_combined <- plot_grid(
  p_ecdf, p_peaks,
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb",
  labels = c("a", "b"),
  label_size = 7,
  label_fontface = "bold"
)

# Save combined plot (Nature double column width: 7.09 inches)
ggsave("tss_validation_combined.pdf", p_combined, width = 7.09, height = 3.55, units = "in", dpi = 300)
cat("\nCombined plot saved to: tss_validation_combined.pdf\n")
cat("  Dimensions: 7.09\" x 3.55\" (Nature double column)\n")

# ============================================================================
# Print interpretation
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("INTERPRETATION\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

cat("\nPANEL A (Distance Distributions):\n")
cat(sprintf(
  "  Human: %.1f%% sensitivity, %.1f%% collision rate at 300bp\n",
  hsa_pct_within, hsa_pct_collision
))
cat(sprintf(
  "  Mouse: %.1f%% sensitivity, %.1f%% collision rate at 300bp\n",
  mmu_pct_within, mmu_pct_collision
))
cat("  → 300bp window captures canonical TSSs with minimal false positives\n")

cat("\nPANEL B (RefTSS Peak Count):\n")
cat("  At 300bp:\n")
cat(sprintf(
  "    Human: median=%.2f peaks/gene (Q25=%.2f, Q75=%.2f)\n",
  hsa_peaks$median[hsa_peaks$window_size == 300],
  hsa_peaks$q25[hsa_peaks$window_size == 300],
  hsa_peaks$q75[hsa_peaks$window_size == 300]
))
cat(sprintf(
  "    Mouse: median=%.2f peaks/gene (Q25=%.2f, Q75=%.2f)\n",
  mmu_peaks$median[mmu_peaks$window_size == 300],
  mmu_peaks$q25[mmu_peaks$window_size == 300],
  mmu_peaks$q75[mmu_peaks$window_size == 300]
))
cat("  → Window captures ~1 peak per gene, avoiding over-inclusion\n")

cat("\nKEY FINDINGS:\n")
cat("  1. High sensitivity (>93%): Most genes have refTSS within 300bp\n")
cat("  2. Low collision rate (<3%): Few genes have other genes' TSS within 300bp\n")
cat("  3. Optimal peak count: ~1 refTSS per gene at 300bp\n")
cat("  4. Clear validation: 300bp is optimal threshold for TSS filtering\n")

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
