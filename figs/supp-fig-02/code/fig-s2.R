###########################
### 0. Define Variables ###
###########################

# All paths must be confined to the fig-s2 folder.
# Derive base_dir as the parent of this script (fig-s2),
# falling back to getwd() when not available.
.script_args <- commandArgs(trailingOnly = FALSE)
.script_file <- sub('^--file=', '', .script_args[grep('^--file=', .script_args)])
base_dir <- if (length(.script_file) == 1 && nzchar(.script_file)) {
  normalizePath(file.path(dirname(.script_file), ".."), mustWork = FALSE)
} else {
  getwd()
}
plot_dir <- file.path(base_dir, "plots")
tsv_dir  <- file.path(base_dir, "tables")
log_path <- file.path(base_dir, "run.log")

# Input data: prefer local ./data/, then shared ../data/ (figs/data).
local_data_root  <- file.path(base_dir, "data")
shared_data_root <- normalizePath(file.path(base_dir, "..", "data"), mustWork = FALSE)
repo_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = FALSE)
repo_data_raw <- file.path(repo_root, "data", "raw")
repo_data_processed <- file.path(repo_root, "data", "processed")

resolve_data_path <- function(...) {
  sub <- file.path(...)
  candidates <- c(
    file.path(local_data_root, sub),
    file.path(shared_data_root, sub),
    file.path(repo_data_raw, sub),
    file.path(repo_data_processed, sub)
  )
  for (p in candidates) {
    if (!is.na(p) && file.exists(p)) return(p)
  }
  return(NA_character_)
}

exists_path <- function(p) is.character(p) && length(p) == 1 && !is.na(p) && file.exists(p)

alphagenome_human_json <- resolve_data_path("alphagenome", "human_scores.json")
alphagenome_mouse_json <- resolve_data_path("alphagenome", "mouse_scores.json")
tusco_human_tsv        <- resolve_data_path("tusco", "tusco_human.tsv")
tusco_mouse_tsv        <- resolve_data_path("tusco", "tusco_mouse.tsv")
# Fallback to processed TUSCO locations by species
if (is.na(tusco_human_tsv)) {
  cand <- file.path(repo_data_processed, "tusco", "hsa", "tusco_human.tsv")
  if (file.exists(cand)) tusco_human_tsv <- cand
}
if (is.na(tusco_mouse_tsv)) {
  cand <- file.path(repo_data_processed, "tusco", "mmu", "tusco_mouse.tsv")
  if (file.exists(cand)) tusco_mouse_tsv <- cand
}

# Output (only log-scale figure requested)
output_file_log    <- file.path(plot_dir, "fig-s2-log.pdf")

# Plot parameters
params <- list(width = 6, height = 4)

# Visual style constants (match Fig. 1d)
LOG_OFFSET <- 1e-4

##################
### 1. Libraries ###
##################

# Open execution log
log_con <- file(log_path, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
cat(sprintf("[INFO] fig-s2 run started at %s\n", format(Sys.time(), tz = "UTC")))
cat(sprintf("[INFO] Working directory: %s\n", base_dir))

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

############################################
### 2. Helpers: read TUSCO lists and JSON ###
############################################

read_tusco_gene_ids <- function(tsv_file) {
  # Robustly read a simple TSV with comment lines starting with '#'
  # Expect first column to be Ensembl Gene ID
  raw_lines <- readLines(tsv_file)
  raw_lines <- raw_lines[!startsWith(raw_lines, "#")]
  if (length(raw_lines) == 0) return(character(0))
  dt <- fread(text = raw_lines, sep = "\t", header = FALSE)
  gene_ids <- dt[[1]]
  gene_ids <- gsub("\\..*$", "", gene_ids) # drop version suffix
  unique(gene_ids)
}

compute_scores_from_json <- function(json_path) {
  # Robust per-record parsing to avoid column recycling/misalignment
  records <- fromJSON(json_path, simplifyVector = FALSE)
  n <- length(records)
  gene_id <- character(n)
  expression_score <- numeric(n)
  splicing_score <- numeric(n)

  for (i in seq_len(n)) {
    rec <- records[[i]]
    gid <- rec$gene_id %||% NA_character_
    gene_id[i] <- gsub("\\..*$", "", gid)

    # Expression: median across tissues
    expr_vals <- tryCatch(as.numeric(unlist(rec$tissue_expression, use.names = FALSE)), error = function(e) numeric(0))
    expression_score[i] <- if (length(expr_vals) > 0) median(expr_vals, na.rm = TRUE) else NA_real_

    # Splicing: median across tissues (may be empty for single-exon genes)
    if (!is.null(rec$tissue_splice_ratios) && length(rec$tissue_splice_ratios) > 0) {
      splice_vals <- tryCatch(as.numeric(unlist(rec$tissue_splice_ratios, use.names = FALSE)), error = function(e) numeric(0))
      splicing_score[i] <- if (length(splice_vals) > 0) suppressWarnings(median(splice_vals, na.rm = TRUE)) else NA_real_
    } else {
      splicing_score[i] <- NA_real_
    }
  }

  data.table(
    gene_id = gene_id,
    expression_score = expression_score,
    splicing_score = splicing_score
  )
}

#############################
### 3. Load and Prepare Data ###
#############################

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tsv_dir,  showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] Ensured output dirs: %s, %s\n", plot_dir, tsv_dir))

missing_inputs <- character(0)

cat("[INFO] Reading TUSCO gene lists...\n")
if (!exists_path(tusco_human_tsv)) missing_inputs <- c(missing_inputs, "tusco_human.tsv (local or shared)")
if (!exists_path(tusco_mouse_tsv)) missing_inputs <- c(missing_inputs, "tusco_mouse.tsv (local or shared)")

if (exists_path(tusco_human_tsv)) cat("[INFO] Using tusco_human.tsv at ", tusco_human_tsv, "\n", sep = "")
if (exists_path(tusco_mouse_tsv)) cat("[INFO] Using tusco_mouse.tsv at ", tusco_mouse_tsv, "\n", sep = "")

tusco_human <- if (exists_path(tusco_human_tsv)) read_tusco_gene_ids(tusco_human_tsv) else character(0)
tusco_mouse <- if (exists_path(tusco_mouse_tsv)) read_tusco_gene_ids(tusco_mouse_tsv) else character(0)
cat("  Human TUSCO genes:", length(tusco_human), "\n")
cat("  Mouse TUSCO genes:", length(tusco_mouse), "\n\n")

cat("[INFO] Reading AlphaGenome JSONs...\n")
if (!exists_path(alphagenome_human_json)) missing_inputs <- c(missing_inputs, "alphagenome human_scores.json (local or shared)")
if (!exists_path(alphagenome_mouse_json)) missing_inputs <- c(missing_inputs, "alphagenome mouse_scores.json (local or shared)")

if (exists_path(alphagenome_human_json)) cat("[INFO] Using human_scores.json at ", alphagenome_human_json, "\n", sep = "")
if (exists_path(alphagenome_mouse_json)) cat("[INFO] Using mouse_scores.json at ", alphagenome_mouse_json, "\n", sep = "")

human_scores <- if (exists_path(alphagenome_human_json)) compute_scores_from_json(alphagenome_human_json) else data.table()
if (nrow(human_scores) > 0) human_scores[, Species := "Human"]
mouse_scores <- if (exists_path(alphagenome_mouse_json)) compute_scores_from_json(alphagenome_mouse_json) else data.table()
if (nrow(mouse_scores) > 0) mouse_scores[, Species := "Mouse"]

cat("  Human genes in JSON:", nrow(human_scores), "\n")
cat("  Mouse genes in JSON:", nrow(mouse_scores), "\n\n")

if (length(missing_inputs) > 0) {
  cat("[WARN] Missing input files detected (remaining work will be skipped):\n")
  for (p in missing_inputs) cat("  - ", p, "\n", sep = "")
  cat("[INFO] Exiting without generating PDFs/TSVs due to missing inputs.\n")
  cat(sprintf("[INFO] fig-s2 run finished at %s\n", format(Sys.time(), tz = "UTC")))
  sink(type = "message"); sink(type = "output"); close(log_con)
  quit(save = "no", status = 0)
}

# Mark groups: TUSCO vs Other (treat JSON set as single-isoform cohort)
human_scores[, Group := ifelse(gene_id %in% tusco_human, "TUSCO", "Other")]
mouse_scores[, Group := ifelse(gene_id %in% tusco_mouse, "TUSCO", "Other")]

# Combine and pivot to long for Expression and Splicing
combined <- rbindlist(list(human_scores, mouse_scores), use.names = TRUE, fill = TRUE)

long_dt <- melt(
  combined,
  id.vars = c("gene_id", "Species", "Group"),
  measure.vars = list(
    Expression = "expression_score",
    Splicing   = "splicing_score"
  ),
  variable.name = "ScoreType",
  value.name = "Score"
)

# The above "measure.vars" with list ensures two rows per gene; if it fails due to data.table version,
# fallback to manual rbind:
if (!"Score" %in% names(long_dt)) {
  long_dt <- rbindlist(list(
    combined[, .(gene_id, Species, Group, ScoreType = "Expression", Score = expression_score)],
    combined[, .(gene_id, Species, Group, ScoreType = "Splicing",   Score = splicing_score)]
  ), use.names = TRUE, fill = TRUE)
}

# Remove NA scores for plotting, but keep track of counts
cat("[INFO] Non-NA counts by type:\n")
print(na.omit(long_dt)[, .N, by = .(Species, Group, ScoreType)])
long_dt <- long_dt[!is.na(Score)]
long_dt[, Group := factor(Group, levels = c("TUSCO", "Other"))]

# Style: create PlotGroup for fill mapping
long_dt[, PlotGroup := paste0(Species, "_", Group)]

############################################
### 4. Plot (violin + boxplot; Fig 1d style) ###
############################################

custom_colors <- c(
  "Human_TUSCO" = "#A8D5A0",  # Light Green
  "Mouse_TUSCO" = "#1b9e77",  # Dark Green
  "Human_Other" = "#D1D3D4",  # Light Grey
  "Mouse_Other" = "#808080"   # Dark Grey
)

long_dt[, PlotGroup := factor(PlotGroup, levels = names(custom_colors))]

# Remove outliers: keep only 1st–99th percentile per Species × ScoreType
facet_limits <- long_dt[, .(
  qmin = quantile(Score, 0.01, na.rm = TRUE),
  qmax = quantile(Score, 0.99, na.rm = TRUE)
), by = .(Species, ScoreType)]

filtered_dt <- merge(long_dt, facet_limits, by = c("Species", "ScoreType"), all.x = TRUE)
filtered_dt <- filtered_dt[Score >= qmin & Score <= qmax]
filtered_dt[, c("qmin", "qmax") := NULL]

# Log-scale plot uses full data (no outlier removal)
plot_dt_log <- copy(long_dt)
plot_dt_log[, log_score := log10(Score + LOG_OFFSET)]

p_log <- ggplot(plot_dt_log, aes(x = Group, y = log_score, fill = PlotGroup)) +
  geom_violin(trim = TRUE, position = position_dodge(width = 0.9), alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(width = 0.9), fill = "white", alpha = 0.5, linewidth = 0.3,
               aes(group = interaction(Species, Group))) +
  facet_grid(ScoreType ~ Species, scales = "fixed") +
  scale_fill_manual(values = custom_colors, name = "Group",
                    labels = c("Human TUSCO", "Mouse TUSCO", "Human Others", "Mouse Others")) +
  scale_x_discrete(labels = c("TUSCO" = "TUSCO", "Other" = "Other")) +
  labs(x = "", y = "log10 [Score + 1e-04]") +
  theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    # Remove panel/plot backgrounds while keeping a boxed border
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = NA, colour = "black"),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    legend.position = "right"
  )

## Linear-scale version omitted per request. Only log-scale figure is generated.

############################
### 5. Save log figure + TSV ###
############################

# =============================================================================
# 7. OUTPUT GENERATION
# =============================================================================

# Save plot to PDF with standardized dimensions
ggsave(filename = output_file_log, plot = p_log, device = "pdf", 
       width = params$width, height = params$height, units = "in", dpi = 300, bg = "transparent")
message("Saved plot: ", output_file_log)

# Matching TSV with underlying data and minimal metadata.
# fig-s2-log.tsv: per-gene log-scale values + summaries
fig_log_id <- "fig-s2-log"
raw_log <- plot_dt_log[, .(figure_id = fig_log_id, panel_id = ScoreType, gene_id, Species, Group, variable = "log10_score_plus_offset", value = log_score)]
summary_log <- raw_log[, .(
  record_type = "summary",
  n = .N,
  median = median(value, na.rm = TRUE),
  mean = mean(value, na.rm = TRUE),
  sd = sd(value, na.rm = TRUE)
), by = .(figure_id, panel_id, Species, Group)]
summary_log$variable <- "log10_score_plus_offset"; summary_log$value <- NA_real_
summary_log$gene_id <- NA_character_

raw_log[, record_type := "raw"]
setcolorder(raw_log, c("record_type", "figure_id", "panel_id", "Species", "Group", "gene_id", "variable", "value"))
setcolorder(summary_log, names(raw_log))
out_tsv_log <- file.path(tsv_dir, "fig-s2-log.tsv")
fwrite(rbindlist(list(raw_log, summary_log), use.names = TRUE, fill = TRUE), file = out_tsv_log, sep = "\t")
cat("[INFO] Wrote TSV for log figure: ", out_tsv_log, "\n", sep = "")

# === EXTRACTION FOR PAPER VALUES ===
# Extract AlphaGenome score statistics for manuscript updates
# Write inside repo root to avoid escaping workspace
extraction_file <- file.path(repo_root, "manuscript", "extraction_values.tsv")

# Extract summary statistics from the summary_log data
extraction_data_list <- list()

# Get the unique combinations from summary_log
if (nrow(summary_log) > 0) {
  for (i in 1:nrow(summary_log)) {
    row <- summary_log[i]
    metric_name <- paste0(tolower(row$Species), "_", 
                         tolower(row$Group), "_", 
                         tolower(row$panel_id), "_",
                         "median")
    
    extraction_data_list[[paste0(i, "_median")]] <- data.frame(
      figure_id = "fig-s2",
      metric_name = metric_name,
      value = formatC(row$median, format="f", digits=3),
      ci_lower = NA,
      ci_upper = NA,
      notes = paste0("AlphaGenome ", row$panel_id, " score median"),
      stringsAsFactors = FALSE
    )
    
    # Also extract mean and sd for reference
    metric_name_mean <- paste0(tolower(row$Species), "_", 
                              tolower(row$Group), "_", 
                              tolower(row$panel_id), "_",
                              "mean")
                              
    extraction_data_list[[paste0(i, "_mean")]] <- data.frame(
      figure_id = "fig-s2",
      metric_name = metric_name_mean,
      value = formatC(row$mean, format="f", digits=3),
      ci_lower = NA,
      ci_upper = NA,
      notes = paste0("AlphaGenome ", row$panel_id, " score mean"),
      stringsAsFactors = FALSE
    )
    
    metric_name_sd <- paste0(tolower(row$Species), "_", 
                            tolower(row$Group), "_", 
                            tolower(row$panel_id), "_",
                            "sd")
    
    extraction_data_list[[paste0(i, "_sd")]] <- data.frame(
      figure_id = "fig-s2",
      metric_name = metric_name_sd,
      value = formatC(row$sd, format="f", digits=3),
      ci_lower = NA,
      ci_upper = NA,  
      notes = paste0("AlphaGenome ", row$panel_id, " score SD"),
      stringsAsFactors = FALSE
    )
  }
}

# Combine all extraction data
if (length(extraction_data_list) > 0) {
  extraction_data <- rbindlist(extraction_data_list, use.names = TRUE, fill = TRUE)
  
  # Write or append to file
  if (!file.exists(extraction_file)) {
    fwrite(extraction_data, extraction_file, sep = "\t", quote = FALSE)
    cat("[INFO] Created extraction file: ", extraction_file, "\n", sep = "")
  } else {
    fwrite(extraction_data, extraction_file, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)
  }
  
  cat("[INFO] Extracted ", nrow(extraction_data), " AlphaGenome metrics to: ", extraction_file, "\n", sep = "")
} else {
  cat("[INFO] No AlphaGenome data available for extraction\n")
}
# === END EXTRACTION ===

cat(sprintf("[INFO] fig-s2 run finished at %s\n", format(Sys.time(), tz = "UTC")))
sink(type = "message"); sink(type = "output"); close(log_con)
