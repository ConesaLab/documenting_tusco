#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

first_existing_dir <- function(paths) {
  for (p in paths) {
    if (!is.na(p) && !is.null(p) && dir.exists(p)) return(p)
  }
  if (length(paths) > 0) return(paths[[1]]) else return(NA_character_)
}

find_repo_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in 0:12) {
    cand <- normalizePath(file.path(cur), winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(cand, "data")) && dir.exists(file.path(cand, "reviewer_response"))) return(cand)
    cur2 <- normalizePath(file.path(cur, ".."), winslash = "/", mustWork = FALSE)
    if (identical(cur2, cur)) break
    cur <- cur2
  }
  NA_character_
}

strip_version <- function(x) {
  ifelse(is.na(x) | is.null(x), NA_character_, gsub("\\.[0-9]+$", "", as.character(x)))
}

repo_root <- find_repo_root()
if (is.na(repo_root)) stop("Could not locate repo root from: ", getwd())

# Write outputs relative to this script's folder (stable regardless of CWD).
argv <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", argv[grep("^--file=", argv)][1])
out_dir <- if (!is.na(script_path) && script_path != "") {
  normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
plots_dir <- file.path(out_dir, "plots")
tables_dir <- file.path(out_dir, "tables")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

data_base <- first_existing_dir(list(
  file.path(repo_root, "data", "raw", "lrgasp", "tusco_novel_evl")
))
if (!dir.exists(data_base)) stop("Missing data dir: ", data_base)

message("Repo: ", repo_root)
message("Data: ", data_base)
message("Out : ", out_dir)

# ---- Load TUSCO transcript lists ----
load_tusco_transcripts <- function(tsv_path) {
  if (!file.exists(tsv_path)) return(character(0))
  dt <- fread(tsv_path, sep = "\t", header = FALSE,
              col.names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
              colClasses = "character", data.table = TRUE, showProgress = FALSE)
  unique(na.omit(strip_version(dt$transcript)))
}

tusco_tx <- list(
  human = load_tusco_transcripts(file.path(repo_root, "data", "processed", "tusco", "hsa", "tusco_human.tsv")),
  mouse = load_tusco_transcripts(file.path(repo_root, "data", "processed", "tusco", "mmu", "tusco_mouse.tsv"))
)
if (length(tusco_tx$human) == 0 || length(tusco_tx$mouse) == 0) {
  stop("Missing TUSCO transcript lists under data/processed/tusco/{hsa,mmu}/.")
}

sample_to_species <- function(sample) {
  if (grepl("^wtc11_", sample, ignore.case = TRUE)) return("human")
  if (grepl("^es_", sample, ignore.case = TRUE)) return("mouse")
  NA_character_
}

# ---- 1) ISOFORM-LEVEL "coverage" proxy across pipelines ----
# Defined as called transcript exonic length / reference transcript length.
# This is *not* per-read, but is available for all pipelines from SQANTI outputs
# and directly reflects TSS/TTS shifts (and other length differences).
discover_classification_files <- function(base_dir) {
  all <- list.files(base_dir, pattern = "classification\\.txt$", recursive = TRUE, full.names = TRUE)
  all <- all[!grepl("read_stat", all, ignore.case = TRUE)]
  all <- all[file.exists(all)]
  all
}

parse_meta_from_path <- function(path, base_dir) {
  rel <- sub(paste0("^", gsub("([\\^\\$\\.|\\+\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", normalizePath(base_dir, winslash = "/", mustWork = FALSE)), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  pipeline <- if (length(parts) >= 1) parts[[1]] else NA_character_
  eval_type <- if (length(parts) >= 2) parts[[2]] else NA_character_
  sample <- if (length(parts) >= 3) parts[[3]] else NA_character_
  data.table(
    pipeline = pipeline,
    eval_type = eval_type,
    sample = sample,
    species = sample_to_species(sample),
    file = path
  )
}

read_sqanti_class <- function(file) {
  hdr <- fread(file, sep = "\t", header = TRUE, nrows = 0, showProgress = FALSE)
  cols <- names(hdr)
  want <- c("isoform","structural_category","associated_transcript","subcategory",
            "length","ref_length","ref_exons","diff_to_TSS","diff_to_TTS")
  sel <- intersect(want, cols)
  dt <- fread(file, sep = "\t", header = TRUE, select = sel, data.table = TRUE, showProgress = FALSE)
  if (!"isoform" %in% names(dt)) setnames(dt, "pbid", "isoform", skip_absent = TRUE)
  dt
}

class_files <- discover_classification_files(data_base)
meta <- rbindlist(lapply(class_files, parse_meta_from_path, base_dir = data_base), fill = TRUE)
meta <- meta[eval_type %in% c("ref_evl","novel_evl") & !is.na(species)]
if (nrow(meta) == 0) stop("No SQANTI classification files found under: ", data_base)

isoform_dt <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  m <- meta[i]
  dt <- read_sqanti_class(m$file)
  dt[, `:=`(
    pipeline = m$pipeline,
    eval_type = m$eval_type,
    sample = m$sample,
    species = m$species
  )]
  dt
}), use.names = TRUE, fill = TRUE)

setnames(isoform_dt, "isoform", "pbid")
isoform_dt[, associated_transcript_nov := strip_version(associated_transcript)]
isoform_dt[, label := fifelse(subcategory == "reference_match", "TP", "PTP")]
isoform_dt[, `:=`(
  ref_length = suppressWarnings(as.numeric(ref_length)),
  length = suppressWarnings(as.numeric(length)),
  ref_exons = suppressWarnings(as.numeric(ref_exons)),
  diff_to_TSS = suppressWarnings(as.numeric(diff_to_TSS)),
  diff_to_TTS = suppressWarnings(as.numeric(diff_to_TTS))
)]

isoform_dt <- isoform_dt[
  structural_category %in% c("full-splice_match","incomplete-splice_match") &
    !is.na(ref_length) & ref_length > 0 &
    !is.na(ref_exons) & ref_exons >= 2
]

isoform_dt <- isoform_dt[
  (species == "human" & associated_transcript_nov %in% tusco_tx$human) |
    (species == "mouse" & associated_transcript_nov %in% tusco_tx$mouse)
]

isoform_dt[, isoform_cov := length / ref_length]
isoform_dt <- isoform_dt[!is.na(isoform_cov) & is.finite(isoform_cov)]
isoform_dt[, ends_abs_sum := abs(diff_to_TSS) + abs(diff_to_TTS)]

summ_stats <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(list(n = 0L, mean = NA_real_, sd = NA_real_, median = NA_real_, q1 = NA_real_, q3 = NA_real_, iqr = NA_real_))
  q <- as.numeric(stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE, type = 7))
  list(
    n = length(x),
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    median = q[[2]],
    q1 = q[[1]],
    q3 = q[[3]],
    iqr = q[[3]] - q[[1]]
  )
}

isoform_summary <- isoform_dt[, c(
  list(metric = "isoform_cov"),
  summ_stats(isoform_cov)
), by = .(pipeline, eval_type, sample, species, label)]

isoform_summary_pdr <- isoform_dt[, c(
  list(metric = "isoform_cov"),
  summ_stats(isoform_cov)
), by = .(pipeline, eval_type, sample, species)][, label := "TP+PTP"]

isoform_summary2 <- rbindlist(list(isoform_summary, isoform_summary_pdr), use.names = TRUE, fill = TRUE)
fwrite(isoform_summary2, file.path(tables_dir, "isoform_cov_summary.tsv"), sep = "\t")

# Handy subset used for the reviewer response text: compare Iso-Seq vs Bambu
# on the same PacBio cDNA samples.
key_isoform_subset <- isoform_summary2[
  pipeline %in% c("isoseq_sq3","bambu_sq3") &
    sample %in% c("wtc11_cdna_pacbio","es_cdna_pacbio") &
    label %in% c("TP","PTP","TP+PTP")
][order(eval_type, sample, pipeline, label)]
fwrite(key_isoform_subset, file.path(tables_dir, "key_isoform_cov_isoseq_vs_bambu.tsv"), sep = "\t")

# ---- 2) PER-READ coverage for Iso-Seq only (available read->pbid map) ----
read_stat_files <- list(
  wtc11_cdna_pacbio = file.path(data_base, "wtc11_cdna_pacbio.transcriptome.read_stat.with_len.txt"),
  es_cdna_pacbio = file.path(data_base, "es_cdna_pacbio.transcriptome.read_stat.with_len.txt")
)

load_read_stats <- function(file) {
  if (!file.exists(file)) return(data.table())
  fread(file, sep = "\t", header = TRUE, select = c("id","pbid","read_length"),
        data.table = TRUE, showProgress = FALSE)
}

isoseq_class <- isoform_dt[pipeline == "isoseq_sq3", .(
  pbid,
  species,
  sample,
  eval_type,
  label,
  associated_transcript_nov,
  ref_length,
  ref_exons
)]
isoseq_class <- unique(isoseq_class[!is.na(pbid) & !is.na(ref_length) & ref_length > 0])

reads_all <- rbindlist(lapply(names(read_stat_files), function(samp) {
  rs <- load_read_stats(read_stat_files[[samp]])
  if (nrow(rs) == 0) return(data.table())
  rs[, sample := samp]
  rs[, species := sample_to_species(samp)]
  rs
}), use.names = TRUE, fill = TRUE)

reads_all[, read_length := suppressWarnings(as.numeric(read_length))]
reads_all <- reads_all[!is.na(read_length) & is.finite(read_length)]

reads_merged <- merge(
  reads_all,
  isoseq_class,
  by = c("pbid","sample","species"),
  all.x = FALSE, all.y = FALSE,
  allow.cartesian = TRUE
)
reads_merged <- reads_merged[!is.na(ref_exons) & ref_exons >= 2]
reads_merged[, f := read_length / ref_length]
reads_merged <- reads_merged[!is.na(f) & is.finite(f)]

read_summary <- reads_merged[, c(
  list(metric = "per_read_cov"),
  summ_stats(f)
), by = .(eval_type, sample, species, label)]

read_summary_pdr <- reads_merged[, c(
  list(metric = "per_read_cov"),
  summ_stats(f)
), by = .(eval_type, sample, species)][, label := "TP+PTP"]

read_summary2 <- rbindlist(list(read_summary, read_summary_pdr), use.names = TRUE, fill = TRUE)
fwrite(read_summary2, file.path(tables_dir, "isoseq_per_read_cov_summary.tsv"), sep = "\t")

# Handy subset for quoting in responses (TP/PTP/PDR on the two Iso-Seq cDNA samples)
key_read_subset <- read_summary2[
  sample %in% c("wtc11_cdna_pacbio","es_cdna_pacbio") &
    label %in% c("TP","PTP","TP+PTP")
][order(eval_type, sample, label)]
fwrite(key_read_subset, file.path(tables_dir, "key_isoseq_per_read_cov.tsv"), sep = "\t")

# Stratify when coverage is higher vs lower: by reference length bins
reads_merged[, ref_len_bin := cut(
  ref_length,
  breaks = c(0, 1000, 2000, 3000, 5000, Inf),
  include.lowest = TRUE,
  right = TRUE,
  labels = c("<=1kb","1-2kb","2-3kb","3-5kb",">5kb")
)]

read_by_lenbin <- reads_merged[, c(summ_stats(f)), by = .(species, sample, label, ref_len_bin)]
fwrite(read_by_lenbin, file.path(tables_dir, "isoseq_per_read_cov_by_ref_len_bin.tsv"), sep = "\t")

# ---- Plots ----
plot_theme <- theme_classic(base_family = "Helvetica", base_size = 9) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25)
  )

if (nrow(reads_merged) > 0) {
  reads_plot <- copy(reads_merged)
  reads_plot[, Species := ifelse(species == "human", "Human", "Mouse")]
  # sample is perfectly confounded with Species here (wtc11_* is Human, es_* is Mouse),
  # so faceting by both produces empty panels; facet by Species only.

  p_scatter <- ggplot(reads_plot, aes(x = ref_length, y = f, color = label)) +
    geom_point(alpha = 0.08, size = 0.35) +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 2.0)) +
    facet_wrap(~Species, nrow = 1) +
    scale_color_manual(values = c(TP = "#6BAED6", PTP = "#FD8D3C")) +
    labs(x = "Reference transcript exonic length (ref_length, log10)", y = "Per-read coverage f = Lread/Lexonic", color = "Label") +
    plot_theme

  ggsave(file.path(plots_dir, "isoseq_per_read_cov_vs_ref_length.pdf"), p_scatter, width = 6.8, height = 2.6, useDingbats = FALSE)

  p_box <- ggplot(reads_plot, aes(x = ref_len_bin, y = f, fill = label)) +
    geom_violin(
      aes(group = interaction(ref_len_bin, label)),
      position = position_dodge(width = 0.8),
      trim = TRUE,
      scale = "width",
      width = 0.85,
      alpha = 0.65,
      linewidth = 0.25
    ) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
    coord_cartesian(ylim = c(0, 2.0)) +
    facet_wrap(~Species, nrow = 1) +
    scale_fill_manual(values = c(TP = "#6BAED6", PTP = "#FD8D3C")) +
    labs(x = "Reference length bin", y = "Per-read coverage f", fill = "Label") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  ggsave(file.path(plots_dir, "isoseq_per_read_cov_by_ref_length_bin.pdf"), p_box, width = 6.8, height = 2.6, useDingbats = FALSE)
}

if (nrow(isoform_dt) > 0) {
  iso_plot <- copy(isoform_dt)
  iso_plot[, Species := ifelse(species == "human", "Human", "Mouse")]
  iso_plot[, pipeline := factor(pipeline, levels = c("bambu_sq3","stringtie_sq3","flair_sq3","isoseq_sq3"))]
  iso_plot <- iso_plot[!is.na(pipeline)]

  p_iso_cov <- ggplot(iso_plot, aes(x = pipeline, y = isoform_cov, fill = label)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
    coord_cartesian(ylim = c(0, 2.0)) +
    facet_grid(eval_type ~ Species) +
    scale_fill_manual(values = c(TP = "#1F77B4", PTP = "#FF7F0E")) +
    labs(x = "Pipeline", y = "Isoform length ratio (called length / ref_length)", fill = "Label") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(plots_dir, "isoform_length_ratio_by_pipeline.pdf"), p_iso_cov, width = 7.6, height = 4.6, useDingbats = FALSE)

  p_end_diff <- ggplot(iso_plot[label == "PTP" & is.finite(ends_abs_sum)], aes(x = pipeline, y = ends_abs_sum)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25, fill = "grey85") +
    scale_y_continuous(trans = "log10") +
    facet_grid(eval_type ~ Species) +
    labs(x = "Pipeline", y = "|diff_to_TSS| + |diff_to_TTS| (log10 bp)") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(plots_dir, "ptp_end_shift_by_pipeline.pdf"), p_end_diff, width = 7.6, height = 4.6, useDingbats = FALSE)
}

message("Wrote:")
message("- ", file.path(tables_dir, "isoseq_per_read_cov_summary.tsv"))
message("- ", file.path(tables_dir, "isoseq_per_read_cov_by_ref_len_bin.tsv"))
message("- ", file.path(tables_dir, "isoform_cov_summary.tsv"))
message("- ", file.path(tables_dir, "key_isoform_cov_isoseq_vs_bambu.tsv"))
message("- ", file.path(tables_dir, "key_isoseq_per_read_cov.tsv"))
message("- ", file.path(plots_dir, "isoseq_per_read_cov_vs_ref_length.pdf"))
message("- ", file.path(plots_dir, "isoseq_per_read_cov_by_ref_length_bin.pdf"))
message("- ", file.path(plots_dir, "isoform_length_ratio_by_pipeline.pdf"))
message("- ", file.path(plots_dir, "ptp_end_shift_by_pipeline.pdf"))
