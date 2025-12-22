#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# --------------------------------------------------------------------------------------
# Resolve paths robustly based on this script's location
# --------------------------------------------------------------------------------------
argv <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", argv[grep("^--file=", argv)][1])
if (is.na(script_path) || script_path == "") {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
}

FIG_DIR <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
FIGS_DIR <- normalizePath(file.path(FIG_DIR, ".."), winslash = "/", mustWork = FALSE)
REPO_DIR <- normalizePath(file.path(FIGS_DIR, ".."), winslash = "/", mustWork = FALSE)

PLOT_DIR <- file.path(FIG_DIR, "plots")
TSV_DIR <- file.path(FIG_DIR, "tables")
RUN_LOG <- file.path(FIG_DIR, "run.log")

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TSV_DIR, showWarnings = FALSE, recursive = TRUE)

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(..., collapse = ""), "\n")
  cat(msg, file = RUN_LOG, append = TRUE)
  message(paste(..., collapse = ""))
}

first_existing_dir <- function(paths) {
  for (p in paths) {
    if (!is.na(p) && !is.null(p) && dir.exists(p)) return(p)
  }
  if (length(paths) > 0) return(paths[[1]]) else return(NA_character_)
}

strip_version <- function(x) {
  ifelse(is.na(x) | is.null(x), NA_character_, gsub("\\.[0-9]+$", "", as.character(x)))
}

sample_to_species <- function(sample) {
  if (grepl("^wtc11_", sample, ignore.case = TRUE)) return("human")
  if (grepl("^es_", sample, ignore.case = TRUE)) return("mouse")
  NA_character_
}

DATA_BASE <- first_existing_dir(list(
  file.path(REPO_DIR, "data", "raw", "lrgasp", "tusco_novel_evl")
))
if (!dir.exists(DATA_BASE)) {
  log_message("[error] Missing data dir: ", DATA_BASE, " (skipping fig-s6)")
  quit(status = 0)
}

TUSCO_TSV <- list(
  human = file.path(REPO_DIR, "data", "processed", "tusco", "hsa", "tusco_human.tsv"),
  mouse = file.path(REPO_DIR, "data", "processed", "tusco", "mmu", "tusco_mouse.tsv")
)

load_tusco_transcripts <- function(tsv_path) {
  if (!file.exists(tsv_path)) return(character(0))
  dt <- fread(tsv_path, sep = "\t", header = FALSE,
              col.names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
              colClasses = "character", data.table = TRUE, showProgress = FALSE)
  unique(na.omit(strip_version(dt$transcript)))
}

tusco_tx <- list(
  human = load_tusco_transcripts(TUSCO_TSV$human),
  mouse = load_tusco_transcripts(TUSCO_TSV$mouse)
)
if (length(tusco_tx$human) == 0 || length(tusco_tx$mouse) == 0) {
  log_message("[error] Missing TUSCO transcript lists under data/processed/tusco/{hsa,mmu}/ (skipping fig-s6)")
  quit(status = 0)
}

READ_STATS_FILES <- list(
  wtc11_cdna_pacbio = file.path(DATA_BASE, "wtc11_cdna_pacbio.transcriptome.read_stat.with_len.txt"),
  es_cdna_pacbio = file.path(DATA_BASE, "es_cdna_pacbio.transcriptome.read_stat.with_len.txt")
)

load_read_stats <- function(file) {
  if (!file.exists(file)) return(data.table())
  fread(file, sep = "\t", header = TRUE, select = c("id", "pbid", "read_length"),
        data.table = TRUE, showProgress = FALSE)
}

# --------------------------------------------------------------------------------------
# Load Iso-Seq SQANTI classification (needed for ref_length + TP/PTP label)
# --------------------------------------------------------------------------------------
discover_classification_files <- function(base_dir) {
  patterns <- c("classification\\.txt$", "isoforms_classification\\.txt$")
  all <- character(0)
  # Limit search scope to Iso-Seq pipeline directories for speed.
  search_roots <- c(
    file.path(base_dir, "isoseq_sq3"),
    file.path(base_dir, "isoseq_sq3_old")
  )
  search_roots <- search_roots[dir.exists(search_roots)]
  if (length(search_roots) == 0) search_roots <- base_dir
  for (root in search_roots) {
    for (pat in patterns) {
      all <- c(all, list.files(root, pattern = pat, recursive = TRUE, full.names = TRUE))
    }
  }
  all <- unique(all)
  all <- all[!grepl("read_stat", all, ignore.case = TRUE)]
  all[file.exists(all)]
}

parse_meta_from_path <- function(path, base_dir) {
  base <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  rel <- sub(paste0("^", gsub("([\\^\\$\\.|\\+\\(\\)\\[\\]\\{\\}\\\\])", "\\\\\\1", base), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
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
  want <- c("isoform", "structural_category", "associated_transcript", "subcategory",
            "ref_length", "ref_exons")
  sel <- intersect(want, cols)
  dt <- fread(file, sep = "\t", header = TRUE, select = sel, data.table = TRUE, showProgress = FALSE)
  if (!"isoform" %in% names(dt)) setnames(dt, "pbid", "isoform", skip_absent = TRUE)
  dt
}

class_files <- discover_classification_files(DATA_BASE)
meta <- rbindlist(lapply(class_files, parse_meta_from_path, base_dir = DATA_BASE), fill = TRUE)

# Focus on Iso-Seq; allow a legacy folder name fallback if present.
meta <- meta[pipeline %in% c("isoseq_sq3", "isoseq_sq3_old") & eval_type == "ref_evl" & !is.na(species)]
if (nrow(meta) == 0) {
  log_message("[error] No Iso-Seq SQANTI classification files found under: ", DATA_BASE, " (skipping fig-s6)")
  quit(status = 0)
}

class_dt <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  m <- meta[i]
  dt <- read_sqanti_class(m$file)
  dt[, `:=`(pipeline = m$pipeline, eval_type = m$eval_type, sample = m$sample, species = m$species)]
  dt
}), use.names = TRUE, fill = TRUE)

setnames(class_dt, "isoform", "pbid")
class_dt[, associated_transcript_nov := strip_version(associated_transcript)]
class_dt[, label := fifelse(subcategory == "reference_match", "TP", "PTP")]
class_dt[, `:=`(
  ref_length = suppressWarnings(as.numeric(ref_length)),
  ref_exons = suppressWarnings(as.numeric(ref_exons))
)]

class_dt <- class_dt[
  structural_category %in% c("full-splice_match", "incomplete-splice_match") &
    !is.na(ref_length) & ref_length > 0 &
    !is.na(ref_exons) & ref_exons >= 2
]

class_dt <- class_dt[
  (species == "human" & associated_transcript_nov %in% tusco_tx$human) |
    (species == "mouse" & associated_transcript_nov %in% tusco_tx$mouse)
]

isoseq_class <- unique(class_dt[, .(pbid, sample, species, label, ref_length, ref_exons)])
if (nrow(isoseq_class) == 0) {
  log_message("[error] Iso-Seq classification data empty after filtering (skipping fig-s6)")
  quit(status = 0)
}

# --------------------------------------------------------------------------------------
# Merge reads -> isoforms and compute per-read coverage
# --------------------------------------------------------------------------------------
reads_all <- rbindlist(lapply(names(READ_STATS_FILES), function(samp) {
  rs <- load_read_stats(READ_STATS_FILES[[samp]])
  if (nrow(rs) == 0) return(data.table())
  rs[, sample := samp]
  rs[, species := sample_to_species(samp)]
  rs
}), use.names = TRUE, fill = TRUE)

if (nrow(reads_all) == 0) {
  log_message("[error] Missing Iso-Seq read_stat.with_len files (skipping fig-s6)")
  quit(status = 0)
}

reads_all[, read_length := suppressWarnings(as.numeric(read_length))]
reads_all <- reads_all[!is.na(read_length) & is.finite(read_length)]

reads_merged <- merge(
  reads_all,
  isoseq_class,
  by = c("pbid", "sample", "species"),
  all.x = FALSE, all.y = FALSE,
  allow.cartesian = TRUE
)
reads_merged[, f := read_length / ref_length]
reads_merged <- reads_merged[!is.na(f) & is.finite(f)]

if (nrow(reads_merged) == 0) {
  log_message("[error] No merged per-read records after filtering (skipping fig-s6)")
  quit(status = 0)
}

reads_merged[, ref_len_bin := cut(
  ref_length,
  breaks = c(0, 1000, 2000, 3000, 5000, Inf),
  include.lowest = TRUE,
  right = TRUE,
  labels = c("<=1kb", "1-2kb", "2-3kb", "3-5kb", ">5kb")
)]

plot_dt <- reads_merged[!is.na(ref_len_bin) & label %in% c("TP", "PTP")]
plot_dt[, Species := ifelse(species == "human", "Human", "Mouse")]
plot_dt[, Group := factor(label, levels = c("TP", "PTP"))]

n_by_bin <- plot_dt[, .(n = .N), by = .(Species, ref_len_bin, Group)]

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

plot_summary <- plot_dt[, c(summ_stats(f)), by = .(species, sample, Species, Group, ref_len_bin)]
fwrite(plot_summary, file.path(TSV_DIR, "fig-s6.tsv"), sep = "\t")

# Colors defined in figs/figure-04/code/fig-s5.R
group_colors <- c(
  TP = "#6BAED6",
  PTP = "#FD8D3C"
)

plot_theme <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.35),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.box.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
    legend.key = element_rect(fill = "white", colour = "black", linewidth = 0.25),
    legend.key.size = grid::unit(3, "mm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )

# Match the Figure S5 style: violin distributions + boxplot overlay (no points).
dodge <- position_dodge(width = 0.8)
p <- ggplot(plot_dt, aes(x = ref_len_bin, y = f, fill = Group)) +
  geom_violin(
    aes(group = interaction(ref_len_bin, Group, drop = TRUE)),
    # Use position_dodge (not dodge2) for polygons to avoid self-crossing paths
    # that can render as long diagonal "random" lines in PDF viewers.
    position = dodge,
    trim = TRUE,
    scale = "width",
    width = 0.85,
    alpha = 0.7,
    linewidth = 0,
    color = NA
  ) +
  geom_boxplot(
    aes(group = interaction(ref_len_bin, Group, drop = TRUE)),
    position = dodge,
    width = 0.18,
    linewidth = 0.25,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  geom_point(
    data = n_by_bin,
    aes(x = ref_len_bin, y = 1.95, size = n, fill = Group),
    inherit.aes = FALSE,
    position = dodge,
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.9
  ) +
  facet_wrap(~Species, nrow = 1) +
  scale_fill_manual(values = group_colors, name = "Group", labels = c("TP", "PTP"), drop = FALSE) +
  scale_size_continuous(
    name = "Reads",
    range = c(1.2, 6.5),
    trans = "sqrt"
  ) +
  coord_cartesian(ylim = c(0, 2.0)) +
  labs(
    x = "Reference transcript exonic length (bin)",
    y = expression(italic(f) == L[read] / L[exonic]),
    fill = "Group",
    size = "Reads"
  ) +
  guides(
    fill = guide_legend(
      title = "Group",
      keyheight = grid::unit(2.5, "mm"),
      keywidth = grid::unit(4, "mm"),
      order = 1
    ),
    size = guide_legend(
      override.aes = list(shape = 21, fill = "grey85", color = "black"),
      order = 2
    )
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

out_pdf <- file.path(PLOT_DIR, "fig-s6.pdf")
ggsave(out_pdf, p, width = 6.8, height = 2.5, useDingbats = FALSE)

log_message("[ok] Wrote: ", out_pdf)
log_message("[ok] Wrote: ", file.path(TSV_DIR, "fig-s6.tsv"))
