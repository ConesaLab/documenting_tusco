###########################
### Figure S5: Combined Per-Read Coverage with "All" Bin
### Uses original TP/PTP classification logic + length stratification
###########################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  # rtracklayer is used via rtracklayer::import only if GTFs exist
})

# Suppress data.table NSE notes for R CMD check/lint
utils::globalVariables(c(
  "exon_len", ".", "transcript_id", "associated_transcript_nov", "associated_transcript",
  "structural_category", "label", "subcategory", "pbid", "read_length",
  "exonic_length", "id", "species", "Species", "Group", "PlotGroup",
  "ref_length", "diff_to_TSS", "diff_to_TTS", "ref_len_bin",
  "record_type", "density", "scaled_density", "f_value", "bw",
  "mean", "sd", "median", "q1", "q3", "iqr",
  "lower_hinge", "upper_hinge", "lower_whisker", "upper_whisker",
  "min", "max", "n"
))

# ---- 0. Paths ----
argv <- commandArgs(trailingOnly = FALSE)
script_path <- tryCatch(
  {
    sub("^--file=", "", argv[grep("^--file=", argv)][1])
  },
  error = function(e) NA_character_
)
if (is.na(script_path) || script_path == "") {
  script_path <- file.path(getwd(), "figs", "supp-fig-05", "code", "fig-s5.R")
}
FIG_DIR <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
FIGS_DIR <- normalizePath(file.path(FIG_DIR, ".."), winslash = "/", mustWork = FALSE)
REPO_DIR <- normalizePath(file.path(FIGS_DIR, ".."), winslash = "/", mustWork = FALSE)

first_existing <- function(paths) {
  for (p in paths) {
    if (!is.na(p) && !is.null(p) && file.exists(p)) {
      return(p)
    }
  }
  if (length(paths) > 0) {
    return(paths[[1]])
  } else {
    return(NA_character_)
  }
}
first_existing_dir <- function(paths) {
  for (p in paths) {
    if (!is.na(p) && !is.null(p) && dir.exists(p)) {
      return(p)
    }
  }
  if (length(paths) > 0) {
    return(paths[[1]])
  } else {
    return(NA_character_)
  }
}

DATA_BASE <- first_existing_dir(list(
  file.path(REPO_DIR, "data", "raw", "lrgasp", "tusco_novel_evl"),
  file.path(FIG_DIR, "data", "lrgasp", "tusco_novel_evl")
))
PLOT_DIR <- file.path(FIG_DIR, "plots")
TSV_DIR <- file.path(FIG_DIR, "tables")
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)
if (!dir.exists(TSV_DIR)) dir.create(TSV_DIR, recursive = TRUE)

READ_STATS_FILES <- list(
  wtc11_cdna_pacbio = file.path(DATA_BASE, "wtc11_cdna_pacbio.transcriptome.read_stat.with_len.txt"),
  es_cdna_pacbio    = file.path(DATA_BASE, "es_cdna_pacbio.transcriptome.read_stat.with_len.txt")
)

# Search for classification files across multiple pipelines
PIPELINES <- c("isoseq_sq3", "isoseq_sq3_old", "stringtie_sq3", "flair_sq3", "bambu_sq3")
find_class_files <- function(eval_type = c("ref_evl", "novel_evl"), sample) {
  eval_type <- match.arg(eval_type)
  out <- character(0)
  for (p in PIPELINES) {
    base <- file.path(DATA_BASE, p, eval_type, sample)
    cand <- c(
      file.path(base, "sqanti3_out", "isoforms_classification.txt"),
      file.path(base, "isoforms_classification.txt")
    )
    cand <- cand[file.exists(cand)]
    if (length(cand) == 0 && dir.exists(base)) {
      extra <- list.files(base, pattern = "^isoforms_classification\\.txt$", recursive = TRUE, full.names = TRUE)
      if (length(extra)) cand <- c(cand, extra)
    }
    out <- c(out, cand)
  }
  unique(out)
}

CLASS_FILES_REF <- list(
  wtc11_cdna_pacbio = find_class_files("ref_evl", "wtc11_cdna_pacbio"),
  es_cdna_pacbio    = find_class_files("ref_evl", "es_cdna_pacbio")
)
CLASS_FILES_NOVEL <- list(
  wtc11_cdna_pacbio = find_class_files("novel_evl", "wtc11_cdna_pacbio"),
  es_cdna_pacbio    = find_class_files("novel_evl", "es_cdna_pacbio")
)

TUSCO_TSV <- list(
  human = first_existing(list(
    file.path(REPO_DIR, "data", "processed", "tusco", "hsa", "tusco_human.tsv"),
    file.path(FIG_DIR, "data", "tusco", "tusco_human.tsv")
  )),
  mouse = first_existing(list(
    file.path(REPO_DIR, "data", "processed", "tusco", "mmu", "tusco_mouse.tsv"),
    file.path(FIG_DIR, "data", "tusco", "tusco_mouse.tsv")
  ))
)

TUSCO_GTF <- list(
  human = first_existing(list(
    file.path(REPO_DIR, "data", "processed", "tusco", "hsa", "tusco_human.gtf"),
    file.path(FIG_DIR, "data", "tusco", "tusco_human.gtf")
  )),
  mouse = first_existing(list(
    file.path(REPO_DIR, "data", "processed", "tusco", "mmu", "tusco_mouse.gtf"),
    file.path(FIG_DIR, "data", "tusco", "tusco_mouse.gtf")
  ))
)

OUTPUT_PDF <- file.path(PLOT_DIR, "fig-s5.pdf")
OUTPUT_TSV <- file.path(TSV_DIR, "fig-s5.tsv")

strip_version <- function(x) {
  ifelse(is.na(x) | is.null(x), NA_character_, gsub("\\.[0-9]+$", "", as.character(x)))
}

# ---- 1. Load TUSCO transcript sets ----
load_tusco_transcripts <- function(tsv_path) {
  if (!file.exists(tsv_path)) {
    return(character(0))
  }
  dt <- fread(tsv_path,
    sep = "\t", header = FALSE,
    col.names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    colClasses = "character", data.table = TRUE, showProgress = FALSE
  )
  tx <- unique(na.omit(strip_version(dt$transcript)))
  return(tx)
}

cat("Loading TUSCO transcript lists...\n")
tusco_tx <- list(
  human = load_tusco_transcripts(TUSCO_TSV$human),
  mouse = load_tusco_transcripts(TUSCO_TSV$mouse)
)

# ---- 2. Parse GTF to compute exonic lengths per transcript ----
compute_exonic_lengths <- function(gtf_path, transcripts_of_interest) {
  if (!file.exists(gtf_path)) {
    return(integer(0))
  }
  gr <- rtracklayer::import(gtf_path)
  df <- as.data.frame(gr)
  if (!"type" %in% names(df)) {
    if (!"type" %in% names(mcols(gr))) stop("Cannot find feature type column in GTF: ", gtf_path)
    df$type <- mcols(gr)$type
  }
  if (!"transcript_id" %in% colnames(df)) {
    if ("transcript_id" %in% colnames(mcols(gr))) {
      df$transcript_id <- as.character(mcols(gr)$transcript_id)
    } else {
      stop("transcript_id not found in GTF attributes: ", gtf_path)
    }
  }
  exon_df <- df[df$type == "exon", c("seqnames", "start", "end", "transcript_id")]
  exon_df$transcript_id <- strip_version(exon_df$transcript_id)
  toi <- unique(transcripts_of_interest[!is.na(transcripts_of_interest)])
  exon_df <- exon_df[exon_df$transcript_id %in% toi, , drop = FALSE]
  if (nrow(exon_df) == 0) {
    return(integer(0))
  }
  exon_dt <- as.data.table(exon_df)
  exon_dt[, exon_len := as.integer(end - start + 1L)]
  len_map <- exon_dt[, .(exonic_length = sum(exon_len, na.rm = TRUE)), by = transcript_id]
  setnames(len_map, "transcript_id", "tx")
  lens <- len_map$exonic_length
  names(lens) <- len_map$tx
  return(lens)
}

cat("Computing exonic lengths...\n")
exonic_lengths <- list(
  human = compute_exonic_lengths(TUSCO_GTF$human, tusco_tx$human),
  mouse = compute_exonic_lengths(TUSCO_GTF$mouse, tusco_tx$mouse)
)

# ---- 2b. Compute exon counts to remove mono-exon transcripts ----
compute_exon_counts <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    return(data.table(transcript_id = character(), exon_count = integer()))
  }
  gr <- rtracklayer::import(gtf_path)
  df <- as.data.frame(gr)
  if (!"type" %in% names(df)) {
    if (!"type" %in% names(mcols(gr))) stop("Cannot find feature type column in GTF: ", gtf_path)
    df$type <- mcols(gr)$type
  }
  if (!"transcript_id" %in% colnames(df)) {
    if ("transcript_id" %in% colnames(mcols(gr))) {
      df$transcript_id <- as.character(mcols(gr)$transcript_id)
    } else {
      stop("transcript_id not found in GTF attributes: ", gtf_path)
    }
  }
  exon_df <- df[df$type == "exon", c("transcript_id"), drop = FALSE]
  if (nrow(exon_df) == 0) {
    return(data.table(transcript_id = character(), exon_count = integer()))
  }
  exon_df$transcript_id <- strip_version(exon_df$transcript_id)
  as.data.table(exon_df)[, .(exon_count = .N), by = transcript_id]
}

exon_counts <- list(
  human = compute_exon_counts(TUSCO_GTF$human),
  mouse = compute_exon_counts(TUSCO_GTF$mouse)
)
multi_exon_tx <- list(
  human = unique(exon_counts$human[exon_count >= 2, transcript_id]),
  mouse = unique(exon_counts$mouse[exon_count >= 2, transcript_id])
)

# ---- 3. Load SQANTI classification with ORIGINAL TP/PTP rules ----
# Now also returns ref_length for length binning
load_sqanti_classification <- function(paths, tusco_tx_set) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) {
    return(data.table(pbid = character(), associated_transcript_nov = character(), label = character(), ref_length = numeric()))
  }

  usecols <- c(
    "isoform", "structural_category", "associated_transcript", "subcategory",
    "ref_length", "diff_to_TSS", "diff_to_TTS", "ref_exons"
  )
  lst <- lapply(paths, function(p) fread(p, sep = "\t", select = usecols, header = TRUE, data.table = TRUE, showProgress = FALSE))
  df <- rbindlist(lst, use.names = TRUE, fill = TRUE)
  setnames(df, "isoform", "pbid")
  df[, associated_transcript_nov := strip_version(associated_transcript)]

  # Keep only mappings to TUSCO transcripts
  df <- df[associated_transcript_nov %in% tusco_tx_set]
  if (nrow(df) == 0) {
    return(df[0])
  }

  # Keep only FSM/ISM and mono-exon_by_intron_retention
  df <- df[structural_category %in% c("full-splice_match", "incomplete-splice_match") |
    subcategory == "mono-exon_by_intron_retention"]
  if (nrow(df) == 0) {
    return(df[0])
  }

  # Ensure required columns exist and are numeric
  suppressWarnings({
    if (!"ref_length" %in% names(df)) df[, ref_length := NA_real_] else df[, ref_length := as.numeric(ref_length)]
    if (!"diff_to_TSS" %in% names(df)) df[, diff_to_TSS := NA_real_] else df[, diff_to_TSS := as.numeric(diff_to_TSS)]
    if (!"diff_to_TTS" %in% names(df)) df[, diff_to_TTS := NA_real_] else df[, diff_to_TTS := as.numeric(diff_to_TTS)]
    if (!"ref_exons" %in% names(df)) df[, ref_exons := NA_real_] else df[, ref_exons := as.numeric(ref_exons)]
  })

  # ---- ORIGINAL TP/PTP CLASSIFICATION RULES ----
  df[, label := NA_character_]

  # Rule 1: SQANTI reference_match is TP
  df[subcategory == "reference_match", label := "TP"]

  # Rule 2: FSM mono-exon with both ends within 50bp is TP
  df[
    is.na(label) & structural_category == "full-splice_match" &
      !is.na(ref_exons) & ref_exons == 1 &
      !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
      abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50,
    label := "TP"
  ]

  # Rule 3: For ref_length > 3000, FSM with TSS and TTS within 100bp is TP
  df[
    is.na(label) & !is.na(ref_length) & ref_length > 3000 &
      structural_category == "full-splice_match" &
      !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
      abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100,
    label := "TP"
  ]

  # Everything else is PTP
  df[is.na(label), label := "PTP"]

  # Collapse to transcript-level label: if any TP for a transcript, label TP else PTP
  tx_lab <- df[, .(label = ifelse(any(label == "TP", na.rm = TRUE), "TP", "PTP")), by = associated_transcript_nov]

  # Get median ref_length per transcript for binning
  tx_len <- df[!is.na(ref_length), .(ref_length = median(ref_length, na.rm = TRUE)), by = associated_transcript_nov]

  df <- unique(merge(df[, .(pbid, associated_transcript_nov)], tx_lab, by = "associated_transcript_nov", all.x = TRUE))
  df <- merge(df, tx_len, by = "associated_transcript_nov", all.x = TRUE)
  df
}

# ---- 4. Load read stats ----
load_read_stats_filter <- function(file_path, keep_pbids) {
  usecols <- c("id", "pbid", "read_length")
  dt <- fread(file_path, sep = "\t", header = TRUE, select = usecols, data.table = TRUE, showProgress = FALSE)
  if (length(keep_pbids)) dt <- dt[pbid %in% keep_pbids]
  dt <- dt[!is.na(pbid) & !is.na(read_length)]
  return(dt)
}

sample_species <- c(
  wtc11_cdna_pacbio = "human",
  es_cdna_pacbio    = "mouse"
)

build_reads_dt <- function(class_files_map) {
  all_reads <- list()
  for (sample in names(class_files_map)) {
    class_paths <- as.character(class_files_map[[sample]])
    class_paths <- class_paths[file.exists(class_paths)]
    if (length(class_paths) == 0) {
      cat("[WARN] Missing classification for:", sample, "\n")
      next
    }
    species <- sample_species[[sample]]
    tusco_tx_set <- tusco_tx[[species]]
    tx_len_map <- exonic_lengths[[species]]

    class_df <- load_sqanti_classification(class_paths, tusco_tx_set)
    if (nrow(class_df) == 0) {
      cat("[WARN] No TUSCO TP/PTP in:", sample, "\n")
      next
    }

    keep_pbids <- unique(class_df$pbid)

    read_stats_path <- READ_STATS_FILES[[sample]]
    if (!file.exists(read_stats_path)) {
      cat("[WARN] Missing read stats:", read_stats_path, "\n")
      next
    }
    reads_dt <- load_read_stats_filter(read_stats_path, keep_pbids)
    if (nrow(reads_dt) == 0) {
      cat("[WARN] No reads for pbids in:", sample, "\n")
      next
    }

    merged <- merge(reads_dt, class_df, by = "pbid", all.x = TRUE, allow.cartesian = TRUE)
    # Remove mono-exon transcripts for this species
    if (length(multi_exon_tx[[species]])) {
      merged <- merged[associated_transcript_nov %in% multi_exon_tx[[species]]]
    } else {
      merged <- merged[0]
    }
    merged[, exonic_length := as.numeric(tx_len_map[associated_transcript_nov])]
    merged <- merged[!is.na(exonic_length)]
    merged[, f := as.numeric(read_length) / as.numeric(exonic_length)]
    merged[, species := species]

    all_reads[[length(all_reads) + 1L]] <- merged[, .(id, pbid, f, label, species, associated_transcript_nov, exonic_length)]
  }
  if (length(all_reads) == 0) {
    return(data.table())
  }
  rbindlist(all_reads, use.names = TRUE, fill = TRUE)
}

# Build datasets from BOTH ref_evl AND novel_evl
cat("Building read datasets from ref_evl...\n")
reads_dt_ref <- build_reads_dt(CLASS_FILES_REF)
cat("Building read datasets from novel_evl...\n")
reads_dt_novel <- build_reads_dt(CLASS_FILES_NOVEL)

# Combine ref and novel
reads_dt_all <- rbindlist(list(reads_dt_ref, reads_dt_novel), use.names = TRUE, fill = TRUE)

if (nrow(reads_dt_all) == 0) {
  message("No per-read records for either ref or novel; skipping fig-s5 plot.")
  quit(status = 0)
}

cat("Total reads after combining ref+novel:", nrow(reads_dt_all), "\n")

# ---- 5. Create "All" bin and length bins ----
# Use exonic_length for binning (computed from TUSCO GTF)
reads_all_bin <- copy(reads_dt_all)
reads_all_bin[, ref_len_bin := "All"]

# Create length bins based on exonic_length
reads_dt_all[, ref_len_bin := cut(
  exonic_length,
  breaks = c(0, 1000, 2000, 3000, 5000, Inf),
  include.lowest = TRUE,
  right = TRUE,
  labels = c("<=1kb", "1-2kb", "2-3kb", "3-5kb", ">5kb")
)]

# Combine "All" bin with length-binned data
plot_dt <- rbind(reads_all_bin, reads_dt_all[!is.na(ref_len_bin)], use.names = TRUE, fill = TRUE)

# Set factor levels (All first)
plot_dt[, ref_len_bin := factor(ref_len_bin,
  levels = c("All", "<=1kb", "1-2kb", "2-3kb", "3-5kb", ">5kb")
)]

plot_dt <- plot_dt[!is.na(ref_len_bin) & !is.na(f) & !is.na(label) & !is.na(species)]
plot_dt[, Species := ifelse(tolower(species) == "human", "Human", "Mouse")]
plot_dt[, Group := factor(label, levels = c("TP", "PTP"))]

# Exclude Mouse PTP <=1kb outlier bin (n=12 with f≈3.0, likely artifacts)
outlier_rows <- plot_dt[Species == "Mouse" & Group == "PTP" & ref_len_bin == "<=1kb", .N]
if (outlier_rows > 0) {
  cat("Excluding Mouse PTP <=1kb outlier bin (n=", outlier_rows, ", median f≈3.0)\n", sep = "")
  plot_dt <- plot_dt[!(Species == "Mouse" & Group == "PTP" & ref_len_bin == "<=1kb")]
}

# Also exclude Mouse TP <=1kb since no PTP reads remain in that bin
mouse_tp_1kb <- plot_dt[Species == "Mouse" & Group == "TP" & ref_len_bin == "<=1kb", .N]
if (mouse_tp_1kb > 0) {
  cat("Excluding Mouse TP <=1kb bin (n=", mouse_tp_1kb, ") for consistency\n", sep = "")
  plot_dt <- plot_dt[!(Species == "Mouse" & ref_len_bin == "<=1kb")]
}

cat("Plot data rows:", nrow(plot_dt), "\n")

# ---- 6. Summary statistics ----
n_by_bin <- plot_dt[, .(n = .N), by = .(Species, ref_len_bin, Group)]

# Calculate total n per Species × bin (for x-axis labels)
n_total_by_bin <- plot_dt[, .(n_total = .N), by = .(Species, ref_len_bin)]

# Create x-axis labels with sample sizes
# Format: "All\n(n=5525)"
x_labels_human <- n_total_by_bin[Species == "Human", setNames(
  paste0(ref_len_bin, "\n(n=", format(n_total, big.mark = ","), ")"),
  ref_len_bin
)]
x_labels_mouse <- n_total_by_bin[Species == "Mouse", setNames(
  paste0(ref_len_bin, "\n(n=", format(n_total, big.mark = ","), ")"),
  ref_len_bin
)]

summ_stats <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(list(n = 0L, mean = NA_real_, sd = NA_real_, median = NA_real_, q1 = NA_real_, q3 = NA_real_, iqr = NA_real_))
  }
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

plot_summary <- plot_dt[, c(summ_stats(f)), by = .(Species, Group, ref_len_bin)]
fwrite(plot_summary, OUTPUT_TSV, sep = "\t")
message("Wrote: ", OUTPUT_TSV)

# ---- 7. Plot styling ----
group_colors <- c(
  TP  = "#6BAED6",
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

# ---- 8. Create unified plot ----
dodge <- position_dodge(width = 0.8)
p <- ggplot(plot_dt, aes(x = ref_len_bin, y = f, fill = Group)) +
  geom_violin(
    aes(group = interaction(ref_len_bin, Group, drop = TRUE)),
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
  # Add sized points above violins to show sample size (n)
  geom_point(
    data = plot_summary,
    aes(x = ref_len_bin, y = 1.85, size = n, fill = Group),
    position = dodge,
    shape = 21,
    color = "black",
    stroke = 0.4,
    alpha = 0.9
  ) +
  scale_size_area(
    name = "n (reads)",
    max_size = 8,
    limits = c(0, 9000),
    breaks = c(50, 500, 2000, 5000, 9000),
    labels = scales::comma
  ) +
  facet_wrap(~Species, nrow = 1, scales = "free_x", labeller = labeller(Species = function(x) {
    # Add sample size info to facet labels
    sapply(x, function(sp) {
      total_n <- n_total_by_bin[Species == sp & ref_len_bin == "All", n_total]
      paste0(sp, " (n=", format(total_n, big.mark = ","), ")")
    })
  })) +
  scale_fill_manual(values = group_colors, name = "Group", labels = c("TP", "PTP"), drop = FALSE) +
  coord_cartesian(ylim = c(0, 2.0)) +
  labs(
    x = "Transcript exonic length",
    y = expression(italic(f) == L[read] / L[exonic]),
    fill = "Group"
  ) +
  guides(
    fill = guide_legend(
      title = "Group",
      keyheight = grid::unit(2.5, "mm"),
      keywidth = grid::unit(4, "mm")
    )
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(OUTPUT_PDF, p, width = 6.8, height = 2.5, useDingbats = FALSE)
message("Wrote: ", OUTPUT_PDF)
message("Figure S5 generated successfully with original TP/PTP rules + 'All' bin")
