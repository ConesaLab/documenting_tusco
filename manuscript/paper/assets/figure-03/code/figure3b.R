###############################################################################
# Figure 3b (human + mouse) in one script
#
# Usage (from this folder):
#   Rscript figure3b.R                 # runs both
#   Rscript figure3b.R --species human # runs only human
#   Rscript figure3b.R --species mouse # runs only mouse
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

colors_sirvs <- "#cab2d6"
colors_tusco_human <- "#a8d5a0"
colors_tusco_mouse <- "#1b9e77"

resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
  }
  candidates[[1]]
}

read_tsv_safe <- function(file_path, col_names = TRUE, ...) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  tryCatch(
    read_tsv(file_path, col_names = col_names, show_col_types = FALSE, progress = FALSE, ...),
    error = function(e) stop("Error reading file ", file_path, ": ", e$message)
  )
}

# Lightweight GTF importer: uses rtracklayer if available, otherwise parses attributes
import_gtf_df <- function(gtf_path) {
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    return(as.data.frame(rtracklayer::import(gtf_path)))
  }
  cols <- readr::cols(
    seqnames = readr::col_character(),
    source = readr::col_character(),
    type = readr::col_character(),
    start = readr::col_integer(),
    end = readr::col_integer(),
    score = readr::col_character(),
    strand = readr::col_character(),
    frame = readr::col_character(),
    attribute = readr::col_character()
  )
  df <- readr::read_tsv(
    gtf_path,
    comment = "#",
    col_names = c("seqnames", "source", "type", "start", "end", "score", "strand", "frame", "attribute"),
    col_types = cols,
    progress = FALSE
  )
  extract_attr <- function(attr, key) {
    m <- regmatches(attr, regexpr(paste0(key, " \"[^\"]+\""), attr))
    sub(paste0(key, " \"([^\"]+)\""), "\\1", m)
  }
  df$transcript_id <- NA_character_
  df$gene_id <- NA_character_
  has_tid <- grepl("transcript_id \"", df$attribute, fixed = TRUE)
  df$transcript_id[has_tid] <- extract_attr(df$attribute[has_tid], "transcript_id")
  has_gid <- grepl("gene_id \"", df$attribute, fixed = TRUE)
  df$gene_id[has_gid] <- extract_attr(df$attribute[has_gid], "gene_id")
  df
}

compute_sirv_bigcats <- function(classification_data, transcript_gtf_file, sirv_gtf_file) {
  sirv_gtf_df <- import_gtf_df(sirv_gtf_file)

  if (!"transcript_id" %in% names(sirv_gtf_df)) {
    if ("Parent" %in% names(sirv_gtf_df)) {
      sirv_gtf_df$transcript_id <- sirv_gtf_df$Parent
    } else {
      stop("No 'transcript_id' or 'Parent' attribute in SIRV GTF. Cannot continue.")
    }
  }

  annotation_data_sirv <- sirv_gtf_df %>%
    dplyr::filter(type == "exon") %>%
    dplyr::distinct(transcript_id) %>%
    dplyr::rename(ref_transcript_id = transcript_id)

  class_sirv <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category == "fusion") %>%
        separate_rows(associated_transcript, sep = "_")
    ) %>%
    mutate(
      associated_gene = str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    distinct(isoform, associated_transcript, .keep_all = TRUE)

  sirv_chroms <- unique(sirv_gtf_df$seqnames)
  class_sirv <- class_sirv %>% filter(chrom %in% sirv_chroms)
  SIRV_transcripts <- class_sirv %>% filter(grepl("SIRV", chrom))

  TP_sirv <- SIRV_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
          !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
          abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50) |
        (!is.na(ref_length) & ref_length > 3000 &
          structural_category == "full-splice_match" &
          !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
          abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100)
    )

  PTP_sirv <- SIRV_transcripts %>%
    filter(
      (structural_category %in% c("full-splice_match", "incomplete-splice_match") |
        subcategory == "mono-exon_by_intron_retention"),
      !(associated_transcript %in% TP_sirv$associated_transcript)
    )

  ref_found <- union(TP_sirv$associated_transcript, PTP_sirv$associated_transcript)
  FN_sirv <- annotation_data_sirv %>% filter(!(ref_transcript_id %in% ref_found))

  FP_sirv <- SIRV_transcripts %>%
    filter(structural_category %in% c(
      "novel_in_catalog", "novel_not_in_catalog", "genic",
      "fusion", "antisense", "intergenic", "genic_intron"
    ) &
      subcategory != "mono-exon_by_intron_retention")

  summary_sirv <- data.frame(
    big_category = c("TP", "PTP", "FP", "FN"),
    count = c(nrow(TP_sirv), nrow(PTP_sirv), nrow(FP_sirv), nrow(FN_sirv))
  )
  total_sirv <- sum(summary_sirv$count)
  summary_sirv$percentage <- 100 * summary_sirv$count / if (total_sirv == 0) 1 else total_sirv
  summary_sirv
}

compute_tusco_bigcats <- function(classification_data, tusco_tsv_file, pipeline_prefix = NULL) {
  tusco_df <- read_delim(
    tusco_tsv_file,
    delim = "\t",
    comment = "#",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = cols(.default = "c"),
    trim_ws = TRUE,
    progress = FALSE
  )
  if (ncol(tusco_df) == 1) {
    tusco_df <- read_table2(
      tusco_tsv_file,
      comment = "#",
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = cols(.default = "c"),
      progress = FALSE
    )
  }

  if (!all(c("ensembl", "refseq", "gene_name") %in% colnames(tusco_df))) {
    stop("Tusco TSV must contain columns: ensembl, refseq, gene_name")
  }

  annotation_data_tusco <- tusco_df %>% select(ensembl, refseq, gene_name) %>% distinct()

  patterns <- list(
    ensembl = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Z0-9]+$"
  )

  class_tusco <- classification_data %>%
    filter(structural_category != "fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category == "fusion") %>%
        separate_rows(associated_gene, sep = "_")
    ) %>%
    mutate(associated_gene = str_remove(associated_gene, "\\.\\d+$")) %>%
    mutate(
      id_type = case_when(
        str_detect(associated_gene, patterns$ensembl) ~ "ensembl",
        str_detect(associated_gene, patterns$refseq) ~ "refseq",
        str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
        TRUE ~ "unknown"
      )
    ) %>%
    distinct(isoform, associated_gene, .keep_all = TRUE)

  id_summary_tusco <- class_tusco %>% count(id_type, sort = TRUE)
  known_id_summary <- id_summary_tusco %>% filter(id_type != "unknown")
  if (nrow(known_id_summary) == 0) {
    warning("No known gene ID type found in classifications for TUSCO; returning empty summary.")
    return(data.frame(big_category = c("TP", "PTP", "FP", "FN"), count = c(0, 0, 0, 0), percentage = c(0, 0, 0, 0)))
  }
  top_type <- known_id_summary$id_type[1]
  if (!is.null(pipeline_prefix)) cat("Using TUSCO ID type:", top_type, "for pipeline", pipeline_prefix, "\n")

  TUSCO_transcripts <- class_tusco %>%
    filter(id_type == top_type) %>%
    filter(associated_gene %in% annotation_data_tusco[[top_type]])

  TP_tusco <- TUSCO_transcripts %>%
    filter(
      subcategory == "reference_match" |
        (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
          !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
          abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50) |
        (!is.na(ref_length) & ref_length > 3000 &
          structural_category == "full-splice_match" &
          !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
          abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100)
    )

  PTP_tusco <- TUSCO_transcripts %>%
    filter(
      (structural_category %in% c("full-splice_match", "incomplete-splice_match") |
        subcategory == "mono-exon_by_intron_retention"),
      !(isoform %in% TP_tusco$isoform)
    )

  found_ids_tusco <- union(TP_tusco$associated_gene, PTP_tusco$associated_gene)
  FN_tusco_count <- annotation_data_tusco %>%
    filter(!(!!sym(top_type) %in% found_ids_tusco)) %>%
    nrow()

  FP_tusco <- TUSCO_transcripts %>%
    filter(structural_category %in% c(
      "novel_in_catalog", "novel_not_in_catalog", "genic",
      "fusion", "antisense", "intergenic", "genic_intron"
    ) &
      subcategory != "mono-exon_by_intron_retention")

  summary_tusco <- data.frame(
    big_category = c("TP", "PTP", "FP", "FN"),
    count = c(nrow(TP_tusco), nrow(PTP_tusco), nrow(FP_tusco), FN_tusco_count)
  )
  total_tusco <- sum(summary_tusco$count)
  summary_tusco$percentage <- 100 * summary_tusco$count / if (total_tusco == 0) 1 else total_tusco
  summary_tusco
}

process_pipeline <- function(pipeline_prefix, data_dir, sirv_gtf, tusco_tsv) {
  class_file <- file.path(data_dir, pipeline_prefix, paste0(pipeline_prefix, "_classification.txt"))
  gtf_file <- file.path(data_dir, pipeline_prefix, paste0(pipeline_prefix, "_corrected.gtf"))

  classification_data <- read_tsv_safe(class_file)

  sirv_res <- compute_sirv_bigcats(classification_data, gtf_file, sirv_gtf) %>%
    mutate(Type = "SIRVs", pipeline = pipeline_prefix)

  tusco_res <- compute_tusco_bigcats(classification_data, tusco_tsv, pipeline_prefix = pipeline_prefix) %>%
    mutate(Type = "TUSCO", pipeline = pipeline_prefix)

  bind_rows(sirv_res, tusco_res)
}

p_stars <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 1e-4) return("****")
  if (x < 1e-3) return("***")
  if (x < 1e-2) return("**")
  if (x < 0.05) return("*")
  "ns"
}

mean_sd <- function(x, mult = 1) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  data.frame(y = m, ymin = m - mult * s, ymax = m + mult * s)
}

compute_paired_stats <- function(all_data, categories, method = c("perm", "wilcox", "t")) {
  method <- match.arg(method)
  stats <- lapply(categories, function(cat_name) {
    df_sub <- all_data %>% filter(big_category == cat_name)
    wide <- df_sub %>%
      select(pipeline, Type, percentage) %>%
      pivot_wider(names_from = Type, values_from = percentage)

    if (!all(c("TUSCO", "SIRVs") %in% colnames(wide))) {
      return(data.frame(
        big_category = cat_name,
        p_value = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        mean_diff_tusco_minus_sirvs = NA_real_,
        direction = NA_character_,
        t_statistic = NA_real_,
        df = NA_real_
      ))
    }

    wide2 <- wide %>% filter(!is.na(TUSCO) & !is.na(SIRVs))
    if (nrow(wide2) < 2) {
      return(data.frame(
        big_category = cat_name,
        p_value = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        mean_diff_tusco_minus_sirvs = NA_real_,
        direction = NA_character_,
        t_statistic = NA_real_,
        df = NA_real_
      ))
    }

    diffs <- wide2$TUSCO - wide2$SIRVs
    mean_diff <- mean(diffs)

    if (method == "t") {
      test_out <- t.test(wide2$TUSCO, wide2$SIRVs, paired = TRUE, alternative = "two.sided")
      effect <- mean_diff
      ci_lower <- unname(test_out$conf.int[1])
      ci_upper <- unname(test_out$conf.int[2])
      stat <- unname(test_out$statistic)
      df <- unname(test_out$parameter)
      p_value <- unname(test_out$p.value)
    } else if (method == "wilcox") {
      # Paired Wilcoxon signed-rank test: robust to outliers / non-normality (small n)
      test_out <- suppressWarnings(
        wilcox.test(
          wide2$TUSCO,
          wide2$SIRVs,
          paired = TRUE,
          alternative = "two.sided",
          conf.int = TRUE,
          exact = FALSE
        )
      )
      # Wilcoxon estimate is a (pseudo)median shift; keep mean_diff separately for direction/plot intuition.
      effect <- mean_diff
      ci_lower <- if (!is.null(test_out$conf.int)) unname(test_out$conf.int[1]) else NA_real_
      ci_upper <- if (!is.null(test_out$conf.int)) unname(test_out$conf.int[2]) else NA_real_
      stat <- unname(test_out$statistic)
      df <- NA_real_
      p_value <- unname(test_out$p.value)
    } else {
      # Exact sign-flip permutation test on paired differences (two-sided).
      # This is robust and exact for small n (here n=6 => 64 sign flips).
      n <- length(diffs)
      if (n == 0) {
        p_value <- NA_real_
      } else {
        obs <- mean_diff
        # enumerate all sign flips
        flips <- expand.grid(rep(list(c(-1, 1)), n), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        flip_mat <- as.matrix(flips)
        perm_means <- as.vector(flip_mat %*% diffs) / n
        p_value <- (sum(abs(perm_means) >= abs(obs)) + 1) / (length(perm_means) + 1)
      }
      effect <- mean_diff
      ci_lower <- NA_real_
      ci_upper <- NA_real_
      stat <- NA_real_
      df <- NA_real_
    }

    direction <- if (mean_diff > 0) "TUSCO>SIRVs" else if (mean_diff < 0) "SIRVs>TUSCO" else "equal"

    data.frame(
      big_category = cat_name,
      p_value = p_value,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      mean_diff_tusco_minus_sirvs = mean_diff,
      direction = direction,
      t_statistic = stat,
      df = df
    )
  })
  bind_rows(stats)
}

format_p <- function(p) {
  out <- rep("p=NA", length(p))
  ok <- !is.na(p)
  out[ok] <- paste0("p=", formatC(p[ok], format = "g", digits = 2))
  out
}

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

boxed_grob <- function(grob, padding_mm = 1.5, line_width = 0.4) {
  pad <- grid::unit(padding_mm, "mm")
  
  if (inherits(grob, "gtable")) {
    w <- sum(grob$widths) + 2 * pad
    h <- sum(grob$heights) + 2 * pad
  } else {
    w <- grid::unit(1, "npc") - 2 * pad
    h <- grid::unit(1, "npc") - 2 * pad
  }

  grid::grobTree(
    grid::rectGrob(
      width = w,
      height = h,
      gp = grid::gpar(fill = NA, col = "black", lwd = line_width)
    ),
    grob
  )
}

draw_three_panel <- function(p_left, p_mid, legend_grob, widths = c(1, 1, 0.35)) {
  g1 <- ggplotGrob(p_left)
  g2 <- ggplotGrob(p_mid)
  g3 <- legend_grob
  g <- gtable::gtable_matrix(
    name = "three_panel",
    grobs = matrix(list(g1, g2, g3), nrow = 1),
    widths = grid::unit(widths, "null"),
    heights = grid::unit(1, "null")
  )
  grid::grid.draw(g)
}

run_figure3b <- function(species = c("human", "mouse"), write_files = TRUE, test_method = c("perm", "wilcox", "t")) {
  species <- match.arg(species)
  test_method <- match.arg(test_method)

  sirv_gtf_file <- resolve_path(c(
    file.path("..", "..", "..", "data", "raw", "spike-ins", "lrgasp_sirvs.gtf"),
    file.path("..", "data", "spike-ins", "lrgasp_sirvs.gtf")
  ))

  if (species == "human") {
    pipelines <- c(
      "WTC11_drna_ont",
      "WTC11_cdna_ont",
      "WTC11_cdna_pacbio",
      "WTC11_drna_ont_ls",
      "WTC11_cdna_ont_ls",
      "WTC11_cdna_pacbio_ls"
    )
    tusco_tsv_file <- resolve_path(c(
      file.path("..", "..", "..", "data", "processed", "tusco", "hsa", "tusco_human.tsv"),
      file.path("..", "data", "tusco", "tusco_human.tsv")
    ))
    data_dir <- resolve_path(c(
      file.path("..", "..", "..", "data", "raw", "lrgasp", "human"),
      file.path("..", "data", "lrgasp", "human")
    ), is_dir = TRUE)
    panel_id <- "3b-human"
    panel_title <- "Human"
    tusco_color <- colors_tusco_human
  } else {
    pipelines <- c(
      "ES_drna_ont",
      "ES_cdna_ont",
      "ES_cdna_pacbio",
      "ES_drna_ont_ls",
      "ES_cdna_ont_ls",
      "ES_cdna_pacbio_ls"
    )
    tusco_tsv_file <- resolve_path(c(
      file.path("..", "..", "..", "data", "processed", "tusco", "mmu", "tusco_mouse.tsv"),
      file.path("..", "data", "tusco", "tusco_mouse.tsv")
    ))
    data_dir <- resolve_path(c(
      file.path("..", "..", "..", "data", "raw", "lrgasp", "mouse"),
      file.path("..", "data", "lrgasp", "mouse")
    ), is_dir = TRUE)
    panel_id <- "3b-mouse"
    panel_title <- "Mouse"
    tusco_color <- colors_tusco_mouse
  }
  type_colors <- c("SIRVs" = colors_sirvs, "TUSCO" = tusco_color)

  cat("\n--- Figure 3b:", species, "---\n")
  cat("\n--- Gathering data for all pipelines ---\n")
  all_data <- lapply(pipelines, process_pipeline,
    data_dir = data_dir,
    sirv_gtf = sirv_gtf_file,
    tusco_tsv = tusco_tsv_file
  ) %>%
    bind_rows()

  categories <- c("TP", "PTP", "FP", "FN")

  cat("\n--- Paired test (two-sided):", test_method, "TUSCO vs SIRVs ---\n")
  my_signifs <- compute_paired_stats(all_data, categories, method = test_method) %>%
    mutate(star_label = vapply(p_value, p_stars, character(1)))
  print(my_signifs)

  mean_data <- all_data %>%
    group_by(big_category, Type) %>%
    summarise(
      mean_perc = mean(percentage, na.rm = TRUE),
      sd_perc = sd(percentage, na.rm = TRUE),
      .groups = "drop"
    )

  mean_data$big_category <- factor(mean_data$big_category, levels = categories)
  mean_data$Type <- factor(mean_data$Type, levels = c("TUSCO", "SIRVs"))

  bracket_data <- my_signifs %>%
    mutate(
      big_category = factor(big_category, levels = categories),
      y_position = 110,
      y_line = 105,
      xmin = as.numeric(big_category) - 0.2,
      xmax = as.numeric(big_category) + 0.2,
      dir_short = dplyr::case_when(
        direction == "TUSCO>SIRVs" ~ "T>S",
        direction == "SIRVs>TUSCO" ~ "S>T",
        TRUE ~ ""
      ),
      star_label = ifelse(is.na(star_label), "NA", star_label),
      plot_label = star_label
    )

  all_data$big_category <- factor(all_data$big_category, levels = categories)

  p_single <- ggplot(all_data, aes(x = big_category, y = percentage, fill = Type)) +
    stat_summary(
      fun = mean,
      geom = "bar",
      width = 0.42,
      color = "black",
      linewidth = 0.2,
      position = position_dodge(width = 0.65)
    ) +
    stat_summary(
      fun.data = mean_sd,
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.16,
      linewidth = 0.25,
      position = position_dodge(width = 0.65)
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.65, seed = 42),
      aes(fill = Type),
      size = 1.2,
      alpha = 0.85,
      shape = 21,
      color = "black",
      stroke = 0.25
    ) +
    geom_text(
      data = bracket_data,
      aes(x = big_category, y = y_position, label = plot_label),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = type_colors,
      breaks = c("TUSCO", "SIRVs"),
      labels = c("TUSCO" = "TUSCO evaluation", "SIRVs" = "SIRVs")
    ) +
    coord_cartesian(ylim = c(0, 110), clip = "off") +
    labs(x = NULL, y = "Percentage") +
    theme_classic(base_family = "Helvetica", base_size = 7) +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 7, face = "bold"),
      legend.position = "none",
      legend.title = element_blank(),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      plot.margin = margin(4.5, 3.5, 3.5, 3.5, unit = "mm")
    )

  plot_dir <- base::file.path("..", "plots")
  tsv_dir <- base::file.path("..", "tables")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  if (!dir.exists(tsv_dir)) dir.create(tsv_dir, recursive = TRUE)

  pdf_path <- file.path(plot_dir, paste0("figure3b-", species, ".pdf"))
  tsv_path <- file.path(tsv_dir, paste0("figure3b-", species, ".tsv"))
  if (isTRUE(write_files)) {
    pdf(file = pdf_path, width = 2.5, height = 1.8)
    print(p_single)
    dev.off()
  }

  group_counts <- all_data %>% count(big_category, Type, name = "n")
  raw_block <- all_data %>%
    left_join(group_counts, by = c("big_category", "Type")) %>%
    mutate(record_type = "raw")

  summary_block <- mean_data %>%
    mutate(record_type = "summary")

  stat_block <- my_signifs %>%
    mutate(record_type = "stat", test_method = test_method, p_value_stars = star_label) %>%
    select(big_category, test_method, p_value, ci_lower, ci_upper, mean_diff_tusco_minus_sirvs, direction, t_statistic, df, p_value_stars, record_type)

  tsv_out <- bind_rows(
    raw_block %>% mutate(figure_id = "fig-3", panel_id = panel_id),
    summary_block %>% mutate(figure_id = "fig-3", panel_id = panel_id),
    stat_block %>% mutate(figure_id = "fig-3", panel_id = panel_id)
  )
  if (isTRUE(write_files)) {
    readr::write_tsv(tsv_out, tsv_path)
  }

  extraction_file <- file.path("..", "..", "..", "manuscript", "extraction_values.tsv")
  ptp_pval <- my_signifs %>% filter(big_category == "PTP") %>% pull(p_value)
  fn_pval <- my_signifs %>% filter(big_category == "FN") %>% pull(p_value)
  if (length(ptp_pval) == 0) ptp_pval <- NA_real_
  if (length(fn_pval) == 0) fn_pval <- NA_real_

  extraction_data <- data.frame(
    figure_id = "fig-3b",
    metric_name = c(paste0(species, "_PTP_pvalue"), paste0(species, "_FN_pvalue")),
    value = c(formatC(ptp_pval, format = "e", digits = 3), formatC(fn_pval, format = "e", digits = 3)),
    ci_lower = NA,
    ci_upper = NA,
    notes = "Two-sided paired test (Wilcoxon signed-rank; TUSCO vs SIRVs); direction from estimated shift (TUSCO-SIRVs)",
    stringsAsFactors = FALSE
  )

  if (isTRUE(write_files)) {
    if (!file.exists(extraction_file)) {
      write.table(extraction_data, extraction_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("Created extraction file:", extraction_file, "\n")
    } else {
      write.table(extraction_data, extraction_file, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)
    }

    cat("Plot saved successfully:\n")
    cat(" - PDF:", pdf_path, "\n")
    cat(" - TSV:", tsv_path, "\n")
  }

  list(
    species = species,
    panel_id = panel_id,
    plot = p_single,
    all_data = all_data,
    stats = my_signifs,
    tusco_color = tusco_color,
    pdf_path = pdf_path,
    tsv_path = tsv_path
  )
}

write_figure3b_combined <- function(human_plot, mouse_plot) {
  if (!requireNamespace("gtable", quietly = TRUE)) stop("Package 'gtable' is required (comes with ggplot2).")
  if (!requireNamespace("grid", quietly = TRUE)) stop("Package 'grid' is required.")

  legend_df <- data.frame(Type = c("SIRVs", "TUSCO Human", "TUSCO Mouse"))
  p_leg <- ggplot(legend_df, aes(x = 1, y = Type, fill = Type)) +
    geom_point(shape = 21, size = 3, color = "black", stroke = 0.25) +
    scale_fill_manual(
      values = c("SIRVs" = colors_sirvs, "TUSCO Human" = colors_tusco_human, "TUSCO Mouse" = colors_tusco_mouse),
      breaks = c("SIRVs", "TUSCO Human", "TUSCO Mouse"),
      labels = c("SIRVs" = "SIRVs", "TUSCO Human" = "TUSCO (Human)", "TUSCO Mouse" = "TUSCO (Mouse)")
    ) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 3, color = "black", stroke = 0.25))) +
    theme_void(base_family = "Helvetica", base_size = 7) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 7)
    )
  leg <- extract_legend(p_leg)
  if (is.null(leg)) stop("Failed to extract legend grob.")
  leg_boxed <- boxed_grob(leg, padding_mm = 1, line_width = 0.4)

  out_pdf <- file.path("..", "plots", "figure3b.pdf")
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(file = out_pdf, width = 180 / 25.4, height = 50 / 25.4)
  grid::grid.newpage()
  draw_three_panel(human_plot, mouse_plot, leg_boxed, widths = c(1, 1, 0.4))
  dev.off()

  cat("Combined plot saved successfully:\n")
  cat(" - PDF:", out_pdf, "\n")
}

parse_args <- function(args) {
  out <- list(species = "both")
  if (length(args) == 0) return(out)
  for (i in seq_along(args)) {
    if (args[[i]] %in% c("--species", "-s") && i < length(args)) {
      out$species <- args[[i + 1]]
    }
  }
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$species) || args$species %in% c("both", "all")) {
    # Build plots for combined PDF (human left, mouse right, boxed legend to the right)
    human_res <- run_figure3b("human", write_files = TRUE, test_method = "perm")
    mouse_res <- run_figure3b("mouse", write_files = TRUE, test_method = "perm")
    write_figure3b_combined(human_res$plot, mouse_res$plot)
    return(invisible(NULL))
  }
  if (!args$species %in% c("human", "mouse")) {
    stop("Invalid --species. Use: human, mouse, or omit to run both.")
  }
  run_figure3b(args$species, write_files = TRUE, test_method = "perm")
  invisible(NULL)
}

if (sys.nframe() == 0) {
  main()
}
