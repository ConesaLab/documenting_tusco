#!/usr/bin/env Rscript

# Generate Table S1: Cosine similarity (cosim) between SIRV- and TUSCO-derived
# performance profiles across sequencing pipelines for mouse ES and human WTC11.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(rtracklayer)
})

# Silence R CMD check notes for NSE columns used by dplyr verbs
utils::globalVariables(c(
  "type", "transcript_id", "structural_category", "associated_transcript",
  "associated_gene", "isoform", "subcategory", "ref_exons", "diff_to_TSS",
  "diff_to_TTS", "ref_transcript_id", "ensembl", "refseq", "gene_name",
  "id_type"
))

# ------------------------------------------------------------
# Config: absolute paths (match figure-3 scripts)
# ------------------------------------------------------------
sirv_file        <- file.path("..","..","..","data","raw","spike-ins","lrgasp_sirvs.gtf")

# Human
human_tusco_file <- file.path("..","..","..","data","processed","tusco","hsa","tusco_human.tsv")
human_data_dir   <- file.path("..","..","..","data","raw","lrgasp","human")
human_pipelines  <- c(
  "WTC11_drna_ont",
  "WTC11_cdna_ont",
  "WTC11_cdna_pacbio",
  "WTC11_drna_ont_ls",
  "WTC11_cdna_ont_ls",
  "WTC11_cdna_pacbio_ls"
)

# Mouse
mouse_tusco_file <- file.path("..","..","..","data","processed","tusco","mmu","tusco_mouse.tsv")
mouse_data_dir   <- file.path("..","..","..","data","raw","lrgasp","mouse")
mouse_pipelines  <- c(
  "es_drna_ont",
  "es_cdna_ont",
  "es_cdna_pacbio",
  "es_drna_ont_ls",
  "es_cdna_ont_ls",
  "es_cdna_pacbio_ls"
)

plots_dir <- base::file.path("..", "tables")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
read_tsv_safe <- function(file_path, col_names = TRUE, ...) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  tryCatch({
    readr::read_tsv(file_path, col_names = col_names, ...)
  }, error = function(e) {
    stop("Error reading file ", file_path, ": ", e$message)
  })
}

format_pipeline_name <- function(prefix, cell_type = c("human", "mouse")) {
  cell_type <- match.arg(cell_type)
  parts <- strsplit(prefix, "_")[[1]]
  name_map <- c(drna = "dRNA", cdna = "cDNA", ont = "ONT", pacbio = "PacBio")
  if (cell_type == "human") {
    formatted_parts <- c("WTC11")
  } else {
    formatted_parts <- c("ES")
  }
  for (p in parts[-1]) {
    if (p %in% names(name_map)) {
      formatted_parts <- c(formatted_parts, name_map[p])
    } else if (p == "ls") {
      formatted_parts <- c(formatted_parts, "LS")
    } else {
      formatted_parts <- c(formatted_parts, toupper(p))
    }
  }
  paste(formatted_parts, collapse = " ")
}

# Compute cosine similarity using the same metric definitions as figure 3A
calculate_cosine_similarity <- function(pipeline_prefix, data_dir, sirv_file, tusco_file) {
  class_file <- file.path(data_dir, paste0(pipeline_prefix, "/", pipeline_prefix, "_classification.txt"))
  if (!file.exists(class_file)) stop("Classification file not found: ", class_file)

  classification_data <- read_tsv_safe(class_file)

  # ---- SIRVs ----
  sirv_gtf <- rtracklayer::import(sirv_file)
  sirv_gtf_df <- as.data.frame(sirv_gtf)

  annotation_data_sirv <- sirv_gtf_df %>%
    dplyr::filter(type == "exon") %>%
    dplyr::distinct(transcript_id) %>%
    dplyr::rename(ref_transcript_id = transcript_id)
  rSIRV <- nrow(annotation_data_sirv)

  classification_data_sirv <- classification_data %>% mutate(id_type = "transcript_id")

  classification_data_cleaned_sirv <- classification_data_sirv %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(
      classification_data_sirv %>%
        dplyr::filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_transcript, sep = "_")
    ) %>%
    dplyr::mutate(
      associated_gene = stringr::str_remove(associated_gene, "\\.\\d+$"),
      associated_transcript = stringr::str_remove(associated_transcript, "\\.\\d+$")
    ) %>%
    dplyr::distinct(isoform, associated_transcript, .keep_all = TRUE) %>%
    dplyr::arrange(isoform)

  sirv_chromosomes <- unique(sirv_gtf_df$seqnames)
  classification_data_cleaned_sirv <- classification_data_cleaned_sirv %>%
    dplyr::filter(chrom %in% sirv_chromosomes)

  SIRV_transcripts <- classification_data_cleaned_sirv %>% dplyr::filter(grepl("SIRV", chrom))

  SIRV_RM <- SIRV_transcripts %>%
    dplyr::mutate(mono_thresh = ifelse(!is.na(ref_length) & ref_length > 3000, 100, 50)) %>%
    dplyr::filter(
      subcategory == "reference_match" |
      (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
         !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
         abs(diff_to_TSS) <= mono_thresh & abs(diff_to_TTS) <= mono_thresh) |
      (!is.na(ref_length) & ref_length > 3000 &
       structural_category == "full-splice_match" &
       !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
       abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100)
    )

  TP_sirv <- SIRV_RM
  TP_SIRV <- unique(TP_sirv$associated_transcript)
  PTP_sirv <- SIRV_transcripts %>%
    dplyr::filter((structural_category %in% c("full-splice_match", "incomplete-splice_match") |
                    subcategory == "mono-exon_by_intron_retention") & !associated_transcript %in% TP_sirv$associated_transcript)

  detected_sirv_transcripts <- unique(classification_data_cleaned_sirv$associated_transcript)
  FN_sirv <- annotation_data_sirv %>% dplyr::filter(!(ref_transcript_id %in% detected_sirv_transcripts))

  fsm_ism_count_sirv <- SIRV_transcripts %>%
    dplyr::filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    nrow()

  non_redundant_sensitivity_sirv <- length(TP_SIRV) / rSIRV
  non_redundant_precision_sirv  <- if (nrow(SIRV_transcripts) > 0) nrow(TP_sirv) / nrow(SIRV_transcripts) else NA
  redundant_precision_sirv      <- if (fsm_ism_count_sirv > 0) (nrow(TP_sirv) + nrow(PTP_sirv)) / nrow(SIRV_transcripts) else NA
  positive_detection_rate_sirv  <- if (rSIRV > 0) length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript))) / rSIRV else NA
  false_discovery_rate_sirv     <- if (nrow(SIRV_transcripts) > 0) (nrow(SIRV_transcripts) - nrow(SIRV_RM)) / nrow(SIRV_transcripts) else NA
  unique_tp_ptp_sirv            <- length(unique(c(TP_sirv$associated_transcript, PTP_sirv$associated_transcript)))
  redundancy_sirv               <- if (unique_tp_ptp_sirv > 0) fsm_ism_count_sirv / unique_tp_ptp_sirv else NA

  # ---- TUSCO ----
  tusco_df <- readr::read_delim(
    tusco_file,
    delim = "\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = readr::cols(.default = "c"),
    trim_ws = TRUE
  )
  if (ncol(tusco_df) == 1) {
    tusco_df <- readr::read_table(
      tusco_file,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = readr::cols(.default = "c")
    )
  }
  annotation_data_tusco <- tusco_df %>% dplyr::select(ensembl, refseq, gene_name) %>% dplyr::distinct()
  rTUSCO <- nrow(annotation_data_tusco)

  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Z0-9]+$"
  )

  classification_data_cleaned_tusco <- classification_data %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(
      classification_data %>%
        dplyr::filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_gene, sep = "_")
    ) %>%
    dplyr::mutate(associated_gene = stringr::str_remove(associated_gene, "\\.\\d+$")) %>%
    dplyr::mutate(
      id_type = dplyr::case_when(
        stringr::str_detect(associated_gene, patterns$ensembl)   ~ "ensembl",
        stringr::str_detect(associated_gene, patterns$refseq)    ~ "refseq",
        stringr::str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
        TRUE ~ "unknown"
      )
    ) %>%
    dplyr::distinct(isoform, associated_gene, .keep_all = TRUE) %>%
    dplyr::arrange(isoform)

  id_summary_tusco <- classification_data_cleaned_tusco %>% dplyr::count(id_type, sort = TRUE)
  if (nrow(id_summary_tusco) == 0) stop("No id_type classifications were made for TUSCO.")
  top_id_type_tusco <- id_summary_tusco$id_type[1]

  TUSCO_transcripts <- classification_data_cleaned_tusco %>%
    dplyr::filter(associated_gene %in% annotation_data_tusco[[top_id_type_tusco]])

  TUSCO_RM <- TUSCO_transcripts %>%
    dplyr::filter(
      # TP rule (uniform): RM or FSM mono-exon with both ends within 50bp
      subcategory == "reference_match" |
      (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
         !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
         abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50) |
      # Rule 3: For ref_length > 3000, FSM with TSS and TTS within 100bp is TP
      (!is.na(ref_length) & ref_length > 3000 & 
       structural_category == "full-splice_match" &
       !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
       abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100)
    )

  TP_tusco  <- TUSCO_RM
  TP_TUSCO  <- unique(TP_tusco$associated_gene)
  PTP_tusco <- TUSCO_transcripts %>%
    dplyr::filter((structural_category %in% c("full-splice_match", "incomplete-splice_match") |
                    subcategory == "mono-exon_by_intron_retention") & !isoform %in% TP_tusco$isoform)
  FN_tusco <- annotation_data_tusco %>% dplyr::filter(!(!!rlang::sym(top_id_type_tusco) %in% TUSCO_transcripts$associated_gene))
  fsm_ism_count_tusco <- TUSCO_transcripts %>% dplyr::filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>% nrow()

  non_redundant_sensitivity_tusco <- length(TP_TUSCO) / rTUSCO
  non_redundant_precision_tusco  <- if (nrow(TUSCO_transcripts) > 0) nrow(TP_tusco) / nrow(TUSCO_transcripts) else NA
  redundant_precision_tusco      <- if (fsm_ism_count_tusco > 0) (nrow(TP_tusco) + nrow(PTP_tusco)) / nrow(TUSCO_transcripts) else NA
  positive_detection_rate_tusco  <- if (rTUSCO > 0) length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene))) / rTUSCO else NA
  false_discovery_rate_tusco     <- if (nrow(TUSCO_transcripts) > 0) (nrow(TUSCO_transcripts) - nrow(TUSCO_RM)) / nrow(TUSCO_transcripts) else NA
  unique_tp_ptp_tusco            <- length(unique(c(TP_tusco$associated_gene, PTP_tusco$associated_gene)))
  redundancy_tusco               <- if (unique_tp_ptp_tusco > 0) fsm_ism_count_tusco / unique_tp_ptp_tusco else NA

  # Metrics vectors (percentages)
  metrics_values_sirv <- c(
    non_redundant_sensitivity_sirv * 100,
    non_redundant_precision_sirv * 100,
    if (!is.na(redundancy_sirv) && redundancy_sirv != 0) (1 / redundancy_sirv) * 100 else 0,
    if (!is.na(false_discovery_rate_sirv)) (100 - (false_discovery_rate_sirv * 100)) else 0,
    if (!is.na(positive_detection_rate_sirv)) positive_detection_rate_sirv * 100 else 0,
    if (!is.na(redundant_precision_sirv)) redundant_precision_sirv * 100 else 0
  )
  metrics_values_tusco <- c(
    non_redundant_sensitivity_tusco * 100,
    non_redundant_precision_tusco  * 100,
    if (!is.na(redundancy_tusco) && redundancy_tusco != 0) (1 / redundancy_tusco) * 100 else 0,
    if (!is.na(false_discovery_rate_tusco)) (100 - (false_discovery_rate_tusco * 100)) else 0,
    if (!is.na(positive_detection_rate_tusco)) positive_detection_rate_tusco * 100 else 0,
    if (!is.na(redundant_precision_tusco)) redundant_precision_tusco * 100 else 0
  )

  # Cosine similarity
  if (any(is.na(metrics_values_sirv)) || any(is.na(metrics_values_tusco))) return(NA_real_)
  denom <- sqrt(sum(metrics_values_sirv^2)) * sqrt(sum(metrics_values_tusco^2))
  if (denom == 0) return(NA_real_)
  sum(metrics_values_sirv * metrics_values_tusco) / denom
}

# ------------------------------------------------------------
# Run for both species and assemble the table
# ------------------------------------------------------------
message("Computing cosine similarities for human pipelines...")
human_rows <- lapply(human_pipelines, function(p) {
  cosim <- tryCatch({
    calculate_cosine_similarity(p, human_data_dir, sirv_file, human_tusco_file)
  }, error = function(e) {
    message("  ", p, ": ", e$message)
    NA_real_
  })
  data.frame(
    `Cell type` = "Human WTC11",
    `Pipeline`  = format_pipeline_name(p, "human"),
    `cosim`     = cosim,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

message("Computing cosine similarities for mouse pipelines...")
mouse_rows <- lapply(mouse_pipelines, function(p) {
  cosim <- tryCatch({
    calculate_cosine_similarity(p, mouse_data_dir, sirv_file, mouse_tusco_file)
  }, error = function(e) {
    message("  ", p, ": ", e$message)
    NA_real_
  })
  data.frame(
    `Cell type` = "Mouse ES",
    `Pipeline`  = format_pipeline_name(p, "mouse"),
    `cosim`     = cosim,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

table_df <- bind_rows(mouse_rows, human_rows) %>%
  mutate(cosim = round(as.numeric(cosim), 4))

csv_path <- file.path(plots_dir, "table_s1.csv")
readr::write_csv(table_df, csv_path)

message("Table S1 written (CSV): ", csv_path)
