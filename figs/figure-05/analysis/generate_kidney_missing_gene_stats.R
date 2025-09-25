suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

project_root <- "/Users/tianyuan/Desktop/github_dev/tusco-paper"
isoseq_root <- file.path(project_root, "data/raw/nih/single_sample")
plot_dir <- file.path(project_root, "figs/figure-05/plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
out_tsv <- file.path(plot_dir, "fig-5d_kidney_missing_gene_stats.tsv")

tusco_kidney_tsv <- file.path(project_root, "data/processed/tusco/mmu/tusco_mouse_kidney.tsv")
tp_gtf_path <- file.path(project_root, "data/processed/tusco/mmu/tusco_mouse_kidney.gtf")

message("Loading TUSCO kidney reference…")
load_tusco_reference_simple <- function(tusco_path) {
  tusco_lines <- readLines(tusco_path)
  tusco_data_lines <- tusco_lines[!grepl("^#", tusco_lines)]
  tusco_temp <- tempfile()
  writeLines(tusco_data_lines, tusco_temp)
  tusco_df <- readr::read_delim(
    tusco_temp, delim = "\t",
    col_names = c("ensembl","transcript","gene_name","gene_id_num","refseq","prot_refseq"),
    col_types = readr::cols(.default = "c"), trim_ws = TRUE, show_col_types = FALSE
  )
  unlink(tusco_temp)
  tusco_df %>%
    mutate(
      ensembl = stringr::str_remove(ensembl, "\\.\\d+$"),
      transcript = stringr::str_remove(transcript, "\\.\\d+$"),
      refseq = stringr::str_remove(refseq, "\\.\\d+$")
    ) %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()
}

tusco_ref <- load_tusco_reference_simple(tusco_kidney_tsv)
refseq_to_ensembl <- setNames(
  tusco_ref$ensembl[!is.na(tusco_ref$refseq) & tusco_ref$refseq != ""],
  tusco_ref$refseq[!is.na(tusco_ref$refseq) & tusco_ref$refseq != ""]
)
gene_to_ensembl <- setNames(
  tusco_ref$ensembl[!is.na(tusco_ref$gene_name) & tusco_ref$gene_name != ""],
  tusco_ref$gene_name[!is.na(tusco_ref$gene_name) & tusco_ref$gene_name != ""]
)

infer_tissue <- function(sample_basename) {
  if (startsWith(sample_basename, "B")) return("Brain")
  if (startsWith(sample_basename, "K")) return("Kidney")
  return("Unknown")
}

get_tp_genes_for_sample <- function(sample_dir) {
  cl <- file.path(sample_dir, paste0(basename(sample_dir), "_classification.txt"))
  if (!file.exists(cl)) return(character())
  df <- readr::read_tsv(cl, show_col_types = FALSE, progress = FALSE)
  df %>%
    dplyr::filter(
      subcategory == "reference_match" |
      (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
         !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
         abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50) |
      # Rule 3: For ref_length > 3000, FSM with TSS and TTS within 100bp is TP
      (!is.na(ref_length) & ref_length > 3000 & 
       structural_category == "full-splice_match" &
       !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
       abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100)
    ) %>%
    mutate(id_norm = stringr::str_remove(associated_gene, "\\.\\d+$")) %>%
    mutate(gene = dplyr::case_when(
      grepl("^(ENSMUSG|ENSG)\\d{11}$", id_norm) ~ id_norm,
      !is.na(refseq_to_ensembl[id_norm]) ~ unname(refseq_to_ensembl[id_norm]),
      !is.na(gene_to_ensembl[id_norm]) ~ unname(gene_to_ensembl[id_norm]),
      TRUE ~ NA_character_
    )) %>%
    dplyr::filter(!is.na(gene) & gene != "") %>%
    distinct(gene) %>%
    pull(gene)
}

message("Scanning kidney Iso-Seq AR samples for TP genes…")
sample_dirs <- list.dirs(isoseq_root, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[grepl("/K[0-9]+\\.isoforms$", sample_dirs)]
tp_sets <- lapply(sample_dirs, get_tp_genes_for_sample)
any_tp_genes <- unique(unlist(tp_sets, use.names = FALSE))
if (length(any_tp_genes) == 0) {
  warning("No TP genes found in any kidney sample; all TUSCO kidney genes will be considered missing.")
}

message("Reading TUSCO kidney GTF and computing FN transcript exonic lengths…")
gtf <- readr::read_tsv(
  tp_gtf_path,
  comment = "#",
  col_names = c("chrom","source","feature","start","end","score","strand","frame","attribute"),
  col_types = cols(
    chrom = col_character(), source = col_character(), feature = col_character(),
    start = col_integer(), end = col_integer(), score = col_character(),
    strand = col_character(), frame = col_character(), attribute = col_character()
  )
)

extract_attr <- function(attr, key) {
  m <- stringr::str_match(attr, paste0(key, " \"([^\"]+)\""))
  out <- m[,2]
  out <- stringr::str_remove(out, "\\.\\d+$")
  out
}

gtf <- gtf %>% mutate(
  gene_id = extract_attr(attribute, "gene_id"),
  transcript_id = extract_attr(attribute, "transcript_id")
)

fn_exons <- gtf %>%
  filter(feature == "exon", !is.na(gene_id), !is.na(transcript_id)) %>%
  filter(!(gene_id %in% any_tp_genes))

fn_tx_stats <- fn_exons %>%
  group_by(gene_id, transcript_id) %>%
  summarise(
    exonic_bases = sum(end - start + 1L),
    .groups = "drop"
  ) %>%
  arrange(gene_id, transcript_id)

message("Writing ", nrow(fn_tx_stats), " FN transcript rows to ", out_tsv)
readr::write_tsv(fn_tx_stats, out_tsv)
message("Done.")
