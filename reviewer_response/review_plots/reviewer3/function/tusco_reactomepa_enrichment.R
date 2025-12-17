#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && any(args %in% c("-h", "--help"))) {
  cat(
    "Usage: Rscript tusco_reactomepa_enrichment.R [--outdir DIR] [--human TSV] [--mouse TSV] [--no-install]\n",
    "\nRuns Reactome pathway ORA using ReactomePA (Entrez IDs).\n",
    "\nDefaults:\n",
    "  --outdir reviewer_response/review_plots/function\n",
    "  --human  data/processed/tusco/tusco_human.tsv\n",
    "  --mouse  data/processed/tusco/tusco_mouse.tsv\n",
    sep = ""
  )
  quit(status = 0)
}

parse_flag_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) {
    return(default)
  }
  if (idx == length(args)) {
    stop(sprintf("Missing value for %s", flag))
  }
  args[[idx + 1]]
}

outdir <- parse_flag_value("--outdir", "reviewer_response/review_plots/function")
human_path <- parse_flag_value("--human", "data/processed/tusco/tusco_human.tsv")
mouse_path <- parse_flag_value("--mouse", "data/processed/tusco/tusco_mouse.tsv")
no_install <- any(args %in% "--no-install")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

need_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  if (no_install) {
    stop(sprintf("Missing package '%s' (rerun without --no-install to install).", pkg))
  }
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

need_pkg("readr")
need_pkg("dplyr")
need_pkg("tibble")
need_pkg("ggplot2")
need_pkg("ReactomePA", bioc = TRUE)
need_pkg("enrichplot", bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ReactomePA)
  library(enrichplot)
})

read_tusco_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  df <- readr::read_tsv(
    I(lines),
    col_names = c("ensembl_gene", "ensembl_transcript", "symbol", "entrez_id", "refseq_mrna", "refseq_protein"),
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols(
      ensembl_gene = col_character(),
      ensembl_transcript = col_character(),
      symbol = col_character(),
      entrez_id = col_character(),
      refseq_mrna = col_character(),
      refseq_protein = col_character()
    )
  )
  df <- df %>%
    mutate(entrez_id = trimws(entrez_id)) %>%
    filter(!is.na(entrez_id), entrez_id != "")
  df
}

save_dotplot <- function(enrich_obj, filename, title) {
  out_path <- file.path(outdir, filename)
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    p <- ggplot() +
      theme_void() +
      ggtitle(title) +
      annotate("text", x = 0, y = 0, label = "No significant terms at chosen cutoff", size = 4)
    ggsave(out_path, p, width = 8.5, height = 4.5, units = "in", dpi = 300)
    return(invisible(NULL))
  }
  p <- enrichplot::dotplot(enrich_obj, showCategory = 20) +
    ggtitle(title) +
    theme(plot.title = element_text(size = 12))
  ggsave(out_path, p, width = 10, height = 7, units = "in", dpi = 300)
  invisible(NULL)
}

write_results <- function(enrich_obj, path) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    readr::write_tsv(tibble(), path)
    return(invisible(NULL))
  }
  readr::write_tsv(as_tibble(enrich_obj@result), path)
  invisible(NULL)
}

filter_enrich <- function(enrich_obj, padj_cutoff = 0.05) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(enrich_obj)
  }
  res <- as_tibble(enrich_obj@result) %>% filter(!is.na(p.adjust), p.adjust <= padj_cutoff)
  enrich_obj@result <- as.data.frame(res)
  enrich_obj
}

run_species <- function(label, tusco_path, organism_code) {
  tusco <- read_tusco_tsv(tusco_path)
  genes <- unique(tusco$entrez_id)
  message(sprintf("[%s] genes=%d", label, length(genes)))

  enr_all <- suppressMessages(ReactomePA::enrichPathway(
    gene = genes,
    organism = organism_code,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    readable = TRUE
  ))
  enr_sig <- filter_enrich(enr_all, 0.05)

  write_results(enr_all, file.path(outdir, sprintf("%s_reactomepa_all.tsv", label)))
  write_results(enr_sig, file.path(outdir, sprintf("%s_reactomepa.tsv", label)))
  save_dotplot(enr_sig, sprintf("%s_REACTOMEPA_dotplot.pdf", label), sprintf("%s TUSCO: ReactomePA enrichment (FDR <= 0.05)", toupper(label)))
  invisible(NULL)
}

run_species("human", human_path, "human")
run_species("mouse", mouse_path, "mouse")
message(sprintf("Done. Outputs written to %s", outdir))

