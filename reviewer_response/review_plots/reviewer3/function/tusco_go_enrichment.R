#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && any(args %in% c("-h", "--help"))) {
  cat(
    "Usage: Rscript tusco_go_enrichment.R [--outdir DIR] [--human TSV] [--mouse TSV] [--no-install]\n",
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
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (bioc) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

need_pkg("readr")
need_pkg("dplyr")
need_pkg("tibble")
need_pkg("stringr")
need_pkg("ggplot2")
need_pkg("AnnotationDbi", bioc = TRUE)
need_pkg("clusterProfiler", bioc = TRUE)
need_pkg("enrichplot", bioc = TRUE)
need_pkg("org.Hs.eg.db", bioc = TRUE)
need_pkg("org.Mm.eg.db", bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(enrichplot)
})

read_tusco_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0) {
    stop(sprintf("No data lines found in %s", path))
  }
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
    mutate(entrez_id = str_trim(entrez_id)) %>%
    filter(!is.na(entrez_id), entrez_id != "")
  if (nrow(df) == 0) {
    stop(sprintf("No non-empty Entrez IDs found in %s", path))
  }
  df
}

write_enrich_results <- function(enrich_obj, path) {
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
  ggsave(out_path, p, width = 8.5, height = 6.5, units = "in", dpi = 300)
  invisible(NULL)
}

keyword_summary <- function(enrich_obj, keywords) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(tibble(keyword = keywords, n_terms = 0L, best_padj = NA_real_, best_term = NA_character_))
  }
  res <- as_tibble(enrich_obj@result) %>%
    mutate(desc_lower = str_to_lower(Description))
  out <- lapply(keywords, function(k) {
    hits <- res %>% filter(str_detect(desc_lower, fixed(tolower(k))))
    if (nrow(hits) == 0) {
      return(tibble(keyword = k, n_terms = 0L, best_padj = NA_real_, best_term = NA_character_))
    }
    best <- hits %>% arrange(p.adjust) %>% slice(1)
    tibble(keyword = k, n_terms = nrow(hits), best_padj = best$p.adjust, best_term = best$Description)
  })
  bind_rows(out)
}

run_species <- function(species, tusco_path, orgdb_pkg) {
  tusco <- read_tusco_tsv(tusco_path)
  message(sprintf("[%s] genes=%d", species, nrow(tusco)))

  OrgDb <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))

  gene_ids <- unique(tusco$entrez_id)
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]

  ann <- AnnotationDbi::select(
    OrgDb,
    keys = gene_ids,
    keytype = "ENTREZID",
    columns = c("ENTREZID", "SYMBOL", "GENENAME")
  ) %>%
    as_tibble() %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    rename(entrez_id = ENTREZID, symbol_db = SYMBOL, genename = GENENAME)

  tusco_annot <- tusco %>%
    left_join(ann, by = "entrez_id") %>%
    mutate(symbol = if_else(is.na(symbol) | symbol == "", symbol_db, symbol)) %>%
    select(ensembl_gene, ensembl_transcript, symbol, entrez_id, genename, refseq_mrna, refseq_protein) %>%
    arrange(symbol)

  readr::write_tsv(tusco_annot, file.path(outdir, sprintf("%s_tusco_gene_annotations.tsv", species)))

  ego_bp_all <- suppressMessages(enrichGO(
    gene = gene_ids,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ))
  ego_cc_all <- suppressMessages(enrichGO(
    gene = gene_ids,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ))
  ego_mf_all <- suppressMessages(enrichGO(
    gene = gene_ids,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ))

  ego_bp_sig <- filter_enrich(ego_bp_all, padj_cutoff = 0.05)
  ego_cc_sig <- filter_enrich(ego_cc_all, padj_cutoff = 0.05)
  ego_mf_sig <- filter_enrich(ego_mf_all, padj_cutoff = 0.05)

  write_enrich_results(ego_bp_all, file.path(outdir, sprintf("%s_enrichGO_BP_all.tsv", species)))
  write_enrich_results(ego_cc_all, file.path(outdir, sprintf("%s_enrichGO_CC_all.tsv", species)))
  write_enrich_results(ego_mf_all, file.path(outdir, sprintf("%s_enrichGO_MF_all.tsv", species)))

  write_enrich_results(ego_bp_sig, file.path(outdir, sprintf("%s_enrichGO_BP.tsv", species)))
  write_enrich_results(ego_cc_sig, file.path(outdir, sprintf("%s_enrichGO_CC.tsv", species)))
  write_enrich_results(ego_mf_sig, file.path(outdir, sprintf("%s_enrichGO_MF.tsv", species)))

  save_dotplot(ego_bp_sig, sprintf("%s_GO_BP_dotplot.pdf", species), sprintf("%s TUSCO: GO BP enrichment (FDR <= 0.05)", toupper(species)))
  save_dotplot(ego_cc_sig, sprintf("%s_GO_CC_dotplot.pdf", species), sprintf("%s TUSCO: GO CC enrichment (FDR <= 0.05)", toupper(species)))
  save_dotplot(ego_mf_sig, sprintf("%s_GO_MF_dotplot.pdf", species), sprintf("%s TUSCO: GO MF enrichment (FDR <= 0.05)", toupper(species)))

  keywords <- c(
    "ribosome", "translation", "mitochond", "metabolic", "enzyme",
    "cytoskeleton", "structural", "actin", "microtubule", "extracellular matrix"
  )
  ks <- bind_rows(
    keyword_summary(ego_bp_all, keywords) %>% mutate(ontology = "BP"),
    keyword_summary(ego_cc_all, keywords) %>% mutate(ontology = "CC"),
    keyword_summary(ego_mf_all, keywords) %>% mutate(ontology = "MF")
  ) %>%
    select(ontology, keyword, n_terms, best_padj, best_term)

  readr::write_tsv(ks, file.path(outdir, sprintf("%s_keyword_summary.tsv", species)))

  invisible(list(bp = ego_bp_sig, cc = ego_cc_sig, mf = ego_mf_sig, annotations = tusco_annot))
}

human <- run_species("human", human_path, "org.Hs.eg.db")
mouse <- run_species("mouse", mouse_path, "org.Mm.eg.db")

message(sprintf("Done. Outputs written to %s", outdir))
