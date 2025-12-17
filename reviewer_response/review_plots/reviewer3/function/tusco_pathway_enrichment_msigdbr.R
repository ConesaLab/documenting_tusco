#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && any(args %in% c("-h", "--help"))) {
  cat(
    "Usage: Rscript tusco_pathway_enrichment_msigdbr.R [--outdir DIR] [--human TSV] [--mouse TSV] [--no-install]\n",
    "\nRuns ORA using MSigDB gene sets via msigdbr (Reactome, KEGG, Hallmark).\n",
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
need_pkg("stringr")
need_pkg("ggplot2")
need_pkg("clusterProfiler", bioc = TRUE)
need_pkg("enrichplot", bioc = TRUE)
need_pkg("msigdbr")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
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
    mutate(symbol = str_trim(symbol)) %>%
    filter(!is.na(symbol), symbol != "")
  if (nrow(df) == 0) {
    stop(sprintf("No non-empty gene symbols found in %s", path))
  }
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

run_ora_all <- function(genes, term2gene, term2name = NULL) {
  genes <- unique(genes)
  if (length(genes) < 2) {
    return(NULL)
  }
  suppressMessages(clusterProfiler::enricher(
    gene = genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500
  ))
}

filter_enrich <- function(enrich_obj, padj_cutoff = 0.05) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(enrich_obj)
  }
  res <- as_tibble(enrich_obj@result) %>% filter(!is.na(p.adjust), p.adjust <= padj_cutoff)
  enrich_obj@result <- as.data.frame(res)
  enrich_obj
}

msig_term2gene <- function(species, collection, subcollection = NULL, db_species = NULL) {
  call_args <- list(species = species, collection = collection, subcollection = subcollection)
  if (!is.null(db_species)) {
    call_args$db_species <- db_species
  }
  m <- do.call(msigdbr::msigdbr, call_args)
  if (nrow(m) == 0) {
    stop(sprintf("No MSigDB entries for species=%s collection=%s subcollection=%s", species, collection, subcollection))
  }
  term2gene <- m %>% select(gs_name, gene_symbol) %>% distinct()
  term2name <- m %>% select(gs_name, gs_description) %>% distinct() %>% rename(term = gs_name, name = gs_description)
  list(term2gene = term2gene, term2name = term2name)
}

analyze_species <- function(label, tusco_path, msig_species, db_species = NULL) {
  tusco <- read_tusco_tsv(tusco_path)
  genes <- unique(tusco$symbol)
  message(sprintf("[%s] genes=%d", label, length(genes)))

  if (is.null(db_species)) {
    db_species <- "HS"
  }

  if (db_species == "HS") {
    reactome_collection <- "C2"
    reactome_subcollection <- "CP:REACTOME"
    hallmark_collection <- "H"
    hallmark_subcollection <- NULL
    kegg_collections_supported <- TRUE
  } else if (db_species == "MM") {
    # Mouse-native MSigDB uses M* collection codes
    reactome_collection <- "M2"
    reactome_subcollection <- "CP:REACTOME"
    hallmark_collection <- "MH"
    hallmark_subcollection <- NULL
    kegg_collections_supported <- FALSE
  } else {
    stop(sprintf("Unsupported db_species=%s (expected HS or MM)", db_species))
  }

  # Reactome pathways
  react <- msig_term2gene(msig_species, collection = reactome_collection, subcollection = reactome_subcollection, db_species = db_species)
  enr_react_all <- run_ora_all(genes, react$term2gene, react$term2name)
  enr_react_sig <- filter_enrich(enr_react_all, 0.05)
  write_results(enr_react_all, file.path(outdir, sprintf("%s_msigdb_reactome_all.tsv", label)))
  write_results(enr_react_sig, file.path(outdir, sprintf("%s_msigdb_reactome.tsv", label)))
  save_dotplot(enr_react_sig, sprintf("%s_MSIGDB_REACTOME_dotplot.pdf", label), sprintf("%s TUSCO: Reactome pathway enrichment (FDR <= 0.05)", toupper(label)))

  # KEGG pathways (human DB only)
  if (kegg_collections_supported) {
    kegg_legacy <- msig_term2gene(msig_species, collection = "C2", subcollection = "CP:KEGG_LEGACY", db_species = db_species)
    kegg_medicus <- msig_term2gene(msig_species, collection = "C2", subcollection = "CP:KEGG_MEDICUS", db_species = db_species)
    kegg <- list(
      term2gene = bind_rows(kegg_legacy$term2gene, kegg_medicus$term2gene) %>% distinct(),
      term2name = bind_rows(kegg_legacy$term2name, kegg_medicus$term2name) %>% distinct()
    )
    enr_kegg_all <- run_ora_all(genes, kegg$term2gene, kegg$term2name)
    enr_kegg_sig <- filter_enrich(enr_kegg_all, 0.05)
    write_results(enr_kegg_all, file.path(outdir, sprintf("%s_msigdb_kegg_all.tsv", label)))
    write_results(enr_kegg_sig, file.path(outdir, sprintf("%s_msigdb_kegg.tsv", label)))
    save_dotplot(enr_kegg_sig, sprintf("%s_MSIGDB_KEGG_dotplot.pdf", label), sprintf("%s TUSCO: KEGG pathway enrichment (FDR <= 0.05)", toupper(label)))
  } else {
    readr::write_tsv(tibble(), file.path(outdir, sprintf("%s_msigdb_kegg_all.tsv", label)))
    readr::write_tsv(tibble(), file.path(outdir, sprintf("%s_msigdb_kegg.tsv", label)))
    save_dotplot(NULL, sprintf("%s_MSIGDB_KEGG_dotplot.pdf", label), sprintf("%s TUSCO: KEGG pathway enrichment (FDR <= 0.05)", toupper(label)))
    enr_kegg_sig <- NULL
  }

  # Hallmark
  hm <- msig_term2gene(msig_species, collection = hallmark_collection, subcollection = hallmark_subcollection, db_species = db_species)
  enr_hm_all <- run_ora_all(genes, hm$term2gene, hm$term2name)
  enr_hm_sig <- filter_enrich(enr_hm_all, 0.05)
  write_results(enr_hm_all, file.path(outdir, sprintf("%s_msigdb_hallmark_all.tsv", label)))
  write_results(enr_hm_sig, file.path(outdir, sprintf("%s_msigdb_hallmark.tsv", label)))
  save_dotplot(enr_hm_sig, sprintf("%s_MSIGDB_HALLMARK_dotplot.pdf", label), sprintf("%s TUSCO: Hallmark enrichment (FDR <= 0.05)", toupper(label)))

  invisible(list(reactome = enr_react_sig, kegg = enr_kegg_sig, hallmark = enr_hm_sig))
}

human <- analyze_species("human", human_path, "Homo sapiens", db_species = "HS")
mouse <- analyze_species("mouse", mouse_path, "mouse", db_species = "MM")

message(sprintf("Done. Outputs written to %s", outdir))
