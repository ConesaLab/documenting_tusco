#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

find_repo_root_for_utils <- function(start = getwd(), limit = 8) {
  cur <- tryCatch(normalizePath(start, winslash = "/", mustWork = FALSE), error = function(e) start)
  for (i in seq_len(limit)) {
    if (file.exists(file.path(cur, "scripts", "figure_utils.R"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not find repo root containing scripts/figure_utils.R starting from: ", start)
}

script_path <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
script_dir <- if (length(script_path) == 1 && nzchar(script_path)) dirname(normalizePath(script_path, mustWork = FALSE)) else getwd()
repo_root <- find_repo_root_for_utils(script_dir)
source(file.path(repo_root, "scripts", "figure_utils.R"))

read_gene_table <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  df <- read.delim(con, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df
}

ensure_gene_table <- function(
  out_tsv_gz,
  species_label,
  junction_bed_gz,
  gtf_gz,
  tusco_tsv,
  min_novel_len = 80L,
  min_read_support = 1L
) {
  if (file.exists(out_tsv_gz)) return(invisible(out_tsv_gz))

  py <- Sys.which("python3")
  if (!nzchar(py)) stop("python3 not found on PATH; required to generate gene junction tables.")

  py_code <- "
import csv
import gzip
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

def ensembl_base_id(raw_id: str) -> str:
    return raw_id.split('.', 1)[0]

def parse_gtf_attrs(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_field.strip().split(';'):
        part = part.strip()
        if not part:
            continue
        key, _, rest = part.partition(' ')
        value = rest.strip().strip('\"')
        if key:
            attrs[key] = value
    return attrs

@dataclass
class GeneModel:
    chrom: str
    start: int  # BED-like start (0-based)
    end: int    # BED-like end (inclusive; matches recount junction BED)
    strand: str
    introns: Set[Tuple[int, int]]  # (start, end) using recount BED coordinates
    intron_count: int

def load_gene_models_union_introns(gtf_gz: str) -> Dict[str, GeneModel]:
    exons_by_tx: Dict[str, List[Tuple[int, int]]] = {}
    tx_meta: Dict[str, Tuple[str, str, str]] = {}  # tx -> (gene_id, chrom, strand)
    gene_meta: Dict[str, Tuple[str, str]] = {}     # gene -> (chrom, strand)
    gene_min_start: Dict[str, int] = {}
    gene_max_end: Dict[str, int] = {}

    with gzip.open(gtf_gz, 'rt') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            fields = line.rstrip('\\n').split('\\t')
            if len(fields) != 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _frame, attrs = fields
            if feature != 'exon':
                continue
            attr_map = parse_gtf_attrs(attrs)
            gene_id = ensembl_base_id(attr_map.get('gene_id', ''))
            tx_id = ensembl_base_id(attr_map.get('transcript_id', ''))
            if not gene_id or not tx_id:
                continue
            start_i = int(start)
            end_i = int(end)
            exons_by_tx.setdefault(tx_id, []).append((start_i, end_i))
            tx_meta[tx_id] = (gene_id, chrom, strand)
            gene_meta[gene_id] = (chrom, strand)
            gene_min_start[gene_id] = min(gene_min_start.get(gene_id, start_i), start_i)
            gene_max_end[gene_id] = max(gene_max_end.get(gene_id, end_i), end_i)

    gene_introns: Dict[str, Set[Tuple[int, int]]] = {}
    for tx_id, exons in exons_by_tx.items():
        meta = tx_meta.get(tx_id)
        if not meta:
            continue
        gene_id, _chrom, _strand = meta
        exons_sorted = sorted(exons)
        introns: List[Tuple[int, int]] = []
        for (_s1, e1), (s2, _e2) in zip(exons_sorted, exons_sorted[1:]):
            intron_start_1based = e1 + 1
            intron_end_1based = s2 - 1
            if intron_end_1based < intron_start_1based:
                continue
            bed_start = intron_start_1based - 1
            bed_end = intron_end_1based
            introns.append((bed_start, bed_end))
        if introns:
            gset = gene_introns.setdefault(gene_id, set())
            gset.update(introns)

    models: Dict[str, GeneModel] = {}
    for gene_id, (chrom, strand) in gene_meta.items():
        g_start_1based = gene_min_start.get(gene_id)
        g_end_1based = gene_max_end.get(gene_id)
        if g_start_1based is None or g_end_1based is None:
            continue
        intron_set = gene_introns.get(gene_id, set())
        models[gene_id] = GeneModel(
            chrom=chrom,
            start=g_start_1based - 1,
            end=g_end_1based,
            strand=strand,
            introns=intron_set,
            intron_count=len(intron_set),
        )
    return models

def load_tusco_gene_set(tsv_path: str) -> Set[str]:
    genes: Set[str] = set()
    with open(tsv_path, 'r', newline='') as handle:
        for raw_line in handle:
            line = raw_line.rstrip('\\n')
            if not line or line.startswith('#'):
                continue
            gene_id = line.split('\\t', 1)[0].strip()
            if gene_id:
                genes.add(ensembl_base_id(gene_id))
    return genes

def build_bin_index(models: Dict[str, GeneModel], bin_size: int = 1_000_000) -> Dict[str, Dict[int, List[str]]]:
    bins: Dict[str, Dict[int, List[str]]] = {}
    for gene_id, model in models.items():
        chrom_bins = bins.setdefault(model.chrom, {})
        b0 = model.start // bin_size
        b1 = model.end // bin_size
        for b in range(b0, b1 + 1):
            chrom_bins.setdefault(b, []).append(gene_id)
    return bins

@dataclass
class GeneCounts:
    known_sum: int = 0
    novel_sum: int = 0
    known_obs_junctions: int = 0
    novel_junctions: int = 0
    novel_max: int = 0

def compute_gene_counts(
    junction_bed_gz: str,
    models: Dict[str, GeneModel],
    bin_index: Dict[str, Dict[int, List[str]]],
    *,
    min_novel_len: int,
    min_read_support: int,
    bin_size: int = 1_000_000,
) -> Dict[str, GeneCounts]:
    counts: Dict[str, GeneCounts] = {gid: GeneCounts() for gid in models}

    def candidate_genes(chrom: str, start: int, end: int) -> Iterable[str]:
        chrom_bins = bin_index.get(chrom)
        if not chrom_bins:
            return []
        b0 = start // bin_size
        b1 = end // bin_size
        if b0 == b1:
            return chrom_bins.get(b0, [])
        out: List[str] = []
        for b in range(b0, b1 + 1):
            out.extend(chrom_bins.get(b, []))
        return out

    opener = gzip.open if junction_bed_gz.endswith('.gz') else open
    with opener(junction_bed_gz, 'rt') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            parts = line.rstrip('\\n').split('\\t')
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            read_support = int(parts[4])
            strand = parts[5]
            if read_support < min_read_support:
                continue

            cand = candidate_genes(chrom, start, end)
            if not cand:
                continue

            seen: Optional[Set[str]] = None
            for gid in cand:
                if seen is None:
                    seen = set()
                elif gid in seen:
                    continue
                seen.add(gid)

                model = models[gid]
                if strand != model.strand:
                    continue
                if start < model.start or end > model.end:
                    continue

                gene_counts = counts[gid]
                if (start, end) in model.introns:
                    gene_counts.known_sum += read_support
                    gene_counts.known_obs_junctions += 1
                else:
                    if (end - start) < min_novel_len:
                        continue
                    gene_counts.novel_sum += read_support
                    gene_counts.novel_junctions += 1
                    if read_support > gene_counts.novel_max:
                        gene_counts.novel_max = read_support

    return counts

def safe_ratio(num: float, den: float) -> float:
    if den <= 0:
        return float('nan')
    return num / den

def write_gene_table(out_tsv_gz: str, models: Dict[str, GeneModel], counts: Dict[str, GeneCounts], tusco_genes: Set[str]) -> None:
    os.makedirs(os.path.dirname(out_tsv_gz), exist_ok=True)
    with gzip.open(out_tsv_gz, 'wt', newline='') as handle:
        writer = csv.writer(handle, delimiter='\\t')
        writer.writerow([
            'gene_id','chrom','start','end','strand','is_tusco','annotated_intron_count',
            'known_sum_reads','known_observed_junctions','known_mean_reads_per_intron',
            'novel_sum_reads','novel_junctions','novel_max_reads',
            'novel_sum_over_known_sum','novel_max_over_known_mean'
        ])

        for gene_id, model in sorted(models.items()):
            c = counts[gene_id]
            known_mean = safe_ratio(c.known_sum, model.intron_count)
            writer.writerow([
                gene_id,
                model.chrom,
                model.start,
                model.end,
                model.strand,
                1 if gene_id in tusco_genes else 0,
                model.intron_count,
                c.known_sum,
                c.known_obs_junctions,
                f\"{known_mean:.6g}\" if not math.isnan(known_mean) else 'NA',
                c.novel_sum,
                c.novel_junctions,
                c.novel_max,
                f\"{safe_ratio(c.novel_sum, c.known_sum):.6g}\",
                f\"{safe_ratio(c.novel_max, known_mean):.6g}\" if known_mean > 0 else 'NA',
            ])

def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit('usage: script.py <gtf.gz> <junctions.bed.gz> <tusco.tsv> <out.tsv.gz>')
    gtf_gz, junction_bed_gz, tusco_tsv, out_tsv_gz = sys.argv[1:]
    min_novel_len = int(os.environ.get('TUSCO_MIN_NOVEL_LEN', '80'))
    min_read_support = int(os.environ.get('TUSCO_MIN_READ_SUPPORT', '1'))

    models = load_gene_models_union_introns(gtf_gz)
    tusco_set = load_tusco_gene_set(tusco_tsv)
    bin_index = build_bin_index(models)
    counts = compute_gene_counts(
        junction_bed_gz,
        models,
        bin_index,
        min_novel_len=min_novel_len,
        min_read_support=min_read_support,
    )
    write_gene_table(out_tsv_gz, models, counts, tusco_set)

if __name__ == '__main__':
    main()
"

  py_path <- tempfile(pattern = "splice_gene_table_", fileext = ".py")
  writeLines(py_code, py_path, useBytes = TRUE)
  on.exit(unlink(py_path), add = TRUE)

  dir.create(dirname(out_tsv_gz), recursive = TRUE, showWarnings = FALSE)

  env <- c(
    paste0("TUSCO_MIN_NOVEL_LEN=", as.integer(min_novel_len)),
    paste0("TUSCO_MIN_READ_SUPPORT=", as.integer(min_read_support))
  )

  cmd <- c(py_path, gtf_gz, junction_bed_gz, tusco_tsv, out_tsv_gz)
  status <- system2(py, cmd, env = env)
  if (!identical(status, 0L)) stop("Failed to generate gene table: ", out_tsv_gz)
  invisible(out_tsv_gz)
}

prep <- function(df, species) {
  df$species <- species
  df$is_tusco <- as.integer(df$is_tusco)
  df$annotated_intron_count <- as.integer(df$annotated_intron_count)
  df$known_observed_junctions <- as.integer(df$known_observed_junctions)
  df$novel_junctions <- as.integer(df$novel_junctions)
  df$novel_max_over_known_mean <- suppressWarnings(as.numeric(df$novel_max_over_known_mean))

  df <- df[df$annotated_intron_count >= 1, , drop = FALSE]
  df <- df[(df$known_observed_junctions + df$novel_junctions) >= 1, , drop = FALSE]
  df <- df[is.finite(df$novel_max_over_known_mean) & df$novel_max_over_known_mean >= 0, , drop = FALSE]
  df$novel_max_over_known_mean_plot <- pmax(df$novel_max_over_known_mean, 1e-4)

  df$group <- ifelse(df$is_tusco == 1, "TUSCO", "non-TUSCO (junction genes)")
  df$line_id <- paste0(df$species, " | ", df$group)
  df
}

out_dir <- file.path(repo_root, "reviewer_response", "review_plots", "splice")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hsa_path <- file.path(out_dir, "hsa_gene_junction_read_counts.tsv.gz")
mmu_path <- file.path(out_dir, "mmu_gene_junction_read_counts.tsv.gz")

ensure_gene_table(
  hsa_path,
  species_label = "hsa",
  junction_bed_gz = file.path(repo_root, "src", "tusco_selector", "data", "hsa", "splicing_tss", "recount3.pass1V1.bed.gz"),
  gtf_gz = file.path(repo_root, "src", "tusco_selector", "data", "hsa", "annotation", "gencode.v49.annotation.gtf.gz"),
  tusco_tsv = file.path(repo_root, "data", "processed", "tusco", "tusco_human.tsv")
)
ensure_gene_table(
  mmu_path,
  species_label = "mmu",
  junction_bed_gz = file.path(repo_root, "src", "tusco_selector", "data", "mmu", "splicing_tss", "recount3.pass1.bed.gz"),
  gtf_gz = file.path(repo_root, "src", "tusco_selector", "data", "mmu", "annotation", "gencode.vM38.annotation.gtf.gz"),
  tusco_tsv = file.path(repo_root, "data", "processed", "tusco", "tusco_mouse.tsv")
)

hsa <- prep(read_gene_table(hsa_path), "hsa")
mmu <- prep(read_gene_table(mmu_path), "mmu")
df <- rbind(hsa, mmu)

alpha_df <- data.frame(
  species = c("hsa", "mmu"),
  alpha = c(0.01, 0.05)
)

df$species <- factor(df$species, levels = c("hsa", "mmu"))

color_map <- c(
  "hsa | TUSCO" = TUSCO_COLORS$human_tusco,
  "mmu | TUSCO" = TUSCO_COLORS$mouse_tusco,
  "hsa | non-TUSCO (junction genes)" = TUSCO_COLORS$human_gencode,
  "mmu | non-TUSCO (junction genes)" = TUSCO_COLORS$mouse_gencode
)

p <- ggplot(df, aes(x = novel_max_over_known_mean_plot, color = line_id)) +
  stat_ecdf(geom = "step", linewidth = 0.4) +
  geom_vline(
    data = alpha_df,
    aes(xintercept = alpha),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.25,
    color = "black"
  ) +
  facet_wrap(~species, nrow = 1, scales = "fixed") +
  scale_x_log10(limits = c(1e-4, 1e1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = color_map, breaks = names(color_map)) +
  labs(
    x = "max(novel junction reads) / mean(annotated junction reads)",
    y = "ECDF across genes",
    caption = "Ratios of 0 are plotted at 1e-4 for log scaling."
  ) +
  theme_tusco(base_size = 7) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 6)
  )

out_path <- file.path(out_dir, "splice_ecdf_max_novel_over_known_mean_hsa_mmu.png")
out_path_pdf <- file.path(out_dir, "splice_ecdf_max_novel_over_known_mean_hsa_mmu.pdf")
ggsave(out_path, p, width = 13, height = 4.5, units = "in", dpi = 300)
ggsave(out_path_pdf, p, width = 13, height = 4.5, units = "in", device = cairo_pdf)
