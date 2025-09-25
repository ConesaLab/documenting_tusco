############################################################
#  Figure S4: Better TUSCO bar-plots (combined single PDF)
#  - Local-only mode: reads/writes ONLY within figs/fig-s4
#  - Expects optional inputs under figs/fig-s4/data/
#  - Outputs PDFs to figs/fig-s4/plot/ and TSVs to figs/fig-s4/tsv/
############################################################

# Load packages gracefully; if missing, log and exit without error
safe_load <- function(pkgs){
  missing <- character(0)
  for(p in pkgs){
    ok <- suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE)))
    if(!isTRUE(ok)) missing <- c(missing, p)
  }
  if(length(missing)>0){
    message("Missing required packages: ", paste(missing, collapse=", "))
    message("Skipping plot generation. Please install the missing packages.")
    quit(save = "no", status = 0)
  }
}
safe_load(c("tidyverse","cowplot","patchwork","glue"))

# silence NSE notes for CRAN/lintr
utils::globalVariables(c(
  "%>%","structural_category","associated_gene","subcategory","final_label",
  "eval","sample","n","group","TP","PTP","FP","FN","cat"
))

## ------------------------------------------------------------------
## CONFIGURATION (paths relative to this script's location)
##  - Strictly confine I/O to figs/fig-s4
## ------------------------------------------------------------------
# Resolve script directory robustly (works with Rscript and source)
args <- commandArgs(trailingOnly = FALSE)
script_path <- NA_character_
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
}
if (is.na(script_path) || !nzchar(script_path)) {
  script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
}
if (is.na(script_path) || !nzchar(script_path)) {
  # fallback: assume repo cwd when running interactively
  script_path <- file.path(getwd(), "figs", "supp-fig-04", "code", "fig-s4.R")
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))

# Figure directory (local sandbox root): figs/supp-fig-04
fig_dir  <- normalizePath(file.path(script_dir, ".."))

# Shared data directory under repo (read-only inputs)
root_dir       <- normalizePath(file.path(fig_dir, "..", ".."))
base_data_dir  <- file.path(root_dir, "data", "raw", "lrgasp", "tusco_novel_evl")
tusco_data_dir <- file.path(root_dir, "data", "processed", "tusco")

# Output directories (local writes only)
out_dir  <- file.path(fig_dir, "plots")               # PDFs here
tsv_dir  <- file.path(fig_dir, "tables")              # TSVs here
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tsv_dir,  showWarnings = FALSE, recursive = TRUE)

# Execution log (capture stdout + messages)
log_file <- file.path(fig_dir, "run.log")
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message");
  sink();
  try(close(log_con), silent = TRUE)
}, add = TRUE)

message("[fig-s4] Starting run at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("[fig-s4] Using shared data dir: ", base_data_dir)
message("[fig-s4] Output plot dir: ", out_dir)
message("[fig-s4] Output tsv dir: ", tsv_dir)

# Expand y-axis headroom to reduce visual bar height (no category removed)
# Increase >1 to make bars appear shorter relative to the panel height
y_limit_expand_factor <- 1.25

sample_order <- c(
  "wtc11_captrap_pacbio","wtc11_cdna_pacbio",
  "wtc11_cdna_ont","wtc11_captrap_ont",
  "es_captrap_pacbio","es_cdna_pacbio",
  "es_cdna_ont","es_captrap_ont"
)

pipeline_specs <- tribble(
  ~label,      ~dir,
  "Bambu",     "bambu_sq3",
  "StringTie", "stringtie_sq3",
  "Flair",     "flair_sq3",
  "IsoSeq",    "isoseq_sq3"
)

cat_levels  <- c("RM","Alternative 3'end","Alternative 5'end","Alternative 3'5'end",
                 "ISM","NIC","NNC","Genic Intron","Genic Genomic","Antisense",
                 "Fusion","Intergenic","Missing")
sirv_pal <- c(
  "RM" = '#c4e1f2', "Alternative 3'end"='#02314d',
  "Alternative 5'end"='#7ccdfc',"Alternative 3'5'end"='#0e5a87',
  "ISM"="#FC8D59","NIC"="#78C679","NNC"="#EE6A50",
  "Genic Intron"="#41B6C4","Genic Genomic"="#969696",
  "Antisense"="#66C2A4","Fusion"="goldenrod1",
  "Intergenic"="darksalmon","Missing"="grey70"
)

tp_grp  <- c("RM")
ptp_grp <- c("Alternative 3'end","Alternative 5'end","Alternative 3'5'end","ISM")
fp_grp  <- c("NIC","NNC","Genic Intron","Genic Genomic","Antisense","Fusion","Intergenic")

## ------------------------------------------------------------------
## Nature-style theme (matching figure1c_figure-s1.R)
## ------------------------------------------------------------------
nature_theme <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 7, face = "bold"),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    legend.key.size = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.background = element_rect(color = "black", linewidth = 0.5)
  )

## ------------------------------------------------------------------
## helper that loads + aggregates one classification file
##   - Confined to local files only
## ------------------------------------------------------------------
load_one <- function(fp, tusco_ref){
  if(!file.exists(fp)) return(tibble::tibble())
  d <- readr::read_tsv(fp, show_col_types = FALSE, progress = FALSE)
  if(nrow(d)==0) return(tibble::tibble())
  d <- d %>%
        dplyr::filter(structural_category!="fusion") %>%            # unpack fusions
        dplyr::bind_rows(d %>% dplyr::filter(structural_category=="fusion") %>%
                    tidyr::separate_rows(associated_gene, sep="_")) %>%
        dplyr::mutate(associated_gene=stringr::str_remove(associated_gene,"\\.[0-9]+$")) %>%
        dplyr::filter(associated_gene %in% tusco_ref)               # keep only TUSCO genes
  
  cat_map <- function(subcat,struct){
    if(struct=="full-splice_match" && subcat=="reference_match")       "RM"
    else if(struct=="full-splice_match" && subcat=="alternative_3end") "Alternative 3'end"
    else if(struct=="full-splice_match" && subcat=="alternative_5end") "Alternative 5'end"
    else if(struct=="full-splice_match" && subcat=="alternative_3end5end") "Alternative 3'5'end"
    else if(subcat=="mono-exon_by_intron_retention") "ISM"  # Treat mono-exon_by_intron_retention as ISM (PTP)
    else if(struct=="incomplete-splice_match") "ISM"
    else if(struct=="novel_in_catalog")         "NIC"
    else if(struct=="novel_not_in_catalog")     "NNC"
    else if(struct=="genic_intron")             "Genic Intron"
    else if(struct=="genic")                    "Genic Genomic"
    else if(struct=="antisense")                "Antisense"
    else if(struct=="fusion")                   "Fusion"
    else if(struct=="intergenic")               "Intergenic"
    else NA_character_
  }
  d %>%
    dplyr::mutate(final_label = mapply(cat_map, subcategory, structural_category)) %>%
    dplyr::filter(!is.na(final_label))
}

## ------------------------------------------------------------------
## Main loop – build plots for each pipeline (Ref and Novel panels)
## ------------------------------------------------------------------
plots_for_grid <- list()
tsv_stacked_list <- list()
tsv_stats_list   <- list()
for(i in seq_len(nrow(pipeline_specs))){
  pinfo <- pipeline_specs[i,]
  message("Processing ", pinfo$label)
  
  # determine which sample prefixes really exist for this pipeline
  existing_samples <- sample_order[ sapply(sample_order, function(sp){
      ref_dir   <- file.path(base_data_dir, pinfo$dir, "ref_evl",  sp)
      novel_dir <- file.path(base_data_dir, pinfo$dir, "novel_evl", sp)
      dir.exists(ref_dir) || dir.exists(novel_dir)
  }) ]

  if(length(existing_samples)==0){
    message("No sample directories found for ", pinfo$label, " in ", base_data_dir)
    next
  }

  df_all <- purrr::map_dfr(existing_samples, function(sprefix){
    species <- ifelse(str_detect(sprefix,"^wtc11"),"human","mouse")
    species_dir <- if(species=="human") "hsa" else "mmu"
    tusco_file <- file.path(tusco_data_dir, species_dir,
                            if(species=="human") "tusco_human_multi_exon.tsv" else "tusco_mouse_multi_exon.tsv")
    if(file.exists(tusco_file)){
      tusco_df <- readr::read_tsv(tusco_file, show_col_types = FALSE, progress = FALSE,
                           col_types = readr::cols(.default="c"),
                           comment = "#", col_names = FALSE) %>%
                  setNames(c("ensembl","transcript","gene_name","entrez","refseq","protein"))
      tusco_ref <- tusco_df %>%
                   dplyr::select(ensembl, refseq, gene_name) %>%
                   unlist(use.names = FALSE) %>%
                   unique()
    } else {
      message("  TUSCO file not found: ", tusco_file)
      tusco_df  <- tibble::tibble()
      tusco_ref <- character(0)
    }
    
    # build expected paths
    dbase <- file.path(base_data_dir, pinfo$dir)
    sqanti_dir <- "sqanti3_out"
    ref_fp   <- file.path(dbase,"ref_evl",  sprefix, sqanti_dir, "isoforms_classification.txt")
    novel_fp <- file.path(dbase,"novel_evl",sprefix, sqanti_dir, "isoforms_classification.txt")
    if(pinfo$dir=="stringtie_sq3"){
      ref_fp   <- file.path(dbase,"ref_evl",  sprefix, sqanti_dir,
                            glue("{sprefix}_stringtie_ref_sqanti_classification.txt"))
      novel_fp <- file.path(dbase,"novel_evl",sprefix, sqanti_dir,
                            glue("{sprefix}_stringtie_tusco_sqanti_classification.txt"))
    }
    if(pinfo$dir %in% c("bambu_sq3","flair_sq3")){
      ref_fp   <- file.path(dbase,"ref_evl",  sprefix, sqanti_dir,
                            glue("{sprefix}_sqanti_classification.txt"))
      novel_fp <- file.path(dbase,"novel_evl",sprefix, sqanti_dir,
                            glue("{sprefix}_sqanti_classification.txt"))
    }
    
    # load both evaluation types
    ref_set   <- load_one(ref_fp,   tusco_ref) %>% dplyr::mutate(eval="Ref Evl")
    novel_set <- load_one(novel_fp, tusco_ref) %>% dplyr::mutate(eval="Novel Evl")

    comb <- bind_rows(ref_set, novel_set)

    # decide primary ID type based on what is present in data (Ensembl > RefSeq > gene_name)
    patterns <- list(ensembl="^(ENSG|ENSMUSG)", refseq="^(NM_|NR_|NP_)")
    top_id <- if(any(stringr::str_detect(comb$associated_gene, patterns$ensembl))) "ensembl" 
              else if(any(stringr::str_detect(comb$associated_gene, patterns$refseq))) "refseq" 
              else "gene_name"

    tusco_ids <- tusco_df[[top_id]] %>% stats::na.omit() %>% unique()

    # add Missing rows per evaluation
    if("associated_gene" %in% colnames(comb)){
      miss_rows <- comb %>%
                   dplyr::group_by(eval) %>%
                   dplyr::summarise(detected = dplyr::n_distinct(associated_gene), .groups="drop") %>%
                   dplyr::mutate(n_missing = pmax(length(tusco_ids) - detected, 0))
    } else {
      miss_rows <- tibble::tibble(eval = c("Ref Evl","Novel Evl"), detected = 0, n_missing = length(tusco_ids))
    }
    miss_rows <- miss_rows %>%
                   dplyr::mutate(final_label = "Missing", sample = sprefix) %>%
                   dplyr::select(sample, eval, final_label, n = n_missing)

    comb <- comb %>% dplyr::mutate(sample = sprefix, n = 1)

    # return full set rows; Missing rows will be counted later
    bind_rows(comb, miss_rows)
  })
  
  if(nrow(df_all)==0){
    message("No data for ", pinfo$label)
    next
  }
  
  ## ----------------------------------------------------------------
  ##  If IsoSeq, print FP details for Ref/Novel to console
  ## ----------------------------------------------------------------
  if (identical(as.character(pinfo$dir), "isoseq_sq3")) {
    iso_fp <- df_all %>%
      dplyr::filter(eval %in% c("Ref Evl", "Novel Evl"),
                    final_label %in% fp_grp,
                    !is.na(associated_gene)) %>%
      dplyr::mutate(sample = as.character(sample),
                    associated_gene = as.character(associated_gene)) %>%
      dplyr::group_by(eval, sample, associated_gene, final_label) %>%
      dplyr::summarise(n = sum(n), .groups = "drop") %>%
      dplyr::arrange(eval, sample, dplyr::desc(n), associated_gene)
    message("IsoSeq – false positives (eval, sample, gene, label, n):")
    if (nrow(iso_fp) > 0) {
      print(iso_fp, n = nrow(iso_fp))
    } else {
      message("  None")
    }
  }
  
  ## summarise → stacked barplot frame + TP/PTP/FP/FN frame
  df_cnt <- df_all %>%
              dplyr::group_by(sample, eval, final_label) %>%
              dplyr::summarise(n = sum(n), .groups="drop") %>%
              dplyr::mutate(final_label=factor(final_label, levels=cat_levels))
  
  df_stacked <- df_cnt
  
  # Recompute TP/PTP/FP/FN groups with extended TP rule for multi-exon only
  df_stats <- df_all %>%
    dplyr::mutate(
      mono_exon_close50 = !is.na(ref_exons) & ref_exons == 1 &
                          !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
                          abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50,
      # Rule 3: For ref_length > 3000, FSM with TSS and TTS within 100bp is TP
      long_transcript_close100 = !is.na(ref_length) & ref_length > 3000 & 
                                 structural_category == "full-splice_match" &
                                 !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
                                 abs(diff_to_TSS) <= 100 & abs(diff_to_TTS) <= 100,
      group = dplyr::case_when(
        final_label == "Missing" ~ "FN",
        final_label %in% tp_grp ~ "TP",
        mono_exon_close50 & structural_category == "full-splice_match" ~ "TP",
        long_transcript_close100 ~ "TP",
        final_label %in% ptp_grp ~ "PTP",
        final_label %in% fp_grp  ~ "FP",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::group_by(sample, eval, group) %>%
    dplyr::summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
    dplyr::mutate(total = TP + PTP + FP + FN)
  
  # ensure all expected columns exist even if absent in data
  for(col in c("TP","PTP","FP","FN")){
    if(!(col %in% colnames(df_stats))) df_stats[[col]] <- 0
  }
  if(!("total" %in% colnames(df_stats))){
    df_stats$total <- df_stats$TP + df_stats$PTP + df_stats$FP + df_stats$FN
  }

  # dynamic label offset per evaluation type to avoid overlap with bars
  label_offset_tbl <- df_stats %>%
    dplyr::group_by(eval) %>%
    dplyr::summarise(max_total = max(total, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(label_offset = pmax(12, ceiling(0.30 * max_total)))
  df_stats <- df_stats %>%
    dplyr::left_join(label_offset_tbl, by = "eval") %>%
    dplyr::mutate(label_y = total + label_offset)
  
  # pretty x-axis labels (three lines: species, prep, platform)
  format_sample_label <- function(s){
    parts <- strsplit(s, "_")[[1]]
    if(length(parts) == 3){
      species  <- parts[1]
      prep     <- parts[2]
      platform <- parts[3]
      sp  <- if (species == "wtc11") "WTC11" else toupper(species)  # ES
      pr  <- dplyr::case_when(prep == "cdna" ~ "cDNA",
                               prep == "captrap" ~ "CapTrap",
                               TRUE ~ prep)
      pl  <- dplyr::case_when(platform == "pacbio" ~ "PB",
                               platform == "pb" ~ "PB",
                               platform == "ont" ~ "ONT",
                               TRUE ~ toupper(platform))
      return(paste(sp, pr, pl, sep = "\n"))
    }
    # fallback mapping
    s <- stringr::str_replace_all(s, c("cdna" = "cDNA", "captrap" = "CapTrap", "pacbio" = "PB", "pb" = "PB", "ont" = "ONT"))
    gsub("_", "\n", s)
  }
  sample_labels <- setNames(vapply(sample_order, format_sample_label, character(1)), sample_order)

  build_one_plot <- function(eval_level){
    dst <- df_stacked %>% filter(eval==eval_level)
    if(nrow(dst)==0) return(NULL)
    # Nudge label text slightly higher for panels a,c,d,e,f
    extra_offset <- if ((identical(as.character(pinfo$label), "Bambu") && eval_level == "Novel Evl") ||
                         (identical(as.character(pinfo$label), "StringTie") && eval_level %in% c("Novel Evl","Ref Evl")) ||
                         (identical(as.character(pinfo$label), "Flair") && eval_level %in% c("Novel Evl","Ref Evl"))) 4 else 0
    df_stats_mod <- df_stats %>%
      dplyr::filter(eval == eval_level) %>%
      dplyr::mutate(label_y = label_y + extra_offset)
    ggplot2::ggplot(dst,
           ggplot2::aes(x=factor(sample, levels = sample_order),
               y=n, fill=final_label)) +
      ggplot2::geom_col(width=.55, colour="black", linewidth=.2) +
      ggplot2::scale_fill_manual(values = sirv_pal,
                                 name = "SQANTI category",
                                 limits = names(sirv_pal),
                                 drop = FALSE) +
      {
        ylim_max <- df_stats_mod %>%
          dplyr::summarise(mx = max(label_y, na.rm = TRUE)) %>% dplyr::pull(mx)
        if(!is.finite(ylim_max) || is.na(ylim_max)) ylim_max <- 10
        ylim_max <- ylim_max * y_limit_expand_factor
        ggplot2::scale_y_continuous(limits = c(0, ylim_max), expand = ggplot2::expansion(mult = c(0,0)))
      } +
      ggplot2::scale_x_discrete(labels = sample_labels) +
      nature_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.85),
                     plot.margin = ggplot2::margin(5,5,12,5),
                     plot.caption = ggplot2::element_text(face = "bold", hjust = .5, size = 6)) +
      ggplot2::labs(x=NULL, y="Count",
           caption = glue::glue("{pinfo$label} – {ifelse(eval_level=='Ref Evl','Reference','Novel')}") ) +
      ggplot2::geom_text(data = df_stats_mod,
                ggplot2::aes(x = sample, y = label_y,
                    label = glue::glue("TP={TP}\nPTP={PTP}\nFP={FP}\nFN={FN}")),
                size = 2.4, fontface = "plain", lineheight = .9,
                hjust = 0.5, inherit.aes = FALSE) +
      ggplot2::coord_cartesian(clip = "off")
  }

  plot_ref   <- build_one_plot("Ref Evl")
  plot_novel <- build_one_plot("Novel Evl")

  for(pp in list(plot_novel, plot_ref)){
    if(!is.null(pp)){
       plots_for_grid[[length(plots_for_grid)+1]] <- pp
    }
  }

  # Accumulate TSV data
  if(nrow(df_stacked) > 0){
    tsv_stacked_list[[length(tsv_stacked_list)+1]] <- df_stacked %>%
      dplyr::mutate(pipeline = as.character(pinfo$label)) %>%
      dplyr::select(pipeline, sample, eval, final_label, n)
  }
  if(nrow(df_stats) > 0){
    tsv_stats_list[[length(tsv_stats_list)+1]] <- df_stats %>%
      dplyr::mutate(pipeline = as.character(pinfo$label)) %>%
      dplyr::select(pipeline, sample, eval, TP, PTP, FP, FN, total, label_y)
  }
}

# ------------------------------------------------------------------
#  Build composite figure with legend (single output PDF)
# ------------------------------------------------------------------
if(length(plots_for_grid)>0){
  # 1) Arrange panels only (no legends)
  plots_no_legend <- lapply(plots_for_grid, function(p) p + ggplot2::theme(legend.position = "none"))
  combined_grid <- cowplot::plot_grid(
    plotlist = plots_no_legend,
    ncol = 2,
    labels = letters[seq_along(plots_no_legend)],
    label_size = 7,
    label_fontface = "bold"
  )

  # 2) Build a standalone legend (2 rows) without relying on ggplot guides
  legend_levels <- names(sirv_pal)
  n_items <- length(legend_levels)
  n_rows <- 2L
  n_cols <- ceiling(n_items / n_rows)
  # Increase spacing between legend columns/rows to avoid text overlap
  legend_col_spacing <- 5  
  legend_row_spacing <- 1.5
  legend_tile_size <- 0.45  # square size for legend swatches (increase to widen)
  legend_layout <- tibble::tibble(
    cat = factor(legend_levels, levels = legend_levels),
    idx = seq_len(n_items),
    row = ceiling(idx / n_cols),
    col = idx - (row - 1L) * n_cols
  ) %>%
  dplyr::mutate(
    # place items on a 0-indexed grid so spacing reflects rows/cols exactly
    row_pos = (row - 1L) * legend_row_spacing,
    col_pos = (col - 1L) * legend_col_spacing
  )
  legend_plot_only <- ggplot2::ggplot(legend_layout, ggplot2::aes(x = col_pos, y = row_pos)) +
    ggplot2::geom_tile(ggplot2::aes(fill = cat), width = legend_tile_size, height = legend_tile_size, color = "black", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = cat), hjust = 0, nudge_x = 0.7, size = 2.2) +
    ggplot2::scale_fill_manual(values = sirv_pal, guide = "none", limits = legend_levels, drop = FALSE) +
    {
      # symmetric padding so legend sits centered in its row
      x_pad <- legend_col_spacing * 0.8
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(x_pad, x_pad)))
    } +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(add = c(0.35, 0.35))) +
    ggplot2::theme_void(base_family = "Helvetica", base_size = 7) +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 6, 3, 6, unit = "pt"),
                   panel.border = ggplot2::element_rect(color = "black", linewidth = 0.5, fill = NA)) +
    ggplot2::coord_fixed(ratio = 1, clip = "off")

  # 3) Save separate PDFs: panels-with-space, legend-only, and final merged
  save_pdf <- function(filename, plot, width, height, units = "mm"){
    ok <- TRUE
    # Prefer Cairo if available
    dev_fun <- NULL
    if("cairo_pdf" %in% getNamespaceExports("grDevices")) dev_fun <- grDevices::cairo_pdf
    tryCatch({
      if(!is.null(dev_fun)){
        ggplot2::ggsave(filename, plot, width = width, height = height, units = units, device = dev_fun)
      } else {
        ggplot2::ggsave(filename, plot, width = width, height = height, units = units, device = "pdf")
      }
    }, error = function(e){ ok <<- FALSE; message("  Failed to save PDF ", basename(filename), ": ", conditionMessage(e)) })
    invisible(ok)
  }
  legend_rel <- 0.12  # portion of total height reserved for the legend in the final figure

  # 3a) Panels with reserved space at the bottom
  empty_spacer <- ggplot2::ggplot() + ggplot2::theme_void()
  panels_with_space <- cowplot::plot_grid(combined_grid, empty_spacer, ncol = 1,
                                          rel_heights = c(1, legend_rel))
  out_panels <- file.path(out_dir, "fig-s4_panels.pdf")
  if (save_pdf(out_panels, panels_with_space, width = 183, height = 170, units = "mm")){
    message("  → saved ", out_panels)
  }

  # 3b) Legend-only PDF sized to the reserved height
  out_legend <- file.path(out_dir, "fig-s4_legend.pdf")
  if (save_pdf(out_legend, legend_plot_only, width = 183, height = 170 * legend_rel, units = "mm")){
    message("  → saved ", out_legend)
  }

  # 3c) Final combined single-page figure
  final_figure <- cowplot::plot_grid(combined_grid, legend_plot_only, ncol = 1, rel_heights = c(1, legend_rel))
  out_pdf <- file.path(out_dir, "fig-s4.pdf")
  if (save_pdf(out_pdf, final_figure, width = 183, height = 170, units = "mm")){
    message("  → saved ", out_pdf)
  }

  # 4) Write TSVs matching each PDF
  combined_stacked <- dplyr::bind_rows(tsv_stacked_list)
  combined_stats   <- dplyr::bind_rows(tsv_stats_list)

  # helper to safely write TSV
  safe_write_tsv <- function(df, path){
    ok <- TRUE
    tryCatch(readr::write_tsv(df, path), error = function(e){ ok <<- FALSE; message("  Failed to write TSV ", basename(path), ": ", conditionMessage(e)) })
    invisible(ok)
  }

  # Panels TSV (stacked + stats, flagged by dataset)
  panels_tsv <- file.path(tsv_dir, "fig-s4_panels.tsv")
  if(nrow(combined_stacked) > 0 || nrow(combined_stats) > 0){
    df_panels <- dplyr::bind_rows(
      combined_stacked %>% dplyr::mutate(dataset = "stacked") %>% dplyr::relocate(dataset),
      combined_stats   %>% dplyr::mutate(dataset = "stats")   %>% dplyr::relocate(dataset)
    ) %>% dplyr::mutate(figure_id = "fig-s4_panels", .before = 1L)
    if(safe_write_tsv(df_panels, panels_tsv)) message("  → wrote ", panels_tsv)
  } else {
    message("  No data for panels TSV; skipping write.")
  }

  # Legend TSV (category-color mapping and grid layout)
  legend_tsv <- file.path(tsv_dir, "fig-s4_legend.tsv")
  legend_tbl <- legend_layout %>%
    dplyr::transmute(figure_id = "fig-s4_legend", category = as.character(cat), row, col, row_pos, col_pos,
                     color = unname(sirv_pal[as.character(cat)]),
                     legend_tile_size = legend_tile_size,
                     legend_row_spacing = legend_row_spacing,
                     legend_col_spacing = legend_col_spacing)
  if(nrow(legend_tbl) > 0){
    if(safe_write_tsv(legend_tbl, legend_tsv)) message("  → wrote ", legend_tsv)
  }

  # Final figure TSV (same data as panels; mark figure_id)
  final_tsv <- file.path(tsv_dir, "fig-s4.tsv")
  if(nrow(combined_stacked) > 0 || nrow(combined_stats) > 0){
    df_final <- dplyr::bind_rows(
      combined_stacked %>% dplyr::mutate(dataset = "stacked") %>% dplyr::relocate(dataset),
      combined_stats   %>% dplyr::mutate(dataset = "stats")   %>% dplyr::relocate(dataset)
    ) %>% dplyr::mutate(figure_id = "fig-s4", .before = 1L)
    if(safe_write_tsv(df_final, final_tsv)) message("  → wrote ", final_tsv)
  }
} else {
  message("No plots generated; check shared data dir: ", base_data_dir)
}

# Ensure TSV exists for every existing PDF (metadata-only if no data)
existing_pdfs <- list.files(out_dir, pattern = "\\.pdf$", full.names = FALSE)
for (pdf_bn in existing_pdfs) {
  tsv_bn <- sub("\\.pdf$", ".tsv", pdf_bn)
  tsv_fp <- file.path(tsv_dir, tsv_bn)
  if (!file.exists(tsv_fp)) {
    meta <- tibble::tibble(
      figure_id = sub("\\.pdf$", "", pdf_bn),
      status = "no_underlying_data_available",
      message = paste0("No underlying data processed; check ", base_data_dir, ". TSV generated with metadata only."),
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    tryCatch({ readr::write_tsv(meta, tsv_fp); message("  → wrote ", tsv_fp) },
             error = function(e) message("  Failed to write metadata TSV for ", pdf_bn, ": ", conditionMessage(e)))
  }
}

# === EXTRACTION FOR PAPER VALUES ===
# Extract TUSCO-novel performance statistics for manuscript updates
# Write inside the repository root (stay within project sandbox)
extraction_file <- file.path(root_dir, "manuscript", "extraction_values.tsv")

# Extract tool performance averages if we have stats data
extraction_data_list <- list()

if (exists("combined_stats") && nrow(combined_stats) > 0) {
  # Get unique tools/pipelines
  tools <- unique(combined_stats$pipeline)
  
  for (tool in tools) {
    tool_data <- combined_stats %>% filter(pipeline == tool)
    
    if (nrow(tool_data) > 0) {
      # Calculate averages for Novel Evl
      novel_data <- tool_data %>% filter(eval == "Novel Evl")
      if (nrow(novel_data) > 0) {
        avg_tp <- mean(novel_data$TP, na.rm = TRUE)
        avg_ptp <- mean(novel_data$PTP, na.rm = TRUE)
        avg_fp <- mean(novel_data$FP, na.rm = TRUE)
        avg_fn <- mean(novel_data$FN, na.rm = TRUE)
        
        extraction_data_list[[paste0(tool, "_novel_tp")]] <- data.frame(
          figure_id = "fig-s4",
          metric_name = paste0(tolower(gsub("[^A-Za-z0-9]", "_", tool)), "_novel_evl_avg_tp"),
          value = formatC(avg_tp, format="f", digits=1),
          ci_lower = NA,
          ci_upper = NA,
          notes = paste0("TUSCO-novel ", tool, " average TP (Novel Evl)"),
          stringsAsFactors = FALSE
        )
        
        extraction_data_list[[paste0(tool, "_novel_ptp")]] <- data.frame(
          figure_id = "fig-s4",
          metric_name = paste0(tolower(gsub("[^A-Za-z0-9]", "_", tool)), "_novel_evl_avg_ptp"),
          value = formatC(avg_ptp, format="f", digits=1),
          ci_lower = NA,
          ci_upper = NA,
          notes = paste0("TUSCO-novel ", tool, " average PTP (Novel Evl)"),
          stringsAsFactors = FALSE
        )
        
        extraction_data_list[[paste0(tool, "_novel_fp")]] <- data.frame(
          figure_id = "fig-s4",
          metric_name = paste0(tolower(gsub("[^A-Za-z0-9]", "_", tool)), "_novel_evl_avg_fp"),
          value = formatC(avg_fp, format="f", digits=1),
          ci_lower = NA,
          ci_upper = NA,
          notes = paste0("TUSCO-novel ", tool, " average FP (Novel Evl)"),
          stringsAsFactors = FALSE
        )
        
        extraction_data_list[[paste0(tool, "_novel_fn")]] <- data.frame(
          figure_id = "fig-s4",
          metric_name = paste0(tolower(gsub("[^A-Za-z0-9]", "_", tool)), "_novel_evl_avg_fn"),
          value = formatC(avg_fn, format="f", digits=1),
          ci_lower = NA,
          ci_upper = NA,
          notes = paste0("TUSCO-novel ", tool, " average FN (Novel Evl)"),
          stringsAsFactors = FALSE
        )
      }
      
      # Calculate averages for Ref Evl
      ref_data <- tool_data %>% filter(eval == "Ref Evl")
      if (nrow(ref_data) > 0) {
        avg_tp_ref <- mean(ref_data$TP, na.rm = TRUE)
        avg_ptp_ref <- mean(ref_data$PTP, na.rm = TRUE)
        avg_fp_ref <- mean(ref_data$FP, na.rm = TRUE)
        avg_fn_ref <- mean(ref_data$FN, na.rm = TRUE)
        
        extraction_data_list[[paste0(tool, "_ref_tp")]] <- data.frame(
          figure_id = "fig-s4",
          metric_name = paste0(tolower(gsub("[^A-Za-z0-9]", "_", tool)), "_ref_evl_avg_tp"),
          value = formatC(avg_tp_ref, format="f", digits=1),
          ci_lower = NA,
          ci_upper = NA,
          notes = paste0("TUSCO-novel ", tool, " average TP (Ref Evl)"),
          stringsAsFactors = FALSE
        )
        
        extraction_data_list[[paste0(tool, "_ref_ptp")]] <- data.frame(
          figure_id = "fig-s4",
          metric_name = paste0(tolower(gsub("[^A-Za-z0-9]", "_", tool)), "_ref_evl_avg_ptp"),
          value = formatC(avg_ptp_ref, format="f", digits=1),
          ci_lower = NA,
          ci_upper = NA,
          notes = paste0("TUSCO-novel ", tool, " average PTP (Ref Evl)"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# Combine all extraction data
if (length(extraction_data_list) > 0) {
  extraction_data <- dplyr::bind_rows(extraction_data_list)
  
  # Write or append to file
  if (!file.exists(extraction_file)) {
    readr::write_tsv(extraction_data, extraction_file)
    message("Created extraction file: ", extraction_file)
  } else {
    readr::write_tsv(extraction_data, extraction_file, append = TRUE, col_names = FALSE)
  }
  
  message("Extracted ", nrow(extraction_data), " TUSCO-novel tool performance metrics to: ", extraction_file)
} else {
  message("No TUSCO-novel data available for extraction")
}
# === END EXTRACTION ===

message("[fig-s4] Completed run at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
