###############################################################################
# Load libraries
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

###############################################################################
# 1) Define file paths and pipelines
###############################################################################
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
  }
  return(candidates[[1]])
}

sirv_gtf_file <- resolve_path(c(
  file.path("..","..","..","data","raw","spike-ins","lrgasp_sirvs.gtf"),
  file.path("..","data","spike-ins","lrgasp_sirvs.gtf")
))
tusco_tsv_file <- resolve_path(c(
  file.path("..","..","..","data","processed","tusco","mmu","tusco_mouse.tsv"),
  file.path("..","data","tusco","tusco_mouse.tsv")
))

pipelines <- c(
  "ES_drna_ont",
  "ES_cdna_ont",
  "ES_cdna_pacbio",
  "ES_drna_ont_ls",
  "ES_cdna_ont_ls",
  "ES_cdna_pacbio_ls"
)

data_dir <- resolve_path(c(
  file.path("..","..","..","data","raw","lrgasp","mouse"),
  file.path("..","data","lrgasp","mouse")
), is_dir = TRUE)

###############################################################################
# 2) Read TSV safely
###############################################################################
read_tsv_safe <- function(file_path, col_names=TRUE, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  tryCatch({
    data <- read_tsv(file_path, col_names=col_names, ...)
    return(data)
  },
  error=function(e) {
    stop("Error reading file ", file_path, ": ", e$message)
  }
  )
}

# Lightweight GTF importer: uses rtracklayer if available, otherwise parses attributes
import_gtf_df <- function(gtf_path) {
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    return(as.data.frame(rtracklayer::import(gtf_path)))
  }
  cols <- readr::cols(
    seqnames = readr::col_character(),
    source   = readr::col_character(),
    type     = readr::col_character(),
    start    = readr::col_integer(),
    end      = readr::col_integer(),
    score    = readr::col_character(),
    strand   = readr::col_character(),
    frame    = readr::col_character(),
    attribute= readr::col_character()
  )
  df <- readr::read_tsv(
    gtf_path,
    comment = "#",
    col_names = c("seqnames","source","type","start","end","score","strand","frame","attribute"),
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

###############################################################################
# 3) Compute SIRVs: (TP, PTP, FP, FN)
###############################################################################
compute_sirv_bigcats <- function(classification_data, transcript_gtf_file, sirv_gtf_file) {

  # Read SIRV GTF
  sirv_gtf_df <- import_gtf_df(sirv_gtf_file)

  # If "transcript_id" is missing, try copying from "Parent"
  if (!"transcript_id" %in% names(sirv_gtf_df)) {
    if ("Parent" %in% names(sirv_gtf_df)) {
      sirv_gtf_df$transcript_id <- sirv_gtf_df$Parent
    } else {
      stop("No 'transcript_id' or 'Parent' attribute in SIRV GTF. Cannot continue.")
    }
  }

  # Distinct SIRV transcripts (from exons)
  annotation_data_sirv <- sirv_gtf_df %>%
    dplyr::filter(type=="exon") %>%
    dplyr::distinct(transcript_id) %>%
    dplyr::rename(ref_transcript_id=transcript_id)

  # Expand classification if "fusion"
  class_sirv <- classification_data %>%
    filter(structural_category!="fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category=="fusion") %>%
        separate_rows(associated_transcript, sep="_")
    ) %>%
    mutate(
      associated_gene=str_remove(associated_gene,"\\.\\d+$"),
      associated_transcript=str_remove(associated_transcript,"\\.\\d+$")
    ) %>%
    distinct(isoform, associated_transcript, .keep_all=TRUE)

  # Keep only SIRV chroms
  sirv_chroms <- unique(sirv_gtf_df$seqnames)
  class_sirv <- class_sirv %>% filter(chrom %in% sirv_chroms)

  SIRV_transcripts <- class_sirv %>% filter(grepl("SIRV", chrom))

  # Define categories
  TP_sirv <- SIRV_transcripts %>%
    filter(
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
  PTP_sirv <- SIRV_transcripts %>%
    filter((structural_category %in% c("full-splice_match","incomplete-splice_match") |
            subcategory == "mono-exon_by_intron_retention"),
           !(associated_transcript %in% TP_sirv$associated_transcript))

  ref_found <- union(TP_sirv$associated_transcript, PTP_sirv$associated_transcript)
  FN_sirv <- annotation_data_sirv %>%
    filter(!(ref_transcript_id %in% ref_found))

  FP_sirv <- SIRV_transcripts %>%
    filter(structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic",
                                      "fusion","antisense","intergenic","genic_intron") &
           subcategory != "mono-exon_by_intron_retention")

  summary_sirv <- data.frame(
    big_category=c("TP","PTP","FP","FN"),
    count=c(nrow(TP_sirv), nrow(PTP_sirv), nrow(FP_sirv), nrow(FN_sirv))
  )
  total_sirv <- sum(summary_sirv$count)
  summary_sirv$percentage <- 100 * summary_sirv$count / if (total_sirv==0) 1 else total_sirv

  return(summary_sirv)
}

###############################################################################
# 4) Compute TUSCO: (TP, PTP, FP, FN)
###############################################################################
compute_tusco_bigcats <- function(classification_data, tusco_tsv_file) {

  # Read TUSCO TSV
  # Robust read of header‑less Tusco annotation (falls back to whitespace if tabs not found)
  tusco_df <- read_delim(
    tusco_tsv_file,
    delim = "\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = cols(.default = "c"),
    trim_ws = TRUE
  )
  if (ncol(tusco_df) == 1) { # Fallback if not tab-delimited
    tusco_df <- read_table2(
      tusco_tsv_file,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = cols(.default = "c")
    )
  }

  if (!all(c("ensembl","refseq","gene_name") %in% colnames(tusco_df))) {
    stop("Tusco TSV must contain columns: ensembl, refseq, gene_name")
  }

  # Distinct annotation lines for genes
  annotation_data_tusco <- tusco_df %>%
    select(ensembl, refseq, gene_name) %>%
    distinct()

  # Patterns for ID types
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
    gene_name = "^[A-Z0-9]+$"
  )

  class_tusco <- classification_data %>%
    filter(structural_category!="fusion") %>%
    bind_rows(
      classification_data %>%
        filter(structural_category=="fusion") %>%
        separate_rows(associated_gene, sep="_")
    ) %>%
    mutate(
      associated_gene=str_remove(associated_gene,"\\.\\d+$")
    ) %>%
    mutate(
      id_type=case_when(
        str_detect(associated_gene, patterns$ensembl) ~ "ensembl",
        str_detect(associated_gene, patterns$refseq)  ~ "refseq",
        str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
        TRUE ~ "unknown"
      )
    ) %>%
    distinct(isoform, associated_gene, .keep_all=TRUE)

  # Determine top ID type used in classification data
  id_summary_tusco <- class_tusco %>%
    count(id_type, sort=TRUE)
    
  if (nrow(id_summary_tusco)==0 || (nrow(id_summary_tusco) > 0 && id_summary_tusco$id_type[1] == "unknown" && sum(id_summary_tusco$n) == id_summary_tusco$n[1])) {
     # Fallback or error if no dominant known ID type or only unknown IDs
    warning("No dominant known gene ID type (ensembl, refseq, gene_name) found in classifications for TUSCO. Results might be affected.")
    # Attempt to proceed or handle as error, here we'll try to use 'gene_name' as a default if available, or skip TUSCO processing for this sample
    if("gene_name" %in% id_summary_tusco$id_type[id_summary_tusco$id_type != "unknown"]) {
        top_type <- "gene_name"
    } else if ("ensembl" %in% id_summary_tusco$id_type[id_summary_tusco$id_type != "unknown"]) {
        top_type <- "ensembl"
    } else if ("refseq" %in% id_summary_tusco$id_type[id_summary_tusco$id_type != "unknown"]) {
        top_type <- "refseq"
    } else {
        # If no valid IDs, return empty summary
        return(data.frame(big_category=c("TP","PTP","FP","FN"), count=c(0,0,0,0), percentage=c(0,0,0,0)))
    }
  } else {
    # Filter out 'unknown' before selecting top_type if other types exist
    known_id_summary <- id_summary_tusco %>% filter(id_type != "unknown")
    if (nrow(known_id_summary) > 0) {
        top_type <- known_id_summary$id_type[1]
    } else { # Only unknown was present, but didn't meet the strict condition above, this is a fallback
        warning("Only 'unknown' gene ID types found for TUSCO. Returning empty summary.")
        return(data.frame(big_category=c("TP","PTP","FP","FN"), count=c(0,0,0,0), percentage=c(0,0,0,0)))
    }
  }
  
  cat("Using TUSCO ID type:", top_type, "for pipeline", 안전하게_파이프라인_접두사_가져오기(classification_data), "\n")


  # Keep only transcripts whose associated gene is in TUSCO annotation, using the determined top_type
  TUSCO_transcripts <- class_tusco %>%
    filter(id_type == top_type) %>%
    filter(associated_gene %in% annotation_data_tusco[[top_type]])

  # Define categories
  TP_tusco <- TUSCO_transcripts %>%
    filter(
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
  PTP_tusco <- TUSCO_transcripts %>%
    filter((structural_category %in% c("full-splice_match","incomplete-splice_match") |
            subcategory == "mono-exon_by_intron_retention"),
           !(isoform %in% TP_tusco$isoform))

  # For FN, consider genes from TUSCO annotation not found in TP or PTP
  ref_ids_tusco <- annotation_data_tusco[[top_type]]
  found_ids_tusco <- union(TP_tusco$associated_gene, PTP_tusco$associated_gene)
  
  FN_tusco_df <- annotation_data_tusco %>%
    filter(!(!!sym(top_type) %in% found_ids_tusco))
  FN_tusco_count <- nrow(FN_tusco_df)

  FP_tusco <- TUSCO_transcripts %>%
    filter(structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic",
                                      "fusion","antisense","intergenic","genic_intron") &
           subcategory != "mono-exon_by_intron_retention")

  summary_tusco <- data.frame(
    big_category=c("TP","PTP","FP","FN"),
    count=c(nrow(TP_tusco), nrow(PTP_tusco), nrow(FP_tusco), FN_tusco_count)
  )
  total_tusco <- sum(summary_tusco$count)
  summary_tusco$percentage <- 100 * summary_tusco$count / if (total_tusco==0) 1 else total_tusco

  return(summary_tusco)
}

# Helper function to safely get pipeline prefix from classification_data for logging
안전하게_파이프라인_접두사_가져오기 <- function(classification_data) {
  if ("pipeline" %in% colnames(classification_data) && nrow(classification_data) > 0) {
    return(as.character(classification_data$pipeline[1]))
  } else if (exists("pipeline_prefix", envir = parent.frame())) {
    return(get("pipeline_prefix", envir = parent.frame()))
  }
  return("Unknown Pipeline")
}

###############################################################################
# 5) Process a single pipeline => returns pipeline-level results
###############################################################################
process_pipeline <- function(pipeline_prefix, data_dir, sirv_gtf, tusco_tsv) {
  class_file <- file.path(data_dir, paste0(pipeline_prefix, "/", pipeline_prefix, "_classification.txt"))
  gtf_file   <- file.path(data_dir, paste0(pipeline_prefix, "/", pipeline_prefix, "_corrected.gtf"))

  classification_data <- read_tsv_safe(class_file)
  if (nrow(classification_data) > 0) {
    classification_data$pipeline <- pipeline_prefix
  }

  # SIRVs
  sirv_res <- compute_sirv_bigcats(classification_data, gtf_file, sirv_gtf) %>%
    mutate(Type="SIRVs", pipeline=pipeline_prefix)

  # TUSCO
  tusco_res <- compute_tusco_bigcats(classification_data, tusco_tsv) %>%
    mutate(Type="TUSCO", pipeline=pipeline_prefix)

  bind_rows(sirv_res, tusco_res)
}

###############################################################################
# 6) Main
###############################################################################
cat("\n--- Gathering data for all pipelines ---\n")

all_data <- lapply(pipelines, function(p) {
                    res <- process_pipeline(p, 
                                            data_dir = data_dir, 
                                            sirv_gtf = sirv_gtf_file, 
                                            tusco_tsv = tusco_tsv_file)
                    res$pipeline <- p
                    return(res)
                  }
                  ) %>%
  bind_rows()

cat("\n--- Single test (one-sided paired t-test: TUSCO > SIRVs) for each category ---\n")
categories <- c("TP","PTP","FP","FN")

p_values <- c()
for (cat_name in categories) {
  cat("\nCategory:", cat_name, "\n")
  df_sub <- all_data %>% filter(big_category==cat_name)
  wide <- df_sub %>%
    select(pipeline, Type, percentage) %>%
    pivot_wider(names_from=Type, values_from=percentage)
  
  if ("TUSCO" %in% colnames(wide) && "SIRVs" %in% colnames(wide) && 
      sum(!is.na(wide$TUSCO)) > 1 && sum(!is.na(wide$SIRVs)) > 1) {
    # One-sided paired t-test: TUSCO > SIRVs
    test_out <- t.test(wide$TUSCO, wide$SIRVs, paired=TRUE, alternative="greater")
    pval <- test_out$p.value
    cat("p-value = ", pval, "\n")
  } else {
    cat("Not enough data or TUSCO/SIRVs columns missing for category:", cat_name, ". Skipping t-test.\n")
  }
  p_values <- c(p_values, pval)
}

# Build data to show significance label
# We'll do (TUSCO vs. SIRVs) for each category
my_signifs <- data.frame(
  big_category = categories,
  p_value = p_values
)

# Convert p-values -> stars
p_stars <- function(x) {
  if (x < 0.001) return("***")
  else if (x < 0.01) return("**")
  else if (x < 0.05) return("*")
  else return("ns")
}
my_signifs$star_label <- sapply(my_signifs$p_value, p_stars)

cat("\n--- Summaries: ---\n")
print(my_signifs)

# Compute mean ± SD for each category × Type
mean_data <- all_data %>%
  group_by(big_category, Type) %>%
  summarise(
    mean_perc = mean(percentage, na.rm = TRUE),
    sd_perc   = sd(percentage, na.rm = TRUE),
    .groups   = "drop"
  )

# Ensure proper ordering on the x-axis
mean_data$big_category <- factor(mean_data$big_category, levels = c("TP","PTP","FP","FN"))
mean_data$Type         <- factor(mean_data$Type, levels = c("TUSCO","SIRVs"))

# Create significance annotation data
bracket_data <- my_signifs %>%
  mutate(
    cat_index     = as.numeric(factor(big_category, levels = c("TP","PTP","FP","FN"))),
    y_position    = 110, # Adjust this based on your data range
    p_value_stars = star_label
  )

# Ensure correct ordering of categories
all_data$big_category <- factor(all_data$big_category, levels = c("TP", "PTP", "FP", "FN"))
bracket_data$big_category <- factor(bracket_data$big_category, levels = c("TP", "PTP", "FP", "FN"))

# Adjust positions for significance annotations and horizontal lines
bracket_data <- bracket_data %>%
  mutate(
    y_position = 110,         # Adjust for significance text
    y_line = 105,             # Position for horizontal lines
    xmin = as.numeric(big_category) - 0.2,
    xmax = as.numeric(big_category) + 0.2
  )

# Generate the barplot
p_single <- ggplot(all_data, aes(x = big_category, y = percentage, fill = Type)) +
  
  # 1. Bars for mean with reduced thickness and narrower width
  stat_summary(
    fun    = mean,
    geom   = "bar",
    width  = 0.4,             # Reduced width for narrower bars
    color  = "black",
    size   = 0.1,             # Reduced border thickness
    position = position_dodge(width = 0.6)
  ) +
  
  # 2. Error bars for SD with minimum thickness
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom     = "errorbar",
    width    = 0.2,
    size     = 0.1,             # Further reduced error bar thickness
    position = position_dodge(width = 0.6)
  ) +
  
  # 3. Jittered raw data points with further reduced size
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size     = 0.2,             # Further reduced sample point size
    alpha    = 0.7,
    shape    = 16,
    color    = "black"
  ) +
  
  # 4. Significance labels with reduced size
  geom_text(
    data       = bracket_data,
    aes(
      x          = big_category,
      y          = y_position,
      label      = p_value_stars
    ),
    inherit.aes = FALSE,
    size        = 3,             # Reduced size for smaller font
    fontface    = "bold"
  ) +
  
  # 5. Horizontal significance lines with further reduced thickness
  geom_segment(
    data       = bracket_data,
    aes(
      x    = xmin,
      xend = xmax,
      y    = y_line,
      yend = y_line
    ),
    inherit.aes = FALSE,
    size        = 0.15,          # Further reduced line thickness
    color       = "black"
  ) +
  
  # 6. Adjust colors and theme
  scale_fill_manual(values = c("TUSCO" = "#1b9e77", "SIRVs" = "#cab2d6")) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +  # Set y-axis max to 100
  labs(
    x     = NULL,
    y     = "Percentage"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 7),
    axis.title.y = element_text(size = 7, face = "bold"),
    legend.position = "none",       # Remove the legend
    legend.title = element_blank(),
    plot.margin  = margin(15, 10, 10, 10)  # Adjusted top margin
  )

# Do not display the plot; rely on PDF writer below
# print(p_single)

## Helper to resolve preferred paths: try absolute figs/data, then repo-relative figs/data
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
  }
  return(candidates[[1]])
}

## Outputs: write only under this figure folder
plot_dir <- base::file.path("..", "plots")
tsv_dir  <- base::file.path("..", "tables")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(tsv_dir))  dir.create(tsv_dir,  recursive = TRUE)

pdf_path <- file.path(plot_dir, "figure3b-mouse.pdf")
tsv_path <- file.path(tsv_dir,  "figure3b-mouse.tsv")

# Save PDF using base device to avoid Cairo/OpenMP issues
pdf(file = pdf_path, width = 2.5, height = 1.8)
print(p_single)
dev.off()

# Prepare and write TSV capturing underlying data and minimal metadata
group_counts <- all_data %>% count(big_category, Type, name = "n")
raw_block <- all_data %>%
  left_join(group_counts, by = c("big_category", "Type")) %>%
  mutate(record_type = "raw")

summary_block <- mean_data %>%
  mutate(record_type = "summary")

stat_block <- bracket_data %>%
  select(big_category, y_position, y_line, p_value_stars) %>%
  mutate(record_type = "stat")

tsv_out <- bind_rows(
  raw_block %>% mutate(figure_id = "fig-3", panel_id = "3b-mouse"),
  summary_block %>% mutate(figure_id = "fig-3", panel_id = "3b-mouse"),
  stat_block %>% mutate(figure_id = "fig-3", panel_id = "3b-mouse")
)
readr::write_tsv(tsv_out, tsv_path)

# === EXTRACTION FOR PAPER VALUES ===
# This section extracts key values for manuscript updates
extraction_file <- file.path("..", "..", "..", "manuscript", "extraction_values.tsv")

# Get p-values for PTP and FN (indices 2 and 4 in categories vector)
ptp_pval <- p_values[2]  # PTP p-value
fn_pval <- p_values[4]   # FN p-value

# Prepare extraction data
extraction_data <- data.frame(
  figure_id = "fig-3b",
  metric_name = c("mouse_PTP_pvalue", "mouse_FN_pvalue"),
  value = c(formatC(ptp_pval, format="e", digits=3), formatC(fn_pval, format="e", digits=3)),
  ci_lower = NA,
  ci_upper = NA,
  notes = "One-sided paired t-test (TUSCO > SIRVs)",
  stringsAsFactors = FALSE
)

# Write or append to file
if (!file.exists(extraction_file)) {
  write.table(extraction_data, extraction_file, sep="\t", 
              row.names=FALSE, quote=FALSE)
  cat("Created extraction file:", extraction_file, "\n")
} else {
  write.table(extraction_data, extraction_file, sep="\t", 
              row.names=FALSE, quote=FALSE, append=TRUE, col.names=FALSE)
}
cat("Extracted p-values: PTP =", formatC(ptp_pval, format="e", digits=3), 
    ", FN =", formatC(fn_pval, format="e", digits=3), "\n")
# === END EXTRACTION ===

cat("Plot saved successfully:\n")
cat(" - PDF:", pdf_path, "\n")
cat(" - TSV:", tsv_path, "\n")
