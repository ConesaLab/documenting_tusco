# -----------------------------
# Full Code: TUSCO Metrics & Plotting for RIN Correlation
# -----------------------------

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)  # for get_legend()
})

# Mitigate OpenMP SHM issues in restricted environments
Sys.setenv(OMP_NUM_THREADS = "1", OMP_PROC_BIND = "FALSE", OMP_WAIT_POLICY = "PASSIVE", KMP_INIT_AT_FORK = "0")

# Path helpers
# - Detect script/figure/repo roots robustly
# - Resolve a list of candidate paths (first that exists wins)
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
  }
  return(candidates[[1]])
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

# Detect paths relative to this script so it works from any cwd
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0 || !nzchar(script_path)) {
  # Fallback for interactive runs
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
} else {
  script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
}
fig_dir  <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
repo_root <- normalizePath(file.path(fig_dir, "..", ".."), winslash = "/", mustWork = FALSE)

# (Optional) Debugging options
# options(error = recover)
# traceback()

# -----------------------------------------------------------------------------
# 1. Helper Functions and TUSCO Metrics Calculation Functions
# -----------------------------------------------------------------------------

# A helper function to read TSV files safely (if needed)
read_tsv_safe <- function(file_path, col_names = TRUE, ...) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }

  tryCatch({
    data <- read_tsv(file_path, col_names = col_names, ...)
    message("Successfully read file: ", file_path)
    return(data)
  }, error = function(e) {
    warning("Error reading file ", file_path, ": ", e$message)
    return(NULL)
  })
}

# Function to calculate TUSCO metrics (TP, PTP, FP, FN) for a single classification file
calculate_tusco_metrics <- function(classification_data,
                                    tusco_annot_file) {
  # 1) Read and prepare the TUSCO TSV annotation
  tusco_df <- read_delim(
    tusco_annot_file,
    delim = "\t",
    col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
    col_types = cols(.default = "c"),
    trim_ws = TRUE,
    comment = "#" # Assuming comments might start with #
  )
  if (ncol(tusco_df) < 3) { # Check if parsing failed and returned fewer than expected core columns
    tusco_df <- read_table2( # Fallback to read_table2 if tab delimiting fails
      tusco_annot_file,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = cols(.default = "c"),
      comment = "#"
    )
  }

  if (!all(c("ensembl","refseq","gene_name") %in% colnames(tusco_df))) {
    stop("Tusco TSV must contain columns: ensembl, refseq, gene_name. Found: ", paste(colnames(tusco_df), collapse=", "))
  }

  annotation_data_tusco <- tusco_df %>%
    dplyr::select(ensembl, refseq, gene_name) %>%
    dplyr::distinct()

  # Total number of TUSCO reference genes (for sensitivity calculations)
  rTUSCO <- nrow(annotation_data_tusco)

  # Define regex patterns for gene ID types
  patterns <- list(
    ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
    refseq    = "^(NM_|NR_|NP_)\\d{6,}(\\.\\d+)?$", # Adjusted RefSeq to allow versions
    gene_name = "^[A-Za-z0-9][A-Za-z0-9.-]*[A-Za-z0-9]$" # Adjusted for more complex gene names
  )

  # 2) Clean up classification data for TUSCO
  classification_data_cleaned_tusco <- classification_data %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(
      classification_data %>%
        dplyr::filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_gene, sep = "_")
    ) %>%
    dplyr::mutate(
      associated_gene = stringr::str_remove(associated_gene, "\\.\\d+$")
    ) %>%
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

  id_summary_tusco <- classification_data_cleaned_tusco %>%
    dplyr::count(id_type, sort = TRUE)

  if (nrow(id_summary_tusco) == 0 || all(id_summary_tusco$id_type == "unknown")) {
     warning("No gene IDs matched the known patterns for TUSCO or all were unknown. Check patterns and data.")
    # Create a placeholder for top_id_type_tusco if no known IDs are found to avoid error later
    # This part might need more robust handling based on expected data
    if (any(c("ensembl", "refseq", "gene_name") %in% colnames(annotation_data_tusco))) {
        top_id_type_tusco <- intersect(c("ensembl", "refseq", "gene_name"), colnames(annotation_data_tusco))[1]
        warning(paste("Defaulting TUSCO ID type to:", top_id_type_tusco))
    } else {
        stop("Cannot determine a valid ID type for TUSCO annotation matching.")
    }
  } else {
     top_id_tusco <- id_summary_tusco %>% dplyr::filter(id_type != "unknown") %>% dplyr::slice_max(n, n = 1)
     if (nrow(top_id_tusco) == 0) { # Handles case where only 'unknown' IDs were found
         top_id_tusco <- id_summary_tusco %>% dplyr::slice_max(n, n = 1) # Fallback to most common, even if unknown
         warning("Only 'unknown' ID types found in classification data. Proceeding with the most frequent one.")
     }
     top_id_type_tusco <- top_id_tusco$id_type[1] # Take the first if multiple have max count
  }
  
  message(paste("Top ID type for TUSCO matching:", top_id_type_tusco))


  # Keep only transcripts that match the TUSCO annotation in the top ID type
  # Ensure the column exists in annotation_data_tusco
  if (!top_id_type_tusco %in% colnames(annotation_data_tusco)) {
    stop(paste("The determined top_id_type_tusco '", top_id_type_tusco, "' is not a column in annotation_data_tusco.",
               "Available columns: ", paste(colnames(annotation_data_tusco), collapse=", ")))
  }
  
  TUSCO_transcripts <- classification_data_cleaned_tusco %>%
    dplyr::filter(id_type == top_id_type_tusco) %>% # Filter by detected ID type first
    dplyr::filter(associated_gene %in% annotation_data_tusco[[top_id_type_tusco]])

  # 3) Define True Positives (TP), Partially True Positives (PTP),
  # False Negatives (FN), and False Positives (FP)
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
  TP_tusco <- TUSCO_RM
  TP_TUSCO <- unique(TP_tusco$associated_gene)

  PTP_tusco <- TUSCO_transcripts %>%
    dplyr::filter((structural_category %in% c("full-splice_match", "incomplete-splice_match") |
                     subcategory == "mono-exon_by_intron_retention") &
                    !(isoform %in% TP_tusco$isoform))

  FN_tusco <- annotation_data_tusco %>%
    dplyr::filter(!( .data[[top_id_type_tusco]] %in% TUSCO_transcripts$associated_gene ))

  FP_tusco <- TUSCO_transcripts %>%
    dplyr::filter(structural_category %in% c(
      "novel_in_catalog", "novel_not_in_catalog", "genic",
      "fusion", "antisense", "intergenic", "genic_intron"
    ) &
    subcategory != "mono-exon_by_intron_retention")

  # 4) Compute raw counts
  TP_count  <- nrow(TP_tusco)
  PTP_count <- nrow(PTP_tusco)
  FP_count  <- nrow(FP_tusco)
  FN_count  <- nrow(FN_tusco)

  # 5) Return a list of the raw counts
  metrics_list <- list(
    TP = TP_count,
    PTP = PTP_count,
    FP = FP_count,
    FN = FN_count
  )

  return(metrics_list)
}

# Function to calculate MANE metrics (TP, PTP, FP, FN) for a single classification file
calculate_mane_metrics <- function(classification_data, mane_gtf_file) {
  # 1) Read and prepare the MANE GTF annotation
  if (!file.exists(mane_gtf_file)) {
    stop("MANE GTF file not found: ", mane_gtf_file)
  }
  mane_gtf_df <- tryCatch({
    import_gtf_df(mane_gtf_file)
  }, error = function(e) {
    stop("Error importing MANE GTF file ", mane_gtf_file, ": ", e$message)
  })

  if (!"transcript_id" %in% names(mane_gtf_df)) {
     warning("MANE GTF does not have a 'transcript_id' column. Attempting to use 'transcript'.")
     if ("transcript" %in% names(mane_gtf_df)) {
        mane_gtf_df$transcript_id <- mane_gtf_df$transcript
     } else {
        stop("MANE GTF must contain a 'transcript_id' or 'transcript' column.")
     }
  }
  
  annotation_data_mane_intermediate <- mane_gtf_df %>%
    dplyr::filter(type == "transcript") %>% 
    dplyr::select(transcript_id) %>%
    dplyr::filter(!is.na(transcript_id)) %>%
    dplyr::distinct() 

  message("---- MANE Debug Info for Sample ----")
  message("MANE Debug: Sample of transcript_id BEFORE version removal (first 5 if available):")
  print(head(annotation_data_mane_intermediate$transcript_id, 5))
  message("MANE Debug: str_detect for '\\\\.\\\\d+$' on these (first 5 if available):")
  print(head(stringr::str_detect(annotation_data_mane_intermediate$transcript_id, "\\.\\d+$"), 5))

  annotation_data_mane <- annotation_data_mane_intermediate %>%
    dplyr::mutate(
      transcript_id_trimmed = trimws(transcript_id), # Add trimws
      transcript_id_no_version = stringr::str_remove(transcript_id_trimmed, "\\.\\d+$")
    )
  
  message("MANE Debug: Sample of transcript_id_trimmed (first 5 if available):")
  print(head(annotation_data_mane$transcript_id_trimmed, 5))
  # The existing debug messages for transcript_id_no_version will now reflect the result of the above mutate

  if (nrow(annotation_data_mane) == 0) {
    stop("No transcript IDs found in the MANE GTF file.")
  }

  # 2) Clean up classification data for MANE matching
  classification_data_cleaned_mane <- classification_data %>%
    dplyr::filter(structural_category != "fusion") %>%
    dplyr::bind_rows(
      classification_data %>%
        dplyr::filter(structural_category == "fusion") %>%
        tidyr::separate_rows(associated_transcript, sep = "_")
    ) %>%
    dplyr::mutate(
      associated_transcript_no_version = stringr::str_remove(associated_transcript, "\\\\.\\\\d+$")
    ) %>%
    dplyr::distinct(isoform, associated_transcript, .keep_all = TRUE)

  MANE_transcripts_matched <- classification_data_cleaned_mane %>%
    dplyr::filter(associated_transcript_no_version %in% annotation_data_mane$transcript_id_no_version)

  # 3) Define TP, PTP, FN, FP (isoform-based counts for TP, PTP, FP; reference-based for FN)
  TP_isoforms_df <- MANE_transcripts_matched %>%
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
    )
  TP_count_isoform_based <- nrow(TP_isoforms_df)

  PTP_isoforms_df <- MANE_transcripts_matched %>%
    dplyr::filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
    dplyr::anti_join(TP_isoforms_df, by="isoform")
  PTP_count_isoform_based <- nrow(PTP_isoforms_df)
  
  all_mane_ref_ids_no_version <- unique(annotation_data_mane$transcript_id_no_version)
  mane_ids_hit_by_tp_isoforms <- unique(TP_isoforms_df$associated_transcript_no_version)
  mane_ids_hit_by_ptp_isoforms <- unique(PTP_isoforms_df$associated_transcript_no_version)
  all_found_mane_ref_ids_by_isoforms <- dplyr::union(mane_ids_hit_by_tp_isoforms, mane_ids_hit_by_ptp_isoforms)
  FN_count_revised <- length(base::setdiff(all_mane_ref_ids_no_version, all_found_mane_ref_ids_by_isoforms))

  FP_isoforms_df_consistent <- MANE_transcripts_matched %>%
    dplyr::filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog",
                                               "genic", "fusion", "antisense", "intergenic", "genic_intron"))
  FP_count_isoform_based <- nrow(FP_isoforms_df_consistent)

  metrics_list <- list(
    TP = TP_count_isoform_based,
    PTP = PTP_count_isoform_based,
    FP = FP_count_isoform_based,
    FN = FN_count_revised
  )

  # Debugging messages
  message("Sample of MANE GTF transcript_id_no_version (first 5 if available):")
  print(head(annotation_data_mane$transcript_id_no_version, 5))
  message("Number of unique MANE GTF transcript_id_no_version: ", length(unique(annotation_data_mane$transcript_id_no_version)))

  message("Sample of classification associated_transcript_no_version (first 5 if available):")
  print(head(classification_data_cleaned_mane$associated_transcript_no_version, 5))
  message("Number of unique classification associated_transcript_no_version: ", length(unique(classification_data_cleaned_mane$associated_transcript_no_version)))
  
  MANE_transcripts_matched <- classification_data_cleaned_mane %>%
    dplyr::filter(associated_transcript_no_version %in% annotation_data_mane$transcript_id_no_version)

  message("Number of rows in MANE_transcripts_matched: ", nrow(MANE_transcripts_matched))
  if (nrow(MANE_transcripts_matched) > 0) {
    message("Sample of matched associated_transcript_no_version (first 5 if available):")
    print(head(MANE_transcripts_matched$associated_transcript_no_version, 5))
  }
  message("------------------------------------")

  return(metrics_list)
}

# Main function: Process all TUSCO classification files in subdirectories
process_all_classifications_tusco <- function(main_data_dir, tusco_annot_file, output_file) {
  classification_dirs <- list.dirs(path = main_data_dir, full.names = TRUE, recursive = FALSE)

  all_metrics <- data.frame(
    Sample = character(),
    TP = integer(),
    PTP = integer(),
    FP = integer(),
    FN = integer(),
    stringsAsFactors = FALSE
  )

  for (dir in classification_dirs) {
    sample_name <- basename(dir)
    message("\nProcessing sample: ", sample_name)

    classification_file <- file.path(dir, paste0(sample_name, "_classification.txt"))

    if (!file.exists(classification_file)) {
      warning("Classification file not found for sample: ", sample_name,
              "\nExpected at: ", classification_file, "\nSkipping this sample.")
      next
    }

    classification_data <- read_tsv_safe(classification_file) # Use read_tsv_safe
    if (is.null(classification_data)) {
        warning("Failed to read classification data for sample: ", sample_name, "\nSkipping this sample.")
        next
    }


    metrics <- calculate_tusco_metrics(classification_data, tusco_annot_file)
    if (is.null(metrics)) {
      warning("Metrics calculation failed for sample: ", sample_name, "\nSkipping this sample.")
      next
    }

    all_metrics <- all_metrics %>%
      add_row(
        Sample = sample_name,
        TP = metrics$TP,
        PTP = metrics$PTP,
        FP = metrics$FP,
        FN = metrics$FN
      )
  }

  # Suppress per-metric file writes to avoid extra outputs under ./plots
  # write_csv(all_metrics, output_file)
  message("\nTUSCO metrics computed for in-memory use.")
  
  return(all_metrics)
}

# Main function: Process all MANE classification files in subdirectories
process_all_classifications_mane <- function(main_data_dir, mane_gtf_file, output_file_mane) {
  classification_dirs <- list.dirs(path = main_data_dir, full.names = TRUE, recursive = FALSE)

  all_metrics_mane <- data.frame(
    Sample = character(),
    TP = integer(),
    PTP = integer(),
    FP = integer(),
    FN = integer(),
    stringsAsFactors = FALSE
  )

  for (dir in classification_dirs) {
    sample_name <- basename(dir)
    message("\\nProcessing sample for MANE: ", sample_name)

    classification_file <- file.path(dir, paste0(sample_name, "_classification.txt"))

    if (!file.exists(classification_file)) {
      warning("Classification file not found for sample (MANE): ", sample_name,
              "\\nExpected at: ", classification_file, "\\nSkipping this sample.")
      next
    }

    classification_data <- read_tsv_safe(classification_file)
    if (is.null(classification_data)) {
        warning("Failed to read classification data for sample (MANE): ", sample_name, "\\nSkipping this sample.")
        next
    }

    metrics_mane <- calculate_mane_metrics(classification_data, mane_gtf_file)
    if (is.null(metrics_mane)) {
      warning("MANE Metrics calculation failed for sample: ", sample_name, "\\nSkipping this sample.")
      next
    }

    all_metrics_mane <- all_metrics_mane %>%
      add_row(
        Sample = sample_name,
        TP = metrics_mane$TP,
        PTP = metrics_mane$PTP,
        FP = metrics_mane$FP,
        FN = metrics_mane$FN
      )
  }

  # Suppress per-metric file writes to avoid extra outputs under ./plots
  # write_csv(all_metrics_mane, output_file_mane)
  message("\\nMANE metrics computed for in-memory use.")

  return(all_metrics_mane)
}

# -----------------------------------------------------------------------------
# 2. Process TUSCO Classification Files
# -----------------------------------------------------------------------------

# Define file paths via repo-rooted data layout (with fallbacks)
# RIN classification parent directory (contains per-sample subfolders)
main_data_dir <- resolve_path(c(
  file.path(repo_root, "data", "raw", "RIN"),
  file.path(fig_dir,    "data", "RIN")
), is_dir = TRUE)

# TUSCO annotation (human)
tusco_annot_file <- resolve_path(c(
  file.path(repo_root, "data", "processed", "tusco", "hsa", "tusco_human.tsv"),
  file.path(repo_root, "data", "raw",       "tusco", "hsa", "tusco_human.tsv"),
  file.path(fig_dir,   "data", "tusco", "tusco_human.tsv")
))

output_file <- NULL

# MANE GTF (prefer compressed under data/raw)
mane_gtf_file <- resolve_path(c(
  file.path(repo_root, "data", "raw", "reference", "human", "mane.v1.4.ensembl_genomic.gtf.gz"),
  file.path(repo_root, "data", "raw", "reference", "human", "mane.v1.4.ensembl_genomic.gtf"),
  file.path(fig_dir,    "data", "reference", "human", "mane.v1.4.ensembl_genomic.gtf.gz"),
  file.path(fig_dir,    "data", "reference", "human", "mane.v1.4.ensembl_genomic.gtf")
))
output_file_mane <- NULL

# Create the main_data_dir if it doesn't exist, with a note to the user
if (!dir.exists(main_data_dir)) {
  # dir.create(main_data_dir, recursive = TRUE, showWarnings = FALSE)
  warning(paste("The directory for TUSCO classification files (main_data_dir):\n",
              main_data_dir, # Use the variable here
              "\ndoes not exist. Please create it and populate it with your data, or update the 'main_data_dir' variable in this script."))
}


# Process the TUSCO classification files and get the metrics summary data frame
tusco_metrics_summary <- process_all_classifications_tusco(
  main_data_dir = main_data_dir,
  tusco_annot_file = tusco_annot_file,
  output_file = output_file
)

# Process the MANE classification files
mane_metrics_summary <- process_all_classifications_mane(
  main_data_dir = main_data_dir,
  mane_gtf_file = mane_gtf_file,
  output_file = output_file_mane
)

# The TUSCO metrics summary (from your classification files) has columns:
#   Sample, TP, PTP, FP, FN
# We now extract the RIN value and sample Group from the Sample names.
tusco_data <- tusco_metrics_summary %>%
  mutate(
    RIN = as.numeric(str_extract(Sample, "(?<=RIN_)[0-9.]+")),
    Group = str_extract(Sample, "^(TS[0-9]+)"), # This pattern might need adjustment for TUSCO sample names
    Dataset = "TUSCO"
  ) %>%
  select(RIN, TP, PTP, FP, FN, Group, Dataset)

message("TUSCO data (from classification metrics):")
print(tusco_data)

# Prepare MANE data
mane_data <- mane_metrics_summary %>%
  mutate(
    RIN = as.numeric(str_extract(Sample, "(?<=RIN_)[0-9.]+")),
    Group = str_extract(Sample, "^(TS[0-9]+)"),
    Dataset = "MANE"
  ) %>%
  select(RIN, TP, PTP, FP, FN, Group, Dataset)

message("MANE data (from classification metrics):")
print(mane_data)

# -----------------------------------------------------------------------------
# 3. Define Hard-Coded Sequin (SIRV) Data (Kept as is from original)
# -----------------------------------------------------------------------------

# TS10 Sequin data
TS10_sequin_RIN <- c(7.7, 7.7, 8.4, 9.3, 9.6, 9.8)
TS10_sequin_TP  <- c(22, 45, 32, 34, 31, 29)
TS10_sequin_PTP <- c(1, 4, 1, 3, 2, 0)
TS10_sequin_FP  <- c(2, 4, 4, 4, 4, 2)
TS10_sequin_FN  <- c(137, 111, 127, 123, 127, 131)

TS10_sequin <- data.frame(
  RIN = TS10_sequin_RIN,
  TP = TS10_sequin_TP,
  PTP = TS10_sequin_PTP,
  FP = TS10_sequin_FP,
  FN = TS10_sequin_FN,
  Group = "TS10",
  Dataset = "Sequin"
)

# TS11 Sequin data
TS11_sequin_RIN <- c(8.7, 8.8, 8.9, 9.3, 9.6, 9.7, 9.7)
TS11_sequin_TP  <- c(34, 38, 30, 31, 26, 32, 40)
TS11_sequin_PTP <- c(6, 5, 3, 1, 1, 2, 7)
TS11_sequin_FP  <- c(0, 4, 1, 3, 0, 3, 6)
TS11_sequin_FN  <- c(120, 117, 127, 128, 133, 126, 113)

TS11_sequin <- data.frame(
  RIN = TS11_sequin_RIN,
  TP = TS11_sequin_TP,
  PTP = TS11_sequin_PTP,
  FP = TS11_sequin_FP,
  FN = TS11_sequin_FN,
  Group = "TS11",
  Dataset = "Sequin"
)

# TS12 Sequin data
TS12_sequin_RIN <- c(7.2, 7.3, 8.2, 8.2, 9.9)
TS12_sequin_TP  <- c(27, 27, 32, 38, 44)
TS12_sequin_PTP <- c(3, 2, 4, 7, 8)
TS12_sequin_FP  <- c(1, 0, 2, 6, 6)
TS12_sequin_FN  <- c(130, 131, 124, 115, 108)

TS12_sequin <- data.frame(
  RIN = TS12_sequin_RIN,
  TP = TS12_sequin_TP,
  PTP = TS12_sequin_PTP,
  FP = TS12_sequin_FP,
  FN = TS12_sequin_FN,
  Group = "TS12",
  Dataset = "Sequin"
)

sequin_data <- rbind(TS10_sequin, TS11_sequin, TS12_sequin)
message("Sequin data (hard-coded):")
print(sequin_data)

# -----------------------------------------------------------------------------
# 4. Combine TUSCO and Sequin Data and Compute Derived Metric
# -----------------------------------------------------------------------------

all_data <- bind_rows(tusco_data, sequin_data, mane_data)

calculate_derived_metrics <- function(df) {
  df <- df %>%
    mutate(TP_TPPTP = ifelse((TP + PTP) > 0, TP / (TP + PTP), NA)) # TP / (TP + PTP)
  return(df)
}

all_data <- calculate_derived_metrics(all_data)

message("Combined data for plotting:")
print(all_data)

# -----------------------------------------------------------------------------
# 5. Plotting: Create Correlation Plots and Save the Final Figure
# -----------------------------------------------------------------------------

group_shapes <- c("TS10" = 16, "TS11" = 17, "TS12" = 18) # Keep or adjust as needed
dataset_colors <- c("TUSCO" = "#a8d5a0", "Sequin" = "#bb80b1", "MANE" = "#e41a1c") # Updated TUSCO color, Added MANE color

custom_theme <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.text   = element_text(size = 7, color = "black"),
    axis.title  = element_text(size = 7),
    axis.line   = element_line(color = "black"),
    plot.title  = element_text(size = 8, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text  = element_text(size = 7)
  )

# Create correlation plot for TUSCO data
p_tusco <- ggplot(subset(all_data, Dataset == "TUSCO"), aes(x = TP_TPPTP, y = RIN)) +
  geom_point(aes(color = Dataset, shape = Group), size = 1, alpha = 0.9) +
  geom_smooth(method = "lm", color = "darkgray", fill = "gray",
              linetype = "solid", alpha = 0.2, se = TRUE, size = 0.5) +
  stat_cor(method = "pearson",
           # cor.coef.name = "R^2", # Use R for Pearson
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), # Show R and p-value
           parse = FALSE, # Changed to FALSE to allow custom label with p.label
           label.x.npc = "left", label.y.npc = "top", # Position R and p-value
           size = 2.5) +
  labs(x = expression("TP / (TP + PTP)"), y = "RIN") +
  scale_color_manual(values = dataset_colors) +
  scale_shape_manual(values = group_shapes) +
  custom_theme +
  theme(legend.position = "none") +
  ggtitle("TUSCO")


# Create correlation plot for Sequin data
p_sequin <- ggplot(subset(all_data, Dataset == "Sequin"), aes(x = TP_TPPTP, y = RIN)) +
  geom_point(aes(color = Dataset, shape = Group), size = 1, alpha = 0.9) +
  geom_smooth(method = "lm", color = "darkgray", fill = "gray",
              linetype = "solid", alpha = 0.2, se = TRUE, size = 0.5) +
  stat_cor(method = "pearson",
           # cor.coef.name = "R^2",
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           parse = FALSE,
           label.x.npc = "left", label.y.npc = "top", # Position R and p-value
           size = 2.5) +
  scale_x_continuous(limits = c(NA, NA)) + # Adjust limits if necessary, e.g. c(0.8,1) from original if data range is similar
  labs(x = expression("TP / (TP + PTP)"), y = "RIN") +
  scale_color_manual(values = dataset_colors) +
  scale_shape_manual(values = group_shapes) +
  custom_theme +
  theme(legend.position = "none") +
  ggtitle("Sequin")

# Create correlation plot for MANE data
p_mane <- ggplot(subset(all_data, Dataset == "MANE"), aes(x = TP_TPPTP, y = RIN)) +
  geom_point(aes(color = Dataset, shape = Group), size = 1, alpha = 0.9) +
  geom_smooth(method = "lm", color = "darkgray", fill = "gray",
              linetype = "solid", alpha = 0.2, se = TRUE, size = 0.5) +
  stat_cor(method = "pearson",
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           parse = FALSE,
           label.x.npc = "left", label.y.npc = "top",
           size = 2.5) +
  labs(x = expression("TP / (TP + PTP)"), y = "RIN") +
  scale_color_manual(values = dataset_colors) +
  scale_shape_manual(values = group_shapes) +
  custom_theme +
  theme(legend.position = "none") +
  ggtitle("MANE")


legend_plot <- ggplot(all_data, aes(x = TP_TPPTP, y = RIN, color = Dataset, shape = Group)) +
  geom_point() +
  scale_color_manual(values = dataset_colors, name = "Dataset") + # Added legend title
  scale_shape_manual(values = group_shapes, name = "Group") +   # Added legend title
  custom_theme +
  guides(color = guide_legend(title="Dataset"), shape = guide_legend(title="Group")) + # Ensure titles appear
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7)
  )
common_legend <- cowplot::get_legend(legend_plot)

combined_plots <- ggarrange(
  p_tusco,
  p_sequin,
  p_mane,
  ncol = 1,
  nrow = 3,
  align = "v"
)

final_plot <- ggarrange(
  combined_plots,
  common_legend,
  ncol = 1,
  nrow = 2,
  heights = c(15, 1.5)
)

# Do not display the plot; rely on ggsave to write the file
# print(final_plot)

# Save the final figure and TSV under this figure folder
plot_dir <- file.path(fig_dir, "plots")
tsv_dir  <- file.path(fig_dir, "tables")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
if (!dir.exists(tsv_dir))  dir.create(tsv_dir,  recursive = TRUE)

pdf_path <- file.path(plot_dir, "figure3c.pdf")
tsv_path <- file.path(tsv_dir,  "figure3c.tsv")

ggsave(
  filename = pdf_path,
  plot = final_plot,
  width = 1.5,
  height = 5,
  dpi = 600
)

if (exists("all_data") && is.data.frame(all_data)) {
  all_out <- all_data %>% mutate(figure_id = "fig-3", panel_id = "3c")
  readr::write_tsv(all_out, tsv_path)
}

# === EXTRACTION FOR PAPER VALUES ===
# Calculate explicit correlations for TUSCO and SIRVs/Sequins
extraction_file <- file.path(repo_root, "manuscript", "extraction_values.tsv")

# Initialize correlation values (ensure numeric type to avoid formatC errors)
tusco_r <- NA_real_
tusco_p <- NA_real_
sequins_r <- NA_real_
sequins_p <- NA_real_

# Calculate correlations if we have the combined data
if (exists("combined_data") && nrow(combined_data) > 0) {
  # TUSCO correlation
  tusco_subset <- combined_data %>% filter(Type == "TUSCO")
  if (nrow(tusco_subset) > 2) {
    cor_test <- cor.test(tusco_subset$RIN_value, tusco_subset$Percentage_match)
    tusco_r <- as.numeric(cor_test$estimate)
    tusco_p <- as.numeric(cor_test$p.value)
  }
  
  # Sequins correlation  
  sequins_subset <- combined_data %>% filter(Type == "Sequins")
  if (nrow(sequins_subset) > 2) {
    cor_test <- cor.test(sequins_subset$RIN_value, sequins_subset$Percentage_match)
    sequins_r <- as.numeric(cor_test$estimate)
    sequins_p <- as.numeric(cor_test$p.value)
  }
}

# Prepare extraction data
extraction_data <- data.frame(
  figure_id = "fig-3c",
  metric_name = c("TUSCO_R_correlation", "TUSCO_p_value", "Sequins_R_correlation", "Sequins_p_value"),
  value = c(
    formatC(tusco_r, format="f", digits=3),
    formatC(tusco_p, format="e", digits=2),
    formatC(sequins_r, format="f", digits=3),
    formatC(sequins_p, format="f", digits=2)
  ),
  ci_lower = NA,
  ci_upper = NA,
  notes = "Pearson correlation: RIN vs TP/(TP+PTP)",
  stringsAsFactors = FALSE
)

# Ensure output directory exists for extraction file
extraction_dir <- dirname(extraction_file)
if (!dir.exists(extraction_dir)) {
  dir.create(extraction_dir, recursive = TRUE, showWarnings = FALSE)
}

# Write or append to file
if (!file.exists(extraction_file)) {
  write.table(extraction_data, extraction_file, sep="\t", 
              row.names=FALSE, quote=FALSE)
  cat("Created extraction file:", extraction_file, "\n")
} else {
  write.table(extraction_data, extraction_file, sep="\t", 
              row.names=FALSE, quote=FALSE, append=TRUE, col.names=FALSE)
}
cat("Extracted correlations: TUSCO R =", formatC(tusco_r, format="f", digits=3), 
    ", p =", formatC(tusco_p, format="e", digits=2),
    "; Sequins R =", formatC(sequins_r, format="f", digits=3),
    ", p =", formatC(sequins_p, format="f", digits=2), "\n")
# === END EXTRACTION ===

message(paste0("Script completed. Output PDF: ", pdf_path))
message(paste0("Output TSV: ", tsv_path, " (if data available)"))
message(paste0("IMPORTANT: Ensure 'main_data_dir' (currently ", shQuote(main_data_dir) ,") points to your TUSCO classification files for RIN analysis."))
