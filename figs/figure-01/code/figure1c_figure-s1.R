# Function to safely load packages
safe_load <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    # Try installing if not found
    # install.packages(package_name, dependencies = TRUE)
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      warning(paste("Package", package_name, "not available and couldn't be installed. Some functionality may be limited."))
      return(FALSE)
    }
  }
  return(TRUE)
}

# Required packages
suppressPackageStartupMessages({
  # Core packages
  safe_load("rtracklayer")
  safe_load("dplyr")
  safe_load("tidyr")
  safe_load("ggplot2")
  safe_load("stringr")
  safe_load("scales")
  safe_load("Biostrings")  # For GC content analysis
  safe_load("BSgenome")    # For genome sequences
  safe_load("GenomicRanges") # For GRanges used in sequence extraction
  safe_load("IRanges")       # For IRanges used in sequence extraction

  # Optional packages for layout
  has_cowplot <- safe_load("cowplot")
  has_gridExtra <- safe_load("gridExtra")
  has_patchwork <- safe_load("patchwork")
})

###############################################################
# Custom color palette & Define Dataset Groups
###############################################################
palette_colors <- c(
  # TUSCO group - Greens
  "Human (TUSCO)"   = "#a8d5a0",  # dark green
  "Mouse (TUSCO)"   = "#1b9e77",  # light green

  # GENCODE/MANE/RefSeq group - Oranges/Reds
  "Human (GENCODE)" = "#fdbf6f",  # light orange
  "Mouse (GENCODE)" = "#e66101",  # darker orange
  "Human (MANE)"    = "#e41a1c",  # NEW: Red
  "Human (RefSeq)"  = "#a6761d",  # NEW: Golden Brown
  "Mouse (RefSeq)"  = "#663d00",  # NEW: Dark Brown

  # SIRVs/ERCCs/Sequin group - Purples (not plotted in main figures but kept for potential use)
  "SIRVs"           = "#cab2d6",  # light purple
  "ERCCs"           = "#6a3d9a",  # darker purple
  "Sequins"         = "#413aa3"   # medium purple 
)

# Datasets for the legend (excluding RefSeq as per image)
human_datasets_legend <- c("Human (TUSCO)", "Human (GENCODE)", "Human (MANE)")
mouse_datasets_legend <- c("Mouse (TUSCO)", "Mouse (GENCODE)")
spikein_datasets_legend <- c("SIRVs", "ERCCs", "Sequins") # Keep all spike-ins

# Combined groups for plotting (still include RefSeq in plots if data exists)
human_datasets_plot <- c("Human (TUSCO)", "Human (GENCODE)", "Human (MANE)", "Human (RefSeq)")
mouse_datasets_plot <- c("Mouse (TUSCO)", "Mouse (GENCODE)", "Mouse (RefSeq)")
spikein_datasets_plot <- c("SIRVs", "ERCCs", "Sequins")

human_plus_spikeins_plot <- c(human_datasets_plot, spikein_datasets_plot)
mouse_plus_spikeins_plot <- c(mouse_datasets_plot, spikein_datasets_plot)

# Define the specific order for the legend based on the image (filled by row)
# Updated legend order for the 6-panel figure (excluding spike-ins)
# legend_order <- c(
#     "Human (BUGSI)", "Human (GENCODE)", "Human (MANE)", "Human (RefSeq)",
#     "Mouse (BUGSI)", "Mouse (GENCODE)", "Mouse (RefSeq)"
# )
# New order based on target image (2 rows)
legend_order <- c(
    "Human (TUSCO)", "Mouse (TUSCO)", "Human (GENCODE)", "Mouse (GENCODE)",
    "Human (RefSeq)", "Mouse (RefSeq)", "Human (MANE)"
)

# Combine all dataset levels relevant for plotting factor ordering (includes RefSeq)
all_plot_datasets_ordered <- unique(c(human_datasets_plot, mouse_datasets_plot, spikein_datasets_plot))

# Define the order for exon count categories
exon_levels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "10+")

###############################################################
# Define input files (support repo-local and absolute figs/data)
###############################################################
# Prefer new data layout under data/raw, then legacy figs/data
data_base_candidates <- c("data/raw", "figs/data")

project_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(6)) {
    if (file.exists(file.path(cur, ".git")) || file.exists(file.path(cur, "config", "project.yml"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  return(normalizePath(start, winslash = "/", mustWork = FALSE))
}

resolve_data <- function(path_in) {
  # Accept full paths like 'figs/data/...' or a relative fragment like 'reference/human/...'
  rel <- sub("^figs/data/", "", path_in)
  rel <- sub("^data/raw/", "", rel)
  # Try as-is first
  if (file.exists(path_in) || dir.exists(path_in)) return(path_in)
  # Try candidate bases relative to CWD and repo root
  for (base in data_base_candidates) {
    cand <- file.path(base, rel)
    if (file.exists(cand) || dir.exists(cand)) return(cand)
    root_cand <- file.path(project_root(), base, rel)
    if (file.exists(root_cand) || dir.exists(root_cand)) return(root_cand)
  }
  return(path_in) # fallback
}

# TUSCO annotations (GTF)
tusco_human_gtf   <- resolve_data("figs/data/tusco/tusco_human.gtf")
tusco_mouse_gtf   <- resolve_data("figs/data/tusco/tusco_mouse.gtf")

# Reference annotations (GTF/GTF.GZ)
human_gencode_gtf <- resolve_data("figs/data/reference/human/gencode.v49.annotation.gtf.gz")
mouse_gencode_gtf <- resolve_data("figs/data/reference/mouse/gencode.vM38.annotation.gtf.gz")
mane_gtf          <- resolve_data("figs/data/reference/human/mane.v1.4.ensembl_genomic.gtf.gz")  # human only
human_refseq_gtf  <- resolve_data("figs/data/reference/human/refseq/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz")
mouse_refseq_gtf  <- resolve_data("figs/data/reference/mouse/refseq/GCF_000001635.27_GRCm39_genomic.gtf.gz")

# Spike-ins
ercc_sirv_gtf     <- resolve_data("figs/data/spike-ins/lrgasp_sirvs4.gtf")
sequin_gtf        <- resolve_data("figs/data/spike-ins/rnasequin_annotation_2.4.gtf")

# Genome sequences for GC content calculation
human_genome_fasta <- resolve_data("figs/data/reference/human/lrgasp_grch38_sirvs.fasta")
mouse_genome_fasta <- resolve_data("figs/data/reference/mouse/mm39_SIRV.fa")
# Support both local and absolute figs/data
sequin_genome_fasta <- resolve_data("figs/data/reference/rnasequin_decoychr_2.4.fa")

###############################################################
# Functions for data processing (largely unchanged, minor safety checks)
###############################################################
import_and_standardize <- function(gtf_file, dataset_type, expected_species = NULL) {
  message("Importing ", dataset_type, " from ", gtf_file)

  if (!file.exists(gtf_file)) {
    warning("File not found: ", gtf_file)
    return(NULL)
  }

  tryCatch({
    gtf <- rtracklayer::import(gtf_file)
    df <- as.data.frame(gtf)

    # Basic check for empty GTF data
    if(nrow(df) == 0) {
        warning("GTF file is empty or could not be parsed: ", gtf_file)
        return(NULL)
    }

    # Standardize gene_id and transcript_id
    if (dataset_type == "TUSCO") {
      if(!"ensembl" %in% colnames(df)) {
        warning("TUSCO GTF expected to have 'ensembl' column for gene_id: ", gtf_file)
        # Attempt to find alternative if possible, e.g., gene_id or transcript_id if present
        if("gene_id" %in% colnames(df)) df$ensembl <- df$gene_id
        else if("transcript_id" %in% colnames(df)) df$ensembl <- df$transcript_id
        else return(NULL) # Cannot proceed without an identifier
      }
      df <- df %>%
        mutate(gene_id = ensembl,
               transcript_id = ensembl) # TUSCO assumes single isoform per 'ensembl' ID
    } else {
        # Ensure gene_id exists
        if (!"gene_id" %in% colnames(df)) {
            # Try alternatives
            if ("gene" %in% colnames(df)) df$gene_id <- df$gene
            else if ("Name" %in% colnames(df) && any(grepl("gene", df$type, ignore.case = TRUE))) {
                # Heuristic: Try to extract from Name attribute if type is gene-like
                gene_entries <- df %>% filter(grepl("gene", type, ignore.case=TRUE))
                if(nrow(gene_entries) > 0 && "Name" %in% colnames(gene_entries)) {
                   # This requires more complex logic to map back to transcripts/exons, skipping for now
                   warning("Complex gene_id inference needed for ", gtf_file, ". Requires 'gene_id' column.")
                   return(NULL)
                } else {
                   warning("No 'gene_id' column or suitable alternative found in non-TUSCO GTF: ", gtf_file)
                   return(NULL)
                }
            } else {
                warning("No 'gene_id' column or suitable alternative found in non-TUSCO GTF: ", gtf_file)
                return(NULL)
            }
        }
        # Ensure transcript_id exists
        if (!"transcript_id" %in% colnames(df)) {
            if ("transcript" %in% colnames(df)) df$transcript_id <- df$transcript
            else if ("rna_id" %in% colnames(df)) df$transcript_id <- df$rna_id # Handle NCBI GTF
            else if ("transcript_name" %in% colnames(df)) df$transcript_id <- df$transcript_name
            else {
                warning("No 'transcript_id' column or suitable alternative found in non-TUSCO GTF: ", gtf_file)
                return(NULL)
            }
        }
    }

    # Minimal check for necessary columns after standardization
    if (!all(c("gene_id", "transcript_id", "type", "width", "start", "end", "seqnames") %in% colnames(df))) {
      warning("Standardized dataframe missing essential columns for ", dataset_type, " from ", gtf_file)
      return(NULL)
    }

    return(df)
  }, error = function(e) {
    warning("Error importing ", dataset_type, " from ", gtf_file, ": ", e$message)
    return(NULL)
  })
}

extract_features <- function(df, dataset_name) {
  if (is.null(df) || nrow(df) == 0) {
    warning("Empty or NULL dataframe provided for feature extraction: ", dataset_name)
    return(NULL)
  }

  # Ensure essential columns exist
  required_cols <- c("type", "transcript_id", "gene_id", "width", "start", "end")
  if (!all(required_cols %in% colnames(df))) {
      missing_cols <- setdiff(required_cols, colnames(df))
      warning("Dataframe for ", dataset_name, " is missing required columns: ", paste(missing_cols, collapse=", "))
      return(NULL)
  }

  # Filter exons
  exons <- df %>% filter(type == "exon") %>% distinct(transcript_id, start, end, .keep_all = TRUE) # Ensure unique exons per transcript

  if (nrow(exons) == 0) {
    warning("No features of type 'exon' found in dataset: ", dataset_name)
    return(NULL) # Return NULL if no exons, can't calculate features
  }

  # Calculate exon lengths
  exons <- exons %>% mutate(exon_length = width)

  # Transcript length (sum of unique exon lengths per transcript)
  transcript_lengths <- exons %>%
    group_by(transcript_id) %>%
    summarize(
      transcript_length = sum(exon_length, na.rm = TRUE),
      gene_id = dplyr::first(gene_id), # Take the first gene_id associated with the transcript
      .groups = "drop"
    ) %>% filter(transcript_length > 0) # Remove transcripts with zero length

  if (nrow(transcript_lengths) == 0) {
      warning("No valid transcript lengths calculated for: ", dataset_name)
      # Don't return NULL entirely, maybe other features are valid
  }

  # Log median transcript length for this dataset
  if (nrow(transcript_lengths) > 0) {
      med_len <- stats::median(transcript_lengths$transcript_length, na.rm = TRUE)
      n_tx <- nrow(transcript_lengths)
      # Use scales::comma if available for nicer formatting
      pretty_med <- tryCatch(scales::comma(round(med_len)), error = function(...) round(med_len))
      pretty_n <- tryCatch(scales::comma(n_tx), error = function(...) n_tx)
      message("[Median transcript length] ", dataset_name, ": ", pretty_med, " bp (n=", pretty_n, " transcripts)")
  }

  # Exon counts per transcript
  exon_counts <- exons %>%
    group_by(transcript_id) %>%
    summarize(
      exon_count = dplyr::n(),
      gene_id = dplyr::first(gene_id),
      .groups = "drop"
    )

  if (nrow(exon_counts) == 0) {
      warning("No valid exon counts calculated for: ", dataset_name)
  }

  # Introns: gaps between consecutive, ordered exons within the same transcript
  introns <- exons %>%
    filter(!is.na(transcript_id) & !is.na(start) & !is.na(end)) %>% # Ensure valid coordinates
    arrange(transcript_id, start) %>%
    group_by(transcript_id) %>%
    mutate(next_exon_start = dplyr::lead(start)) %>%
    filter(!is.na(next_exon_start)) %>% # Only consider exons that have a following exon
    mutate(intron_length = (next_exon_start - 1) - end) %>% # Calculate gap
    filter(intron_length > 0) %>% # Keep only positive gaps (introns)
    ungroup() %>%
    select(transcript_id, intron_length, gene_id) # Keep relevant columns

  if (nrow(introns) == 0) {
      # This is expected for single-exon transcripts, so maybe just message, not warn
      message("No introns calculated for: ", dataset_name, " (potentially many single-exon transcripts)")
  }

  # Transcripts per gene (based on transcripts with calculated lengths)
  transcripts_per_gene <- transcript_lengths %>%
    filter(!is.na(gene_id)) %>%
    group_by(gene_id) %>%
    summarize(
      transcript_count = dplyr::n(),
      .groups = "drop"
    )

  # Add dataset name to all resulting data frames
  if(nrow(transcript_lengths) > 0) transcript_lengths$dataset <- dataset_name
  if(nrow(introns) > 0) introns$dataset <- dataset_name
  if(nrow(exons) > 0) exons$dataset <- dataset_name # Keep original exons data if needed elsewhere
  if(nrow(exon_counts) > 0) exon_counts$dataset <- dataset_name
  if(nrow(transcripts_per_gene) > 0) transcripts_per_gene$dataset <- dataset_name

  # Return list, handling cases where some features might be empty
  feature_list <- list(
    transcripts = if(nrow(transcript_lengths) > 0) transcript_lengths else NULL,
    introns = if(nrow(introns) > 0) introns else NULL,
    exons = if(nrow(exons) > 0) exons else NULL, # Full exon data
    exon_counts = if(nrow(exon_counts) > 0) exon_counts else NULL,
    tpg = if(nrow(transcripts_per_gene) > 0) transcripts_per_gene else NULL
  )

  # Remove NULL elements from the list before returning
  feature_list <- feature_list[!sapply(feature_list, is.null)]
  if(length(feature_list) == 0) return(NULL) # Return NULL if nothing was extracted

  return(feature_list)
}


extract_ercc_sirv <- function(df) {
  # Simplified: Assumes df is already standardized
  if (is.null(df) || nrow(df) == 0) {
    warning("Empty dataframe provided for ERCC/SIRV extraction")
    return(list(ERCC = NULL, SIRV = NULL))
  }
  ercc_df <- df %>% filter(grepl("^ERCC", seqnames, ignore.case = TRUE))
  sirv_df <- df %>% filter(grepl("^SIRV", seqnames, ignore.case = TRUE))
  return(list(ERCC = ercc_df, SIRV = sirv_df))
}

# Compute per-transcript exonic GC content using provided genome sequences
calculate_gc_content <- function(feature_list, dataset_name, genome_seq) {
  if (is.null(genome_seq) || length(genome_seq) == 0) {
    warning("Genome sequences not available for GC content: ", dataset_name)
    return(NULL)
  }
  if (is.null(feature_list) || is.null(feature_list$exons)) {
    warning("No exons available for GC content: ", dataset_name)
    return(NULL)
  }

  exons_df <- feature_list$exons

  required_cols <- c("seqnames", "start", "end", "transcript_id")
  missing_cols <- setdiff(required_cols, colnames(exons_df))
  if (length(missing_cols) > 0) {
    warning("Exon data missing required columns (", paste(missing_cols, collapse = ", "), ") for ", dataset_name)
    return(NULL)
  }

  exons_df <- exons_df %>%
    mutate(seqnames = as.character(seqnames)) %>%
    filter(seqnames %in% names(genome_seq)) %>%
    filter(!is.na(start) & !is.na(end) & end >= start)

  if (nrow(exons_df) == 0) {
    warning("No exons overlapping available genome contigs for ", dataset_name)
    return(NULL)
  }

  gr <- GenomicRanges::GRanges(
    seqnames = exons_df$seqnames,
    ranges = IRanges::IRanges(start = exons_df$start, end = exons_df$end),
    strand = if ("strand" %in% colnames(exons_df)) exons_df$strand else "*"
  )

  exon_seqs <- tryCatch({
    Biostrings::getSeq(genome_seq, gr)
  }, error = function(e) {
    warning("getSeq failed for ", dataset_name, ": ", e$message)
    return(NULL)
  })

  if (is.null(exon_seqs)) return(NULL)

  gc_mat <- Biostrings::letterFrequency(exon_seqs, letters = c("G", "C"), as.prob = FALSE)
  gc_counts <- rowSums(gc_mat)
  exon_lengths <- Biostrings::width(exon_seqs)

  gc_df <- data.frame(
    transcript_id = exons_df$transcript_id,
    gene_id = if ("gene_id" %in% colnames(exons_df)) exons_df$gene_id else NA_character_,
    gc_num = gc_counts,
    len = exon_lengths,
    stringsAsFactors = FALSE
  )

  gc_by_tx <- gc_df %>%
    group_by(transcript_id) %>%
    summarize(
      gc_content = if (sum(len) > 0) sum(gc_num) / sum(len) else NA_real_,
      gene_id = dplyr::first(gene_id),
      .groups = "drop"
    ) %>%
    filter(!is.na(gc_content))

  if (nrow(gc_by_tx) == 0) return(NULL)
  gc_by_tx$dataset <- dataset_name
  return(gc_by_tx)
}

###############################################################
# Load datasets and Extract Features
###############################################################

# Load genome sequences for GC content calculation
message("Loading genome sequences...")
human_genome <- NULL
mouse_genome <- NULL
sequin_genome <- NULL

if (file.exists(human_genome_fasta)) {
  tryCatch({
    human_genome <- readDNAStringSet(human_genome_fasta)
    # Normalize FASTA names to first token before first whitespace (e.g., '>chr1 AC:...' -> 'chr1')
    try({ names(human_genome) <- sub("\\s.*$", "", names(human_genome)) }, silent = TRUE)
    message("Loaded human genome sequences: ", length(human_genome), " sequences")
  }, error = function(e) {
    warning("Could not load human genome sequences: ", e$message)
  })
} else {
  warning("Human genome file not found: ", human_genome_fasta)
}

if (file.exists(mouse_genome_fasta)) {
  tryCatch({
    mouse_genome <- readDNAStringSet(mouse_genome_fasta)
    try({ names(mouse_genome) <- sub("\\s.*$", "", names(mouse_genome)) }, silent = TRUE)
    message("Loaded mouse genome sequences: ", length(mouse_genome), " sequences")
  }, error = function(e) {
    warning("Could not load mouse genome sequences: ", e$message)
  })
} else {
  warning("Mouse genome file not found: ", mouse_genome_fasta)
}

if (file.exists(sequin_genome_fasta)) {
  tryCatch({
    sequin_genome <- readDNAStringSet(sequin_genome_fasta)
    message("Loaded sequin sequences: ", length(sequin_genome), " sequences")
  }, error = function(e) {
    warning("Could not load sequin sequences: ", e$message)
  })
} else {
  warning("Sequin genome file not found: ", sequin_genome_fasta)
}

# Helpers for prior TSV-driven TUSCO subsets retained for reference (unused now)
# read_tusco_transcripts <- function(tsv_path) { ... }
# subset_to_transcripts <- function(df, transcript_ids) { ... }

datasets_list <- list()

# Import base reference annotations
human_gencode_df <- import_and_standardize(human_gencode_gtf, "GENCODE", "human")
mouse_gencode_df <- import_and_standardize(mouse_gencode_gtf, "GENCODE", "mouse")
human_mane_df    <- import_and_standardize(mane_gtf,          "MANE",    "human")
human_refseq_df  <- import_and_standardize(human_refseq_gtf,  "RefSeq",  "human")
mouse_refseq_df  <- import_and_standardize(mouse_refseq_gtf,  "RefSeq",  "mouse")

# TUSCO from provided GTF files (standard GTFs with gene_id/transcript_id)
human_tusco_df  <- import_and_standardize(tusco_human_gtf, "GENCODE", "human")
mouse_tusco_df  <- import_and_standardize(tusco_mouse_gtf, "GENCODE", "mouse")

# Assemble datasets list
datasets_list[["Human (TUSCO)"]]  <- extract_features(human_tusco_df,   "Human (TUSCO)")
datasets_list[["Mouse (TUSCO)"]]  <- extract_features(mouse_tusco_df,   "Mouse (TUSCO)")
datasets_list[["Human (GENCODE)"]] <- extract_features(human_gencode_df, "Human (GENCODE)")
datasets_list[["Mouse (GENCODE)"]] <- extract_features(mouse_gencode_df, "Mouse (GENCODE)")
datasets_list[["Human (MANE)"]]    <- extract_features(human_mane_df,    "Human (MANE)")
datasets_list[["Human (RefSeq)"]]  <- extract_features(human_refseq_df,  "Human (RefSeq)")
datasets_list[["Mouse (RefSeq)"]]  <- extract_features(mouse_refseq_df,  "Mouse (RefSeq)")

# Add GC content to each dataset
for (dataset_name in names(datasets_list)) {
  if (!is.null(datasets_list[[dataset_name]])) {
    # Choose appropriate genome sequence based on dataset
    genome_seq <- NULL
    if (grepl("Human", dataset_name, ignore.case = TRUE)) {
      genome_seq <- human_genome
    } else if (grepl("Mouse", dataset_name, ignore.case = TRUE)) {
      genome_seq <- mouse_genome
    } else if (grepl("Sequins", dataset_name, ignore.case = TRUE)) {
      genome_seq <- sequin_genome
    }
    
    gc_content <- calculate_gc_content(datasets_list[[dataset_name]], dataset_name, genome_seq)
    if (!is.null(gc_content)) {
      datasets_list[[dataset_name]][["gc_content"]] <- gc_content
    }
  }
}

# Spike-ins (Load but don't assign to human/mouse groups initially)
ercc_sirv_df <- import_and_standardize(ercc_sirv_gtf, "ERCC_SIRV")
if (!is.null(ercc_sirv_df)) {
  es_split <- extract_ercc_sirv(ercc_sirv_df)
  datasets_list[["ERCCs"]] <- extract_features(es_split$ERCC, "ERCCs")
  datasets_list[["SIRVs"]] <- extract_features(es_split$SIRV, "SIRVs")
}
datasets_list[["Sequins"]] <- extract_features(import_and_standardize(sequin_gtf, "Sequins"), "Sequins")

# Remove NULL entries from the list (failed loads/extractions)
datasets_list <- datasets_list[!sapply(datasets_list, is.null)]

###############################################################
# Combine Features Separately for Human and Mouse (Main Figures - NO Spike-ins)
###############################################################
combine_features_subset <- function(datasets_list, feature_name, subset_datasets_group, all_factor_levels) {
  # Filter the main list for the subset datasets
  subset_list <- datasets_list[names(datasets_list) %in% subset_datasets_group]

  # Extract the specified feature from each dataset in the subset
  feature_list <- lapply(subset_list, function(x) {
      if (!is.null(x) && feature_name %in% names(x)) {
          return(x[[feature_name]])
      } else {
          return(NULL) # Return NULL if feature doesn't exist for this dataset
      }
  })

  # Filter out NULL feature dataframes
  valid_features <- feature_list[!sapply(feature_list, is.null)]

  if (length(valid_features) == 0) {
    warning("No valid data for feature '", feature_name, "' in the subset: ", paste(subset_datasets_group, collapse=", "))
    return(NULL)
  }

  # Combine the valid feature dataframes
  combined_data <- bind_rows(valid_features)

  # Ensure the dataset column is a factor with the specified levels
  combined_data <- combined_data %>%
    mutate(dataset = factor(dataset, levels = intersect(all_factor_levels, unique(dataset))))

  return(combined_data)
}

# Combine Features for Human ONLY (Main Plots)
transcript_lengths_human <- combine_features_subset(datasets_list, "transcripts", human_datasets_plot, all_plot_datasets_ordered)
intron_lengths_human <- combine_features_subset(datasets_list, "introns", human_datasets_plot, all_plot_datasets_ordered)
exon_counts_human <- combine_features_subset(datasets_list, "exon_counts", human_datasets_plot, all_plot_datasets_ordered)
gc_content_human <- combine_features_subset(datasets_list, "gc_content", human_datasets_plot, all_plot_datasets_ordered)

# Combine Features for Mouse ONLY (Main Plots)
transcript_lengths_mouse <- combine_features_subset(datasets_list, "transcripts", mouse_datasets_plot, all_plot_datasets_ordered)
intron_lengths_mouse <- combine_features_subset(datasets_list, "introns", mouse_datasets_plot, all_plot_datasets_ordered)
exon_counts_mouse <- combine_features_subset(datasets_list, "exon_counts", mouse_datasets_plot, all_plot_datasets_ordered)
gc_content_mouse <- combine_features_subset(datasets_list, "gc_content", mouse_datasets_plot, all_plot_datasets_ordered)

# Define final legend order and apply it to combined datasets
final_legend_order <- c(
    "Human (TUSCO)", "Mouse (TUSCO)", "Human (GENCODE)", "Mouse (GENCODE)",
    "Human (RefSeq)", "Mouse (RefSeq)", "Human (MANE)"
)

set_final_order <- function(df, order_levels) {
  if (!is.null(df) && nrow(df) > 0 && "dataset" %in% colnames(df)) {
    present_levels <- intersect(order_levels, unique(df$dataset))
    if (length(present_levels) > 0) {
      df$dataset <- factor(df$dataset, levels = present_levels)
    }
  }
  return(df)
}

transcript_lengths_human <- set_final_order(transcript_lengths_human, final_legend_order)
intron_lengths_human     <- set_final_order(intron_lengths_human, final_legend_order)
exon_counts_human        <- set_final_order(exon_counts_human, final_legend_order)
gc_content_human         <- set_final_order(gc_content_human, final_legend_order)
transcript_lengths_mouse <- set_final_order(transcript_lengths_mouse, final_legend_order)
intron_lengths_mouse     <- set_final_order(intron_lengths_mouse, final_legend_order)
exon_counts_mouse        <- set_final_order(exon_counts_mouse, final_legend_order)
gc_content_mouse         <- set_final_order(gc_content_mouse, final_legend_order)

###############################################################
# Helpers: per‑dataset medians for density plots
###############################################################
# Generic median summarizer that avoids tidy‑eval dependencies
summarize_feature_median <- function(df, value_col, panel_id, feature_label,
                                     value_transform = NULL, figure_id = "fig-s1") {
  if (is.null(df) || nrow(df) == 0 || !(value_col %in% names(df))) return(NULL)
  # Keep only rows with a dataset and non‑missing values
  ok <- !is.na(df[[value_col]]) & !is.na(df$dataset)
  if (!any(ok)) return(NULL)
  tmp <- df[ok, c("dataset", value_col)]
  # Compute per‑dataset median and N using base aggregate (no tidy‑eval)
  med <- stats::aggregate(tmp[[value_col]], by = list(dataset = tmp$dataset),
                          FUN = function(x) stats::median(x, na.rm = TRUE))
  ns  <- stats::aggregate(tmp[[value_col]], by = list(dataset = tmp$dataset),
                          FUN = function(x) sum(!is.na(x)))
  out <- merge(med, ns, by = "dataset")
  # Standardize column names
  names(out)[names(out) == "x.x"] <- "median"
  names(out)[names(out) == "x.y"] <- "n"
  if (!"median" %in% names(out)) names(out)[names(out) == "x"] <- "median"
  if (!"n" %in% names(out))      names(out)[names(out) == "x"] <- "n"
  # Optional value transform (e.g., GC proportion -> percent)
  if (!is.null(value_transform)) out$median <- value_transform(out$median)
  # Order columns and tag figure/panel/feature
  out <- out[order(out$dataset), , drop = FALSE]
  out$figure_id <- figure_id
  out$panel_id  <- panel_id
  out$feature   <- feature_label
  out <- out[, c("figure_id", "panel_id", "feature", "dataset", "n", "median")]
  rownames(out) <- NULL
  return(out)
}

###############################################################
# Nature-style Theme for Plots (Unchanged)
###############################################################
nature_theme <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0), # Titles aligned left
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Default angle 0 for x-axis text
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 7, face = "bold"),
    axis.line = element_line(linewidth = 0.25), # Use linewidth instead of size
    axis.ticks = element_line(linewidth = 0.25),
    legend.key.size = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom", # Default legend position
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt"), # Small top margin for legend box
    legend.box.background = element_rect(color = "black", linewidth = 0.5) # Add black border
  )

###############################################################
# Define Plotting Functions
###############################################################

# 1. Transcript Length Density Plot
plot_transcript_density <- function(data, title, species_colors) {
  if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Data)")) + nature_theme)
  
  # Filter out zero-length transcripts just in case
  data <- data %>% filter(transcript_length > 0)
  if (nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Positive Length Data)")) + nature_theme)

  ggplot(data, aes(x = transcript_length, color = dataset)) +
    geom_density(linewidth = 0.5, key_glyph = "path") + # Use linewidth, path glyph better for density lines
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 5), # Auto breaks
      labels = scales::trans_format("log10", scales::math_format(10^.x)), # Nicer labels
      limits = c(NA, 1e5) # Ensure upper limit consistent if needed
    ) +
    annotation_logticks(sides = "b", linewidth = 0.25, outside = TRUE, # Add log ticks below axis
                        short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    coord_cartesian(clip = "off") + # Allow ticks outside panel
    scale_color_manual(values = species_colors, name = "Dataset") +
    labs(
      title = title,
      x = "Transcript Length (nt, log scale)",
      y = "Density"
    ) +
    nature_theme +
    theme(plot.margin = margin(t = 5, r = 5, b = 15, l = 5, unit = "pt")) # Adjust margin for ticks
}

# 2. Exon Counts Bar Plot
plot_exon_counts <- function(data, title, species_colors) {
  if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Data)")) + nature_theme)

  # Prepare data for plotting: group >10, calculate frequencies
  exon_data_for_plot <- data %>%
    mutate(exon_count_group = ifelse(exon_count > 10, "10+", as.character(exon_count))) %>%
    mutate(exon_count_group = factor(exon_count_group, levels = exon_levels)) %>%
    # Calculate frequency PER DATASET
    group_by(dataset) %>%
    mutate(total_transcripts_in_dataset = n()) %>%
    ungroup() %>%
    group_by(dataset, exon_count_group, total_transcripts_in_dataset) %>%
    summarize(count = n(), .groups = "drop") %>%
    mutate(freq = (count / total_transcripts_in_dataset) * 100) %>%
    ungroup()

  if (nrow(exon_data_for_plot) == 0) return(ggplot() + ggtitle(paste(title, "(Processing Error)")) + nature_theme)

  ggplot(exon_data_for_plot,
         aes(x = exon_count_group, y = freq, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.2, color = "black") + # Add outline
    scale_fill_manual(values = species_colors, name = "Dataset") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) + # Y starts at 0
    labs(
      title = title,
      x = "Number of Exons per Transcript",
      y = "Frequency (%)"
    ) +
    nature_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate labels if many groups
}

# 3. Intron Length Density Plot
plot_intron_density <- function(data, title, species_colors) {
  if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Intron Data)")) + nature_theme)
  
  # Ensure positive intron lengths
  data <- data %>% filter(intron_length > 0)
  if (nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Positive Length Introns)")) + nature_theme)

  ggplot(data, aes(x = intron_length, color = dataset)) +
    geom_density(linewidth = 0.5, key_glyph = "path") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6), # More breaks for introns
      labels = scales::trans_format("log10", scales::math_format(10^.x))
      # Consider setting limits if needed, e.g., limits = c(10, 1e6)
    ) +
    annotation_logticks(sides = "b", linewidth = 0.25, outside = TRUE,
                        short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = species_colors, name = "Dataset") +
    labs(
      title = title,
      x = "Intron Length (nt, log scale)",
      y = "Density"
    ) +
    nature_theme +
    theme(plot.margin = margin(t = 5, r = 5, b = 15, l = 5, unit = "pt")) # Adjust margin for ticks
}

# 4. GC Content Density Plot
plot_gc_content <- function(data, title, species_colors) {
  if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No GC Data)")) + nature_theme)
  
  # Ensure valid GC content values (0-1)
  data <- data %>% filter(gc_content >= 0 & gc_content <= 1 & !is.na(gc_content))
  if (nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Valid GC Content Data)")) + nature_theme)

  ggplot(data, aes(x = gc_content * 100, color = dataset)) +  # Convert to percentage
    geom_density(linewidth = 0.5, key_glyph = "path") +
    scale_x_continuous(
      breaks = seq(20, 80, 10),
      limits = c(20, 80),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_color_manual(values = species_colors, name = "Dataset") +
    labs(
      title = title,
      x = "GC Content (%)",
      y = "Density"
    ) +
    nature_theme
}


###############################################################
# Generate Plots for Human (No Spike-ins)
###############################################################
# Filter palette for Human datasets only
human_palette <- palette_colors[names(palette_colors) %in% human_datasets_plot]

plot_h_transcript_density <- plot_transcript_density(transcript_lengths_human, "", human_palette)
plot_h_exon_counts <- plot_exon_counts(exon_counts_human, "", human_palette)
plot_h_intron_density <- plot_intron_density(intron_lengths_human, "", human_palette)
plot_h_gc_content <- plot_gc_content(gc_content_human, "", human_palette)

###############################################################
# Generate Plots for Mouse (No Spike-ins)
###############################################################
# Filter palette for Mouse datasets only
mouse_palette <- palette_colors[names(palette_colors) %in% mouse_datasets_plot]

plot_m_transcript_density <- plot_transcript_density(transcript_lengths_mouse, "", mouse_palette)
plot_m_exon_counts <- plot_exon_counts(exon_counts_mouse, "", mouse_palette)
plot_m_intron_density <- plot_intron_density(intron_lengths_mouse, "", mouse_palette)
plot_m_gc_content <- plot_gc_content(gc_content_mouse, "", mouse_palette)

###############################################################
# Generate Specific Transcript Length Plot (Human/Mouse TUSCO/GENCODE/RefSeq + ERCCs + SIRVs)
###############################################################
message("Generating specific transcript length plot...")

# Define the target datasets for the specific plot
target_datasets_specific <- c(
  "Human (TUSCO)",
  "Mouse (TUSCO)",
  "Human (GENCODE)",
  "Mouse (GENCODE)",
  "SIRVs",
  "ERCCs"
)

# Extract and combine transcript length data for target datasets
transcript_lengths_specific_list <- lapply(target_datasets_specific, function(ds_name) {
  if (!is.null(datasets_list[[ds_name]]) && "transcripts" %in% names(datasets_list[[ds_name]])) {
    return(datasets_list[[ds_name]][["transcripts"]])
  } else {
    warning("Transcript data not found for: ", ds_name)
    return(NULL)
  }
})

# Filter out NULLs and bind rows
valid_transcript_lengths_specific <- transcript_lengths_specific_list[!sapply(transcript_lengths_specific_list, is.null)]
if (length(valid_transcript_lengths_specific) > 0) {
  transcript_lengths_specific <- bind_rows(valid_transcript_lengths_specific)

  # Ensure dataset is a factor with the correct levels (subset of the full list)
  transcript_lengths_specific <- transcript_lengths_specific %>%
      mutate(dataset = factor(dataset, levels = intersect(target_datasets_specific, unique(dataset))))

  # Filter the palette for the target datasets
  specific_palette <- palette_colors[names(palette_colors) %in% target_datasets_specific]

  # Create the plot (legend removed for fig-1c)
  plot_specific_transcript_density <- plot_transcript_density(
    transcript_lengths_specific,
    "", # No title
    specific_palette
  ) + theme(legend.position = "none")

  # Define output path and save (fig-1c)
  # Save only under this figure folder
  specific_plot_filename <- "figs/figure-01/plots/fig-1c.pdf"
  # Ensure output directory exists
  try({ dir.create(dirname(specific_plot_filename), recursive = TRUE, showWarnings = FALSE) }, silent = TRUE)
  # Using a smaller size, e.g., half-width, adjusted height
  ggsave(
      specific_plot_filename,
      plot_specific_transcript_density,
      width = 85 / 25.4, # Approx half Nature column width in inches
      height = 40 / 25.4, # Reduced height in inches
      dpi = 300,
      device = "pdf"
  )
  message("Saved specific transcript length plot: ", specific_plot_filename)

  # Also write TSV with underlying data + metadata for fig-1c
  try({
    tsv_dir <- "figs/figure-01/tables"
    dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
    fig1c_tsv <- transcript_lengths_specific %>%
      mutate(
        figure_id = "fig-1c",
        panel_id = NA_character_
      ) %>%
      select(figure_id, panel_id, dataset, gene_id, transcript_id, transcript_length)
    # Add per-dataset summary rows as a separate block (n, median)
    summaries <- transcript_lengths_specific %>%
      group_by(dataset) %>%
      summarize(
        n_transcripts = dplyr::n(),
        median_transcript_length = stats::median(transcript_length, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(figure_id = "fig-1c", panel_id = NA_character_) %>%
      select(figure_id, panel_id, dataset, n_transcripts, median_transcript_length)
    # Write two sections by appending; consumers can read both
    out_file <- file.path(tsv_dir, "fig-1c.tsv")
    suppressWarnings(write.table(fig1c_tsv, out_file, sep = "\t", quote = FALSE, row.names = FALSE))
    suppressWarnings(write("\n# summaries\n", out_file, append = TRUE))
    suppressWarnings(write.table(summaries, out_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
    message("Saved data: ", out_file)
  }, silent = TRUE)

} else {
  warning("No valid transcript data found for the specific plot datasets. Skipping generation.")
}

###############################################################
# Arrange and Save Figures (Using Patchwork ONLY)
###############################################################

# Helper function using patchwork for arrangement and legend
arrange_and_save_combined <- function(plots_list, filename_base, fig_width_mm, fig_height_mm, ncol = 2, labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')) {

  # Ensure patchwork is available
  if (!requireNamespace("patchwork", quietly = TRUE)) {
     stop("patchwork package is required for arranging plots. Please install it.")
  }
 
  # Remove individual legends
  plots_no_legend <- lapply(plots_list, function(p) p + theme(legend.position = "none"))

  # Arrange using patchwork - support both 6 and 8 plots
  if (length(plots_no_legend) == 6 && ncol == 2) {
      p1 <- plots_no_legend[[1]] # Human Transcripts
      p2 <- plots_no_legend[[2]] # Mouse Transcripts
      p3 <- plots_no_legend[[3]] # Human Exons
      p4 <- plots_no_legend[[4]] # Mouse Exons
      p5 <- plots_no_legend[[5]] # Human Introns
      p6 <- plots_no_legend[[6]] # Mouse Introns

      # Arrange in 2 columns, 3 rows
      arranged_plots <- (p1 | p2) / (p3 | p4) / (p5 | p6)

  } else if (length(plots_no_legend) == 8 && ncol == 2) {
      p1 <- plots_no_legend[[1]] # Human Transcripts
      p2 <- plots_no_legend[[2]] # Mouse Transcripts
      p3 <- plots_no_legend[[3]] # Human Exons
      p4 <- plots_no_legend[[4]] # Mouse Exons
      p5 <- plots_no_legend[[5]] # Human Introns
      p6 <- plots_no_legend[[6]] # Mouse Introns
      p7 <- plots_no_legend[[7]] # Human GC Content
      p8 <- plots_no_legend[[8]] # Mouse GC Content

      # Arrange in 2 columns, 4 rows
arranged_plots <- (p1 | p2) / (p3 | p4) / (p5 | p6) / (p7 | p8)

  } else {
     stop("Plot arrangement logic currently supports 6 or 8 plots in 2 columns. Provided: ", length(plots_no_legend), " plots.")
  }

  # Add common legend and labels
  final_figure <- arranged_plots +
      plot_layout(guides = 'collect') +
      plot_annotation(tag_levels = 'a') & # Add labels
      theme(legend.position = 'bottom',       # Position collected legend at the bottom
            legend.direction = "horizontal",  # Arrange items horizontally
            legend.box = "horizontal",       # Box layout
            legend.title = element_blank(),   # No legend title
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.4, "lines"), # Slightly smaller key size
            legend.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), # Increased margin around legend content
            legend.box.background = element_rect(color = "black", linewidth = 0.5), # Add black border
            plot.tag = element_text(size = 7, face = "bold") # Label styling
           ) &
      guides(color = "none", fill = guide_legend(nrow = 2)) # Disable COLOR guide, collect+format FILL guide

  # Save the final figure (applies to both patchwork and fallback)
  width_in <- fig_width_mm / 25.4
  height_in <- fig_height_mm / 25.4

  pdf_filename <- paste0(filename_base, ".pdf")

  # Ensure output directory exists
  try({ dir.create(dirname(pdf_filename), recursive = TRUE, showWarnings = FALSE) }, silent = TRUE)

  # Use ggsave for patchwork object directly
  ggsave(pdf_filename, final_figure, width = width_in, height = height_in, dpi = 300, device = "pdf")
  message("Saved: ", pdf_filename)
}

# Define figure dimensions (adjusted for 4 rows instead of 3)
fig_width_mm <- 170 # Full Nature width for 2 columns
fig_height_mm <- 190 # Reduced height for 4 rows + legend

# Combine all plots for the final figure (now including GC content plots)
all_plots_main <- list(
  plot_h_transcript_density, plot_m_transcript_density,
  plot_h_exon_counts, plot_m_exon_counts,
  plot_h_intron_density, plot_m_intron_density,
  plot_h_gc_content, plot_m_gc_content
)

# Arrange and Save Combined 8-Panel Figure (No Spike-ins, including GC content)
arrange_and_save_combined(
  plots_list = all_plots_main,
  filename_base = "figs/figure-01/plots/fig-s1",
  fig_width_mm = fig_width_mm,
  fig_height_mm = fig_height_mm,
  ncol = 2 # Arrange in 2 columns
)

# Prepare and write TSV for fig-s1 capturing data per panel
try({
  tsv_dir <- "figs/figure-01/tables"
  dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)

  # Helper to build exon summary used in the bar plots
  exon_plot_df <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      mutate(exon_count_group = ifelse(exon_count > 10, "10+", as.character(exon_count))) %>%
      mutate(exon_count_group = factor(exon_count_group, levels = exon_levels)) %>%
      group_by(dataset) %>%
      mutate(total_transcripts_in_dataset = n()) %>%
      ungroup() %>%
      group_by(dataset, exon_count_group, total_transcripts_in_dataset) %>%
      summarize(count = n(), .groups = "drop") %>%
      mutate(freq = (count / total_transcripts_in_dataset) * 100) %>%
      ungroup()
  }

  # Build per-panel data blocks
  a <- if (!is.null(transcript_lengths_human) && nrow(transcript_lengths_human) > 0) transcript_lengths_human %>% mutate(panel_id = "a") else NULL
  b <- if (!is.null(transcript_lengths_mouse) && nrow(transcript_lengths_mouse) > 0) transcript_lengths_mouse %>% mutate(panel_id = "b") else NULL
  c <- exon_plot_df(exon_counts_human)
  if (!is.null(c)) c$panel_id <- "c"
  d <- exon_plot_df(exon_counts_mouse)
  if (!is.null(d)) d$panel_id <- "d"
  e <- if (!is.null(intron_lengths_human) && nrow(intron_lengths_human) > 0) intron_lengths_human %>% mutate(panel_id = "e") else NULL
  f <- if (!is.null(intron_lengths_mouse) && nrow(intron_lengths_mouse) > 0) intron_lengths_mouse %>% mutate(panel_id = "f") else NULL
  g <- if (!is.null(gc_content_human) && nrow(gc_content_human) > 0) gc_content_human %>% mutate(panel_id = "g") else NULL
  h <- if (!is.null(gc_content_mouse) && nrow(gc_content_mouse) > 0) gc_content_mouse %>% mutate(panel_id = "h") else NULL

  # Tag figure_id and select columns by panel
  sel_a_b <- function(df) df %>% mutate(figure_id = "fig-s1") %>% select(figure_id, panel_id, dataset, gene_id, transcript_id, transcript_length)
  sel_e_f <- function(df) df %>% mutate(figure_id = "fig-s1") %>% select(figure_id, panel_id, dataset, gene_id, transcript_id, intron_length)
  sel_g_h <- function(df) df %>% mutate(figure_id = "fig-s1") %>% select(figure_id, panel_id, dataset, gene_id, transcript_id, gc_content)
  sel_c_d <- function(df) df %>% mutate(figure_id = "fig-s1") %>% select(figure_id, panel_id, dataset, exon_count_group, total_transcripts_in_dataset, count, freq)

  blocks <- list(
    if (!is.null(a)) sel_a_b(a),
    if (!is.null(b)) sel_a_b(b),
    if (!is.null(c)) sel_c_d(c),
    if (!is.null(d)) sel_c_d(d),
    if (!is.null(e)) sel_e_f(e),
    if (!is.null(f)) sel_e_f(f),
    if (!is.null(g)) sel_g_h(g),
    if (!is.null(h)) sel_g_h(h)
  )
  blocks <- blocks[!sapply(blocks, is.null)]

  if (length(blocks) > 0) {
    fig_s1_tsv <- dplyr::bind_rows(blocks)
    out_file <- file.path(tsv_dir, "fig-s1.tsv")
    suppressWarnings(write.table(fig_s1_tsv, out_file, sep = "\t", quote = FALSE, row.names = FALSE))
    message("Saved data: ", out_file)

    # Build and append per‑dataset medians for density plots
    s1_summaries <- list(
      summarize_feature_median(transcript_lengths_human, "transcript_length", "a", "transcript_length_bp"),
      summarize_feature_median(transcript_lengths_mouse, "transcript_length", "b", "transcript_length_bp"),
      summarize_feature_median(intron_lengths_human,     "intron_length",    "e", "intron_length_bp"),
      summarize_feature_median(intron_lengths_mouse,     "intron_length",    "f", "intron_length_bp"),
      summarize_feature_median(gc_content_human,         "gc_content",       "g", "gc_content_percent", function(x) x * 100),
      summarize_feature_median(gc_content_mouse,         "gc_content",       "h", "gc_content_percent", function(x) x * 100)
    )
    s1_summaries <- s1_summaries[!sapply(s1_summaries, is.null)]
    if (length(s1_summaries) > 0) {
      fig_s1_summaries <- dplyr::bind_rows(s1_summaries)
      suppressWarnings(write("\n# summaries\n", out_file, append = TRUE))
      suppressWarnings(write.table(fig_s1_summaries, out_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))

      # Log medians to console for quick reference
      try({
        for (i in seq_len(nrow(fig_s1_summaries))) {
          msg <- sprintf("[Median %s, panel %s] %s: %s (n=%s)",
                         fig_s1_summaries$feature[i],
                         fig_s1_summaries$panel_id[i],
                         as.character(fig_s1_summaries$dataset[i]),
                         formatC(fig_s1_summaries$median[i], digits = 3, format = "fg"),
                         formatC(fig_s1_summaries$n[i], big.mark = ",", format = "d"))
          message(msg)
        }
      }, silent = TRUE)
    }
  } else {
    warning("No panel data available to write TSV for fig-s1")
  }
}, silent = TRUE)
