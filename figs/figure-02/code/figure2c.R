# Function to safely load packages
safe_load <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      warning(paste("Package", package_name, "not available and couldn't be installed. Some functionality may be limited."))
      return(FALSE)
    }
  }
  return(TRUE)
}

# Required packages
suppressPackageStartupMessages({
  safe_load("rtracklayer")
  safe_load("dplyr")
  safe_load("tidyr")
  safe_load("ggplot2")
  safe_load("stringr")
  safe_load("scales")
})

###############################################################
# Custom color palette for Figure 2c datasets
###############################################################
palette_colors <- c(
  "Human (TUSCO gene set)"   = "#a8d5a0",
  "Mouse (TUSCO gene set)"   = "#1b9e77",
  "Human (GENCODE)" = "#fdbf6f",
  "Mouse (GENCODE)" = "#e66101",
  "SIRVs"           = "#cab2d6",
  "ERCCs"           = "#6a3d9a"
)

###############################################################
# Define input files (support repo-local and absolute figs/data)
###############################################################
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
  rel <- sub("^figs/data/", "", path_in)
  rel <- sub("^data/raw/", "", rel)
  if (file.exists(path_in) || dir.exists(path_in)) return(path_in)
  for (base in data_base_candidates) {
    cand <- file.path(base, rel)
    if (file.exists(cand) || dir.exists(cand)) return(cand)
    root_cand <- file.path(project_root(), base, rel)
    if (file.exists(root_cand) || dir.exists(root_cand)) return(root_cand)
  }
  return(path_in)
}

# TUSCO annotations (GTF)
tusco_human_gtf   <- resolve_data("figs/data/tusco/tusco_human.gtf")
tusco_mouse_gtf   <- resolve_data("figs/data/tusco/tusco_mouse.gtf")

# Reference annotations (GTF/GTF.GZ)
human_gencode_gtf <- resolve_data("figs/data/reference/human/gencode.v49.annotation.gtf.gz")
mouse_gencode_gtf <- resolve_data("figs/data/reference/mouse/gencode.vM38.annotation.gtf.gz")

# Spike-ins
ercc_sirv_gtf     <- resolve_data("figs/data/spike-ins/lrgasp_sirvs4.gtf")

###############################################################
# Functions for data processing
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

    if(nrow(df) == 0) {
        warning("GTF file is empty or could not be parsed: ", gtf_file)
        return(NULL)
    }

    # Standardize gene_id and transcript_id
    if (dataset_type == "TUSCO") {
      if(!"ensembl" %in% colnames(df)) {
        warning("TUSCO GTF expected to have 'ensembl' column for gene_id: ", gtf_file)
        if("gene_id" %in% colnames(df)) df$ensembl <- df$gene_id
        else if("transcript_id" %in% colnames(df)) df$ensembl <- df$transcript_id
        else return(NULL)
      }
      df <- df %>%
        mutate(gene_id = ensembl,
               transcript_id = ensembl)
    } else {
        if (!"gene_id" %in% colnames(df)) {
            if ("gene" %in% colnames(df)) df$gene_id <- df$gene
            else {
                warning("No 'gene_id' column found in non-TUSCO GTF: ", gtf_file)
                return(NULL)
            }
        }
        if (!"transcript_id" %in% colnames(df)) {
            if ("transcript" %in% colnames(df)) df$transcript_id <- df$transcript
            else if ("rna_id" %in% colnames(df)) df$transcript_id <- df$rna_id
            else if ("transcript_name" %in% colnames(df)) df$transcript_id <- df$transcript_name
            else {
                warning("No 'transcript_id' column found in non-TUSCO GTF: ", gtf_file)
                return(NULL)
            }
        }
    }

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

  required_cols <- c("type", "transcript_id", "gene_id", "width", "start", "end")
  if (!all(required_cols %in% colnames(df))) {
      missing_cols <- setdiff(required_cols, colnames(df))
      warning("Dataframe for ", dataset_name, " is missing required columns: ", paste(missing_cols, collapse=", "))
      return(NULL)
  }

  exons <- df %>% filter(type == "exon") %>% distinct(transcript_id, start, end, .keep_all = TRUE)

  if (nrow(exons) == 0) {
    warning("No features of type 'exon' found in dataset: ", dataset_name)
    return(NULL)
  }

  exons <- exons %>% mutate(exon_length = width)

  transcript_lengths <- exons %>%
    group_by(transcript_id) %>%
    summarize(
      transcript_length = sum(exon_length, na.rm = TRUE),
      gene_id = dplyr::first(gene_id),
      .groups = "drop"
    ) %>% filter(transcript_length > 0)

  if (nrow(transcript_lengths) == 0) {
      warning("No valid transcript lengths calculated for: ", dataset_name)
      return(NULL)
  }

  transcript_lengths$dataset <- dataset_name

  return(list(transcripts = transcript_lengths))
}

extract_ercc_sirv <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    warning("Empty dataframe provided for ERCC/SIRV extraction")
    return(list(ERCC = NULL, SIRV = NULL))
  }
  ercc_df <- df %>% filter(grepl("^ERCC", seqnames, ignore.case = TRUE))
  sirv_df <- df %>% filter(grepl("^SIRV", seqnames, ignore.case = TRUE))
  return(list(ERCC = ercc_df, SIRV = sirv_df))
}

###############################################################
# Load datasets and Extract Features
###############################################################
datasets_list <- list()

# Import TUSCO
human_tusco_df  <- import_and_standardize(tusco_human_gtf, "TUSCO", "human")
mouse_tusco_df  <- import_and_standardize(tusco_mouse_gtf, "TUSCO", "mouse")

# Import GENCODE
human_gencode_df <- import_and_standardize(human_gencode_gtf, "GENCODE", "human")
mouse_gencode_df <- import_and_standardize(mouse_gencode_gtf, "GENCODE", "mouse")

# Assemble datasets list
datasets_list[["Human (TUSCO gene set)"]]  <- extract_features(human_tusco_df,   "Human (TUSCO gene set)")
datasets_list[["Mouse (TUSCO gene set)"]]  <- extract_features(mouse_tusco_df,   "Mouse (TUSCO gene set)")
datasets_list[["Human (GENCODE)"]] <- extract_features(human_gencode_df, "Human (GENCODE)")
datasets_list[["Mouse (GENCODE)"]] <- extract_features(mouse_gencode_df, "Mouse (GENCODE)")

# Spike-ins
ercc_sirv_df <- import_and_standardize(ercc_sirv_gtf, "ERCC_SIRV")
if (!is.null(ercc_sirv_df)) {
  es_split <- extract_ercc_sirv(ercc_sirv_df)
  datasets_list[["ERCCs"]] <- extract_features(es_split$ERCC, "ERCCs")
  datasets_list[["SIRVs"]] <- extract_features(es_split$SIRV, "SIRVs")
}

# Remove NULL entries
datasets_list <- datasets_list[!sapply(datasets_list, is.null)]

###############################################################
# Nature-style Theme for Plots
###############################################################
nature_theme <- theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    legend.key.size = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

###############################################################
# Generate Figure 2c: Transcript Length Binned Plot
###############################################################
message("Generating figure 2c...")

# Define the target datasets for the specific plot
target_datasets_specific <- c(
  "Human (TUSCO gene set)",
  "Mouse (TUSCO gene set)",
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

  # Ensure dataset is a factor with the correct levels
  transcript_lengths_specific <- transcript_lengths_specific %>%
      mutate(dataset = factor(dataset, levels = intersect(target_datasets_specific, unique(dataset))))

  # Filter the palette for the target datasets
  specific_palette <- palette_colors[names(palette_colors) %in% target_datasets_specific]

  # Transcript Length Binned Plot function
  plot_transcript_binned <- function(data, title, species_colors) {
    if (is.null(data) || nrow(data) == 0) return(ggplot() + ggtitle(paste(title, "(No Data)")) + nature_theme)

    data <- data %>%
      filter(transcript_length > 0) %>%
      mutate(length_bin = cut(transcript_length,
                              breaks = c(0, 600, 1000, 2000, 4000, Inf),
                              labels = c("<600", "600-1kb", "1kb-2kb", "2kb-4kb", ">4kb"),
                              include.lowest = TRUE, right = FALSE))

    plot_data <- data %>%
      group_by(dataset, length_bin) %>%
      summarize(count = n(), .groups = "drop") %>%
      group_by(dataset) %>%
      mutate(freq = count / sum(count) * 100) %>%
      ungroup()

    ggplot(plot_data, aes(x = length_bin, y = freq, fill = dataset)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, linewidth = 0.2, color = "black") +
      scale_fill_manual(values = species_colors, name = "Dataset") +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      labs(
        title = title,
        x = "Transcript Length (bp)",
        y = "Frequency (%)"
      ) +
      nature_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  # Create the plot (legend removed for fig-2c)
  plot_specific_transcript_density <- plot_transcript_binned(
    transcript_lengths_specific,
    "", # No title
    specific_palette
  ) + theme(legend.position = "none")

  # Define output path and save (fig-2c)
  specific_plot_filename <- "figs/figure-02/plots/fig-2c.pdf"
  try({ dir.create(dirname(specific_plot_filename), recursive = TRUE, showWarnings = FALSE) }, silent = TRUE)
  ggsave(
      specific_plot_filename,
      plot_specific_transcript_density,
      width = 85 / 25.4,
      height = 40 / 25.4,
      dpi = 300,
      device = "pdf"
  )
  message("Saved specific transcript length plot: ", specific_plot_filename)

  # Also write TSV with underlying data + metadata for fig-2c
  try({
    tsv_dir <- "figs/figure-02/tables"
    dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
    fig2c_tsv <- transcript_lengths_specific %>%
      mutate(
        figure_id = "fig-2c",
        panel_id = NA_character_
      ) %>%
      select(figure_id, panel_id, dataset, gene_id, transcript_id, transcript_length)
    summaries <- transcript_lengths_specific %>%
      group_by(dataset) %>%
      summarize(
        n_transcripts = dplyr::n(),
        median_transcript_length = stats::median(transcript_length, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(figure_id = "fig-2c", panel_id = NA_character_) %>%
      select(figure_id, panel_id, dataset, n_transcripts, median_transcript_length)
    out_file <- file.path(tsv_dir, "fig-2c.tsv")
    suppressWarnings(write.table(fig2c_tsv, out_file, sep = "\t", quote = FALSE, row.names = FALSE))
    suppressWarnings(write("\n# summaries\n", out_file, append = TRUE))
    suppressWarnings(write.table(summaries, out_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
    message("Saved data: ", out_file)
  }, silent = TRUE)

} else {
  warning("No valid transcript data found for the specific plot datasets. Skipping generation.")
}


