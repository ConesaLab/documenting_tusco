###########################
### 0. Define Variables ###
###########################

## Data path resolver: prefer data/raw, fallback to legacy figs/data
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

# Paths to TUSCO GTF files (support local and absolute figs/data)
tusco_dir <- resolve_data("figs/data/tusco")
human_gtf_file <- resolve_data("figs/data/tusco/tusco_human.gtf")
# Mouse GTF (fix typo; rely on canonical filename)
mouse_gtf_file <- resolve_data("figs/data/tusco/tusco_mouse.gtf")

# Directories containing expression data (support both bases)
human_dir       <- resolve_data("figs/data/expression/gtex")
mouse_dir       <- resolve_data("figs/data/expression/encode")

# Where to save final plot under this figure folder
figure_base_dir <- file.path("..")
output_file     <- file.path(figure_base_dir, "plots", "fig-1d.pdf")

# Libraries
library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)
library(cowplot)
library(readr)

# Define a small constant for log transformation to avoid log10(0)
LOG_OFFSET <- 1e-4

# Limit threads to avoid OpenMP SHM errors in restricted environments
try({
  data.table::setDTthreads(1)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
}, silent = TRUE)

###########################################
### 1. Identify TUSCO Genes from GTF   ####
###########################################

extract_tusco_genes <- function(gtf_file) {
  # NOTE: Function re-purposed for TUSCO GTFs which contain a single isoform
  # per gene. The Ensembl gene identifier is stored in the metadata column
  # "ensembl" after import with rtracklayer.

  # Prefer rtracklayer if available; otherwise, parse attributes column directly
  df <- tryCatch({
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      as.data.frame(rtracklayer::import(gtf_file))
    } else {
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
      df0 <- readr::read_tsv(
        gtf_file, comment = "#",
        col_names = c("seqnames","source","type","start","end","score","strand","frame","attribute"),
        col_types = cols, progress = FALSE
      )
      # Extract gene_id to emulate ensembl column
      m <- regmatches(df0$attribute, regexpr("gene_id \"[^\"]+\"", df0$attribute))
      df0$gene_id <- sub("gene_id \"([^\"]+)\"", "\\1", m)
      df0
    }
  }, error = function(e) {
    stop("Failed to import GTF ", gtf_file, ": ", e$message)
  })

  if ("ensembl" %in% colnames(df)) {
    gene_ids <- df$ensembl
  } else if ("gene_id" %in% colnames(df)) {
    gene_ids <- df$gene_id
  } else {
    stop("Cannot find Ensembl gene identifiers in ", gtf_file)
  }

  gene_ids <- unique(gsub("\\..*$", "", gene_ids))  # drop version suffix
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]

  return(gene_ids)
}

cat("Extracting TUSCO genes...\n")
human_tusco_genes <- extract_tusco_genes(human_gtf_file)
mouse_tusco_genes <- extract_tusco_genes(mouse_gtf_file)

cat("Human TUSCO genes:", length(human_tusco_genes), "\n")
cat("Mouse TUSCO genes:", length(mouse_tusco_genes), "\n\n")


##################################################################
### 2. Read Human Tissue Data and Compute Median TPM per Gene  ###
##################################################################

human_files <- list.files(human_dir, pattern="\\.gct\\.gz$", full.names=TRUE)

clean_human_tissue_name <- function(filename) {
  tissue_name <- basename(filename)
  tissue_name <- str_replace(tissue_name, "gene_tpm_2017-06-05_v8_", "")
  tissue_name <- str_replace(tissue_name, "\\.gct\\.gz$", "") # Also remove extension
  tissue_name <- str_replace_all(tissue_name, "_", " ")
  tissue_name <- str_to_title(tissue_name)
  return(tissue_name)
}

read_human_gct <- function(file) {
  con <- gzfile(file, "rt")
  lines <- readLines(con)
  close(con)

  if (length(lines) < 4) {
    warning(paste("File", file, "does not have enough lines. Skipping."))
    return(NULL)
  }

  # Find the header line for data columns
  header_line_index <- grep("^Name\\s+Description", lines)
  if (length(header_line_index) == 0) {
     header_line_index <- 3 # Default guess if specific header not found
  }

  # Extract header and data lines based on the found index
  header <- strsplit(lines[header_line_index], "\t")[[1]]
  data_lines <- lines[(header_line_index + 1):length(lines)]
  data_text  <- paste(data_lines, collapse = "\n")


  gct <- tryCatch({
    dt <- fread(text = data_text, sep = "\t", header = FALSE, na.strings = c("NA", "NaN", ""))
    # Ensure the number of columns matches the header length
    if (ncol(dt) != length(header)) {
         warning(paste("Column number mismatch in", file, ". Expected", length(header), "got", ncol(dt), ". Skipping."))
         return(NULL)
    }
    setnames(dt, names(dt), header)
    dt

  }, error = function(e) {
    warning(paste("Error reading file", file, ":", e$message))
    return(NULL)
  })

  if (is.null(gct)) return(NULL)

  gene_col <- "Name"
  desc_col <- "Description"

  # Check if Name and Description columns exist
  if (!all(c(gene_col, desc_col) %in% colnames(gct))) {
      warning(paste("Missing 'Name' or 'Description' column in", file, ". Using first two columns. Actual cols:", paste(colnames(gct), collapse=", ")))
      # Assuming first col is Gene ID, second is Description if standard names aren't present
      gene_col <- colnames(gct)[1]
      desc_col <- colnames(gct)[2]
      if (is.null(gene_col) || is.null(desc_col)) {
         warning(paste("Cannot identify Gene and Description columns in", file, ". Skipping."))
         return(NULL)
      }
  }

  # Identify sample columns (all columns except Gene and Description)
  sample_cols <- setdiff(colnames(gct), c(gene_col, desc_col))
  if (length(sample_cols) == 0) {
      warning(paste("No sample columns found in", file, ". Skipping."))
      return(NULL)
  }

  # Convert sample columns to numeric, coercing errors to NA
  # Using a loop with 'set' for potentially better stability with dynamic columns
  for (col in sample_cols) {
    # Ensure the column exists before trying to modify it
    if (col %in% names(gct)) {
        # Use as.numeric directly; set handles coercion warnings/errors if needed
        set(gct, j = col, value = as.numeric(as.character(gct[[col]])))
    } else {
        warning(paste("Sample column", col, "not found during numeric conversion in file:", file))
    }
  }

  # Melt the data table
  gct_long <- melt(
    gct,
    id.vars        = c(gene_col, desc_col),
    measure.vars   = sample_cols,
    variable.name  = "Sample",
    value.name     = "TPM"
  )

  # Calculate median TPM per gene for this tissue
  gene_medians <- gct_long[
    , .(median_tpm = median(TPM, na.rm = TRUE)), by = c(gene_col) # Group by the actual gene column name
  ]

  # Clean Gene IDs (remove version)
  # Using 'set' for potentially better stability
  set(gene_medians, j = gene_col, value = str_replace(gene_medians[[gene_col]], "\\..*$", ""))
  setnames(gene_medians, gene_col, "GeneID") # Standardize column name

  # Remove rows where median_tpm is NA (can happen if all samples were NA for a gene)
  gene_medians <- gene_medians[!is.na(median_tpm)]

  return(gene_medians)
}

cat("Reading Human GCT files...\n")
human_data_list <- lapply(human_files, function(f) {
  tissue_name <- clean_human_tissue_name(f)
  cat("  Tissue:", tissue_name, "\n")

  dt <- read_human_gct(f)
  if (is.null(dt) || nrow(dt) == 0) {
      warning(paste("No valid data returned for tissue:", tissue_name))
      return(NULL)
  }

  dt[, Tissue := tissue_name]
  return(dt)
})

human_data_list <- human_data_list[!sapply(human_data_list, is.null)]
if (length(human_data_list) == 0) {
    stop("No human data could be read successfully.")
}
human_data      <- rbindlist(human_data_list, fill = TRUE)


# remove tRNAscan genes and ensure GeneID column exists
if ("GeneID" %in% names(human_data)) {
    human_data <- human_data[!grepl("^tRNAscan:", GeneID)]
    cat("Unique Human genes after reading:", uniqueN(human_data$GeneID), "\n\n")
} else {
    warning("Column 'GeneID' not found in human_data after processing GCT files.")
    # Handle error or stop execution if GeneID is critical and missing
}


################################################################
### 3. Read Mouse Tissue Data and Compute Median TPM per Gene ###
################################################################

mouse_files <- list.files(mouse_dir, pattern="\\.tsv$", full.names=TRUE)
# Exclude subdirectories if any accidentally match pattern
mouse_files <- mouse_files[!sapply(mouse_files, dir.exists)]
cat("List of Mouse files being processed:\n")
print(mouse_files)
cat("\n")

clean_mouse_tissue_name <- function(filename) {
  tissue_name <- basename(filename)
  tissue_name <- str_replace(tissue_name, "\\.tsv$", "")
  tissue_name <- str_replace_all(tissue_name, "_", " ")
  tissue_name <- str_to_title(tissue_name)
  return(tissue_name)
}

read_mouse_tsv <- function(file) {
  cat(paste("    Attempting to read:", basename(file), "\n")) # DEBUG
  dt <- NULL # Initialize dt
  tryCatch({
    # Skip lines starting with #, find the header row
    lines <- readLines(file)
    header_line_index <- grep("^Feature ID\\s+Gene Name", lines) # Adjust pattern if header differs
    if(length(header_line_index) == 0) {
        header_line_index <- grep("^Experiment", lines) # Alternative header pattern
         if(length(header_line_index) == 0) {
             # Try finding the first non-comment line as header
             non_comment_lines <- which(!startsWith(lines, "#"))
             if (length(non_comment_lines) > 0) {
                header_line_index <- non_comment_lines[1] -1 # fread uses skip, so index-1
             } else {
                 warning(paste("Cannot find header line in", file))
                 return(NULL)
             }
         } else {
            header_line_index = header_line_index[1] -1 # Use first match
         }
    } else {
       header_line_index = header_line_index[1] -1 # Use first match
    }

    # === TEMPORARILY SIMPLIFIED (NOW REVERTING) ===
    cat("      Attempting fread...\n")
    dt <- fread(file, sep = "\t", header = TRUE, skip = header_line_index, na.strings = c("NA", "NaN", ""))
    # cat("      fread successful for:", basename(file), " Dimensions:", paste(dim(dt), collapse="x"), "\n") # DEBUG

    # === ORIGINAL PROCESSING (RESTORED) ===
    if (is.null(dt) || nrow(dt) == 0) return(NULL)

    # Force a copy to potentially avoid reference-related internal errors
    dt <- copy(dt)
    cat("      Forced copy of data.table for:", basename(file), "\n") # DEBUG

    # Identify potential TPM and Gene ID columns
    tpm_col <- names(dt)[grepl("TPM", names(dt), ignore.case = TRUE)][1]
    gene_id_col <- names(dt)[grepl("Feature ID|Geneid|Gene ID", names(dt), ignore.case = TRUE)][1] # Adjusted pattern

    if (is.na(tpm_col) || is.na(gene_id_col)) {
      warning(paste("Missing required columns (TPM or Gene ID) in", file, ". Found:", paste(names(dt), collapse=", ")))
      return(NULL)
    }

    # Select and rename relevant columns
    cols_to_select <- c(gene_id_col, tpm_col)
    dt_subset <- dt[, .SD, .SDcols = cols_to_select] # Use .SDcols
    setnames(dt_subset, c(gene_id_col, tpm_col), c("FeatureID", "TPM")) # Rename to standard names

    # Data cleaning
    dt_subset <- dt_subset[!grepl("^tRNAscan:", FeatureID)]
    dt_subset[, TPM := as.numeric(TPM)]
    dt_subset <- dt_subset[!is.na(TPM)] # Remove rows with NA TPM

    if (nrow(dt_subset) == 0) {
        warning(paste("No valid data after cleaning for file:", file))
        return(NULL)
    }

    # Aggregate: Calculate median TPM per gene for this tissue
    dt_agg <- dt_subset[, .(median_tpm = median(TPM, na.rm = TRUE)), by = FeatureID]

    # Clean Gene IDs (remove version)
    dt_agg[, FeatureID := str_replace(FeatureID, "\\..*$", "")]

    setnames(dt_agg, "FeatureID", "GeneID") # Standardize column name

    # Remove rows where median_tpm is NA
    dt_agg <- dt_agg[!is.na(median_tpm)]

    return(dt_agg) # Return the aggregated data

    # === Restoring original return value ===
    # return(dt)
    # === END OF RESTORED SECTION ===

  }, error = function(e) {
    warning(paste("Error during fread/initial processing of file", basename(file), ":", e$message))
    return(NULL)
  })

  # If dt is NULL here, it means fread or the tryCatch block failed
  if (is.null(dt)) {
      return(NULL)
  }

  # If we got here in the simplified version, fread worked.
  # We will return the raw data table read by fread just for this test.
  # In the real run, we'd return dt_agg from the commented out section.
  return(dt)
}


cat("Reading Mouse TSV files...\n")
# Initialize an empty list to store results
mouse_data_list <- vector("list", length(mouse_files))
names(mouse_data_list) <- mouse_files # Optional: keep track by filename

# Use a for loop instead of lapply
for (i in seq_along(mouse_files)) {
  f <- mouse_files[i]
  tissue_name <- clean_mouse_tissue_name(f)
  cat("Processing file:", basename(f), "for Tissue:", tissue_name, "\n") # DEBUG

  # Call the reading function
  dt <- read_mouse_tsv(f)

  if (is.null(dt) || nrow(dt) == 0) {
       warning(paste("No valid data returned for tissue:", tissue_name))
       mouse_data_list[[i]] <- NULL # Assign NULL to the list element
  } else {
      dt[, Tissue := tissue_name]
      mouse_data_list[[i]] <- dt # Assign the result to the list element
  }
  cat("Finished processing:", basename(f), "\n\n") # DEBUG
}

mouse_data_list <- mouse_data_list[!sapply(mouse_data_list, is.null)]
if (length(mouse_data_list) == 0) {
    stop("No mouse data could be read successfully.")
}
mouse_data      <- rbindlist(mouse_data_list, fill = TRUE)

if ("GeneID" %in% names(mouse_data)) {
    cat("Unique Mouse genes after reading:", uniqueN(mouse_data$GeneID), "\n\n")
} else {
     warning("Column 'GeneID' not found in mouse_data after processing TSV files.")
     # Handle error or stop execution
}

#####################################################
### 4. Combine Human & Mouse Median Tissue Data   ###
#####################################################
# Add species information BEFORE combining
human_data[, Species := "Human"]
mouse_data[, Species := "Mouse"]

# Combine the data tables
combined_tissue_data <- rbind(human_data, mouse_data, fill = TRUE)

# Ensure no NA GeneIDs or Tissues crept in
combined_tissue_data <- combined_tissue_data[!is.na(GeneID) & !is.na(Tissue)]

cat("Total records in combined_tissue_data:", nrow(combined_tissue_data), "\n")
cat("Unique genes in combined_tissue_data:", uniqueN(combined_tissue_data$GeneID), "\n")
cat("Unique tissues in combined_tissue_data:", uniqueN(combined_tissue_data$Tissue), "\n\n")

#######################################################################
### 5. Calculate Median Across Tissues for TUSCO Genes             ####
#######################################################################

# Identify TUSCO genes within the combined data
combined_tissue_data[, is_tusco := FALSE]
combined_tissue_data[Species == "Human" & GeneID %in% human_tusco_genes, is_tusco := TRUE]
combined_tissue_data[Species == "Mouse" & GeneID %in% mouse_tusco_genes, is_tusco := TRUE]

# Filter for TUSCO genes
tusco_data <- combined_tissue_data[is_tusco == TRUE]

cat("Number of TUSCO gene records (across tissues):", nrow(tusco_data), "\n")

# Calculate the median TPM across all tissues for each TUSCO gene
tusco_median_across_tissues <- tusco_data[
  , .(median_value = median(median_tpm, na.rm = TRUE)), # This is median of (median per tissue)
  by = .(Species, GeneID)
]

# Add a group label for plotting
tusco_median_across_tissues[, Group := "TUSCO"]

cat("Number of unique TUSCO genes with cross-tissue median:", nrow(tusco_median_across_tissues), "\n")
print(summary(tusco_median_across_tissues$median_value))
cat("\n")


#################################################################################################
### 6. Rank Genes per Tissue & Calculate Median Across Tissues for Top 10k Ranks per Tissue ####
#################################################################################################

# Rank genes within each tissue based on median_tpm (descending)
# Handle ties by assigning the average rank ('average' method)
combined_tissue_data[, rank := frank(-median_tpm, ties.method = "average"), by = .(Species, Tissue)]

cat("Ranking complete.\n")

# --- MODIFIED: Filter for top 10,000 ranks per tissue ---
top_10k_ranks_data <- combined_tissue_data[rank <= 10000]
cat("Number of records after filtering for top 10k ranks per tissue:", nrow(top_10k_ranks_data), "\n")

# Calculate the median TPM across tissues for each of the top 10k ranks
rank_median_across_tissues <- top_10k_ranks_data[
  , .(median_value = median(median_tpm, na.rm = TRUE)),
  by = .(Species, rank) # Group by species and the specific rank number
]

# Add a group label for plotting
# --- MODIFIED: Updated group label ---
rank_median_across_tissues[, Group := "Top 10k Genes"] # Changed label

# Rename 'rank' to 'GeneID_or_Rank' conceptually for merging/plotting structure
setnames(rank_median_across_tissues, "rank", "GeneID_or_Rank")

# Rename GeneID in tusco data for consistency before merging
setnames(tusco_median_across_tissues, "GeneID", "GeneID_or_Rank")

cat("Number of top 10k ranks with cross-tissue median:", nrow(rank_median_across_tissues), "\n")
print(summary(rank_median_across_tissues$median_value))
cat("\n")


##########################################
### 7. Prepare Data for Plotting       ###
##########################################

# Combine the TUSCO results and Top 10k Rank-based results
plot_data <- rbind(
  tusco_median_across_tissues[, .(Species, Group, median_value, GeneID_or_Rank)], # Keep GeneID_or_Rank temporarily if needed for debugging, can remove later
  rank_median_across_tissues[, .(Species, Group, median_value, GeneID_or_Rank)],
  fill = TRUE
)

# Remove any rows where median_value might be NA
plot_data <- plot_data[!is.na(median_value)]

# Log-transform the median values for plotting
# Add a small offset to handle zeros before log transformation
plot_data[, log_median_value := log10(median_value + LOG_OFFSET)]

# Create a combined Group for coloring (Species + Original Group)
# --- MODIFIED: Ensure consistent group naming for plotting ---
# Replace space in "Top 10k Genes" for cleaner factor levels/names if needed
plot_data[Group == "Top 10k Genes", Group := "Top_10k_Genes"]
plot_data[, PlotGroup := paste0(Species, "_", Group)]


cat("Plotting data prepared. Records:", nrow(plot_data), "\n")
print(table(plot_data$PlotGroup))


############################################################
### 8. Make Boxplots (TUSCO vs. Top 10k Rank-Based)       ###
############################################################

# --- MODIFIED: Updated group names and colors ---
custom_colors <- c(
  "Human_TUSCO"         = "#A8D5A0",  # Light Green
  "Mouse_TUSCO"         = "#1b9e77",  # Dark Green
  "Human_Top_10k_Genes" = "#D1D3D4",  # Light Grey (like Human Other)
  "Mouse_Top_10k_Genes" = "#808080"   # Darker Grey (like Mouse Other)
)

# --- MODIFIED: Ensure PlotGroup is a factor with updated levels ---
plot_data[, PlotGroup := factor(PlotGroup, levels = names(custom_colors))]

# --- MODIFIED: Updated legend labels and title, changed plot type ---
p <- ggplot(plot_data, aes(x = Group, y = log_median_value, fill = PlotGroup)) +
  # Violin plot layer - Added linewidth
  geom_violin(trim=TRUE, position = position_dodge(width = 0.9), alpha=0.7, linewidth = 0.3) +
  # Boxplot layer inside the violin - Added linewidth
  geom_boxplot(width=0.15, outlier.shape = NA, position = position_dodge(width = 0.9), fill="white", alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~Species, scales = "fixed") + # Keep y-axis scale the same for comparison
  scale_fill_manual(
      values = custom_colors,
      name = "Group", # Legend title
      # --- MODIFIED: Updated legend labels ---
      labels = c("Human TUSCO", "Mouse TUSCO", "Human Top 10k Genes", "Mouse Top 10k Genes")
      ) +
  guides(fill = "none") +
  coord_cartesian(ylim = quantile(plot_data$log_median_value, c(0.01, 0.99), na.rm=TRUE)) + # Zoom in, excluding extreme outliers
  labs(
    x = "",
    y = expression(log[10] ("Median TPM")),
    title = NULL,
    tag = "d"
  ) +
  # --- MODIFIED: Updated x-axis labels ---
  scale_x_discrete(labels = c("TUSCO" = "TUSCO", "Top_10k_Genes" = "Top 10k Genes")) + # Clean up x-axis labels
  # --- MODIFIED: Set base font size to 7pt ---
  theme_classic(base_family = "Helvetica", base_size = 7) +
  # --- MODIFIED: Adjusted line widths ---
  theme(
    axis.line = element_line(color = "black", linewidth = 0.35), # Set axis line width
    axis.ticks = element_line(color = "black", linewidth = 0.35), # Set axis tick width
    # --- MODIFIED: Rely on base_size for these elements ---
    axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5), # Horizontal labels
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5), # Center title
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.35),
    strip.text = element_text(face = "plain"),
    plot.tag = element_text(face = "bold"),
    plot.tag.position = c(0, 1)
  )


################################
### 9. Save the Figure to PDF ###
################################

# --- MODIFIED: Adjust width and height for smaller font size ---
# Ensure output and TSV directories exist
try({ dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE) }, silent = TRUE)
try({ dir.create(file.path(figure_base_dir, "tsv"), recursive = TRUE, showWarnings = FALSE) }, silent = TRUE)

ggsave(
  filename = output_file,
  plot = p,
  device = "pdf",
  width = 85 / 25.4,  # ~85 mm
  height = 40 / 25.4, # ~40 mm
  units = "in",
  dpi = 300
)

cat("Boxplot saved to:", output_file, "\n")

# Also export underlying data to TSV alongside minimal metadata
try({
  tsv_path <- file.path(figure_base_dir, "tsv", "fig-1d.tsv")
  # Keep core columns for plotting replication and metadata
  out_dt <- data.table::copy(plot_data)
  if (!"PlotGroup" %in% names(out_dt)) out_dt[, PlotGroup := paste0(Species, "_", Group)]
  out_dt[, figure_id := "fig-1d"]
  # panel view is a facet per Species; no panel_id needed
  data.table::fwrite(out_dt[, .(figure_id, Species, Group, PlotGroup, median_value, log_median_value)],
                     file = tsv_path, sep = "\t")
  cat("Wrote TSV:", tsv_path, "\n")
}, silent = TRUE)
