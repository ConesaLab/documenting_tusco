# TUSCO Figure Utilities Library
# Common functions and patterns for figure generation scripts

#' Safe library loading with graceful error handling
#' @param packages Character vector of package names to load
load_required_packages <- function(packages) {
  suppressPackageStartupMessages({
    for (pkg in packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop("Required package '", pkg, "' is not installed. Please install it with: install.packages('", pkg, "')")
      }
    }
  })
}

#' Resolve file paths from multiple candidate locations
#' @param candidates Character vector of candidate file paths
#' @param is_dir Logical, whether to check for directories instead of files
#' @return First existing path, or first candidate if none exist
resolve_path <- function(candidates, is_dir = FALSE) {
  for (p in candidates) {
    if (is.na(p) || is.null(p)) next
    if (!is_dir && file.exists(p)) return(p)
    if (is_dir && dir.exists(p)) return(p)
  }
  return(candidates[[1]])
}

#' Resolve data file paths from standard data locations
#' @param ... Path components to combine
#' @param bases Character vector of base directories to search
#' @return Resolved file path
resolve_data_path <- function(..., bases = c("data/raw", "data/processed", "../data", "../../data", "../../../data")) {
  sub_path <- file.path(...)
  for (base in bases) {
    candidate <- file.path(base, sub_path)
    if (file.exists(candidate)) return(normalizePath(candidate))
  }
  stop("Data file not found: ", sub_path, "\nSearched in: ", paste(bases, collapse = ", "))
}

#' Create standard output directories
#' @param base_dir Base directory for outputs (default: current figure directory)
#' @return List with plot_dir and table_dir paths
create_output_dirs <- function(base_dir = ".") {
  plot_dir <- file.path(base_dir, "plots")
  table_dir <- file.path(base_dir, "tables")
  
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  list(plot_dir = plot_dir, table_dir = table_dir)
}

#' Save figure with standardized naming and optional data export
#' @param plot_obj ggplot object or similar plot object
#' @param fig_id Figure identifier (e.g., "fig-3a", "fig-s2")
#' @param plot_dir Directory for plot outputs
#' @param table_dir Directory for table outputs
#' @param width Plot width in inches (default: Nature single column)
#' @param height Plot height in inches
#' @param data Optional data frame to save as TSV
#' @param data_suffix Optional suffix for data file name
save_figure <- function(plot_obj, fig_id, plot_dir = "plots", table_dir = "tables", 
                       width = 3.35, height = 4.0, data = NULL, data_suffix = NULL) {
  
  # Ensure output directories exist
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save plot
  plot_file <- file.path(plot_dir, paste0(fig_id, ".pdf"))
  ggsave(plot_file, plot_obj, width = width, height = height, 
         units = "in", device = "pdf", dpi = 300)
  message("Saved plot: ", plot_file)
  
  # Save data if provided
  if (!is.null(data)) {
    data_name <- if (is.null(data_suffix)) fig_id else paste0(fig_id, "-", data_suffix)
    data_file <- file.path(table_dir, paste0(data_name, ".tsv"))
    readr::write_tsv(data, data_file)
    message("Saved data: ", data_file)
  }
  
  invisible(plot_file)
}

#' Parse command line arguments with defaults
#' @param defaults Named list of default values
#' @return List of parsed arguments
parse_figure_args <- function(defaults = list(
  out_dir = "..",
  width = 3.35,  # Nature single column width in inches
  height = 4.0
)) {
  args <- commandArgs(trailingOnly = TRUE)
  
  result <- defaults
  if (length(args) > 0) result$out_dir <- args[1]
  if (length(args) > 1) result$width <- as.numeric(args[2])
  if (length(args) > 2) result$height <- as.numeric(args[3])
  
  return(result)
}

#' Determine whether a path is absolute (works cross-platform)
is_absolute_path <- function(path) {
  if (!is.character(path) || length(path) == 0) return(FALSE)
  path <- path[1]
  if (!nzchar(path)) return(FALSE)
  startsWith(path, "/") || grepl("^[A-Za-z]:[\\\\/]", path)
}

#' Normalize a path without requiring it to exist
safe_normalize <- function(path) {
  tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE), error = function(e) path)
}

#' Discover the repository root by searching for a Git directory or config file
find_repo_root <- function(start = getwd(), limit = 8) {
  cur <- safe_normalize(start)
  for (i in seq_len(limit)) {
    if (file.exists(file.path(cur, ".git")) || file.exists(file.path(cur, "config", "project.yml"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  safe_normalize(start)
}

#' Resolve the script directory using commandArgs or current working directory
detect_script_dir <- function() {
  argv <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", argv[grep("^--file=", argv)])
  if (length(script_path) == 1 && nzchar(script_path)) {
    return(dirname(safe_normalize(script_path)))
  }
  # Fallback to calling frame for sourced scripts
  script_frame <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.character(script_frame) && length(script_frame) == 1 && nzchar(script_frame)) {
    return(dirname(safe_normalize(script_frame)))
  }
  safe_normalize(getwd())
}

#' Build a unified context describing figure locations, data roots, and outputs
#'
#' @param defaults Default argument list forwarded to parse_figure_args()
#' @param extra_data_roots Optional additional directories to search for inputs
#' @param ensure_dirs Logical, whether to create plot/table directories eagerly
#' @return List with metadata and helper functions for figure scripts
figure_context <- function(
  defaults = list(out_dir = "..", width = 3.35, height = 4.0),
  extra_data_roots = NULL,
  ensure_dirs = TRUE
) {
  params <- parse_figure_args(defaults = defaults)
  script_dir <- detect_script_dir()
  # Resolve figure directory (parent of code/ when applicable)
  maybe_parent <- safe_normalize(file.path(script_dir, ".."))
  figure_dir <- if (basename(script_dir) == "code" && dir.exists(maybe_parent)) maybe_parent else script_dir
  repo_root <- find_repo_root(figure_dir)
  figs_root <- safe_normalize(file.path(figure_dir, ".."))

  data_roots <- unique(c(
    file.path(figure_dir, "data"),
    file.path(figs_root, "data"),
    repo_root,
    file.path(repo_root, "data"),
    file.path(repo_root, "data", "raw"),
    file.path(repo_root, "data", "processed"),
    extra_data_roots
  ))
  data_roots <- data_roots[!is.na(data_roots) & nzchar(data_roots)]

  output_base <- if (!is.null(params$out_dir) && nzchar(params$out_dir)) {
    if (is_absolute_path(params$out_dir)) params$out_dir else file.path(script_dir, params$out_dir)
  } else figure_dir
  output_base <- safe_normalize(output_base)

  plot_dir <- file.path(output_base, "plots")
  table_dir <- file.path(output_base, "tables")
  if (ensure_dirs) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ctx <- list(
    script_dir = script_dir,
    figure_dir = figure_dir,
    output_base = output_base,
    plot_dir = plot_dir,
    table_dir = table_dir,
    log_dir = output_base,
    params = params,
    repo_root = repo_root,
    data_roots = data_roots
  )

  # Resolve candidate paths across known data roots
  ctx$resolve_input <- function(...,
                                candidates = NULL,
                                type = c("auto", "file", "directory"),
                                search_roots = NULL,
                                required = TRUE) {
    type <- match.arg(type)
    parts <- list(...)
    path_candidates <- character()
    if (length(parts) > 0) {
      path_candidates <- c(path_candidates, do.call(file.path, parts))
    }
    if (!is.null(candidates)) {
      path_candidates <- c(path_candidates, candidates)
    }
    if (length(path_candidates) == 0) {
      stop("resolve_input() requires at least one path candidate")
    }
    search_roots_full <- unique(c(search_roots, ctx$figure_dir, ctx$script_dir, ctx$data_roots))
    search_roots_full <- search_roots_full[!is.na(search_roots_full) & nzchar(search_roots_full)]

    for (cand in path_candidates) {
      if (is_absolute_path(cand)) {
        resolved <- safe_normalize(cand)
        exists_ok <- switch(type,
          file = file.exists(resolved),
          directory = dir.exists(resolved),
          auto = file.exists(resolved) || dir.exists(resolved)
        )
        if (exists_ok) return(resolved)
      }

      # Try relative to each search root
      for (root in search_roots_full) {
        candidate_path <- safe_normalize(file.path(root, cand))
        exists_ok <- switch(type,
          file = file.exists(candidate_path),
          directory = dir.exists(candidate_path),
          auto = file.exists(candidate_path) || dir.exists(candidate_path)
        )
        if (exists_ok) return(candidate_path)
      }
    }

    if (required) {
      stop("Unable to resolve input path from candidates: ", paste(path_candidates, collapse = ", "))
    }
    safe_normalize(path_candidates[[1]])
  }

  # Save plot into plots/ with optional extension inference
  ctx$save_plot <- function(plot_obj, filename, width = ctx$params$width, height = ctx$params$height, units = "in", dpi = 300, ...) {
    if (missing(filename) || !nzchar(filename)) stop("Filename must be supplied to save_plot()")
    ext <- tools::file_ext(filename)
    fname <- if (nzchar(ext)) filename else paste0(filename, ".pdf")
    out_path <- file.path(ctx$plot_dir, fname)
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_path, plot_obj, width = width, height = height, units = units, dpi = dpi, ...)
    message("Saved plot: ", out_path)
    invisible(out_path)
  }

  # Write a tibble/data.frame into tables/
  ctx$write_table <- function(data, filename, ...) {
    if (missing(filename) || !nzchar(filename)) stop("Filename must be supplied to write_table()")
    ext <- tools::file_ext(filename)
    fname <- if (nzchar(ext)) filename else paste0(filename, ".tsv")
    out_path <- file.path(ctx$table_dir, fname)
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(data, out_path, ...)
    message("Saved table: ", out_path)
    invisible(out_path)
  }

  # Open a sink that mirrors stdout/messages into a log file
  ctx$start_log <- function(name = "run.log", split = TRUE, append = FALSE) {
    log_path <- file.path(ctx$log_dir, name)
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    con <- file(log_path, open = if (append) "at" else "wt")
    sink(con, split = split)
    sink(con, type = "message")
    message("Logging to ", log_path)
    invisible(function() {
      sink(type = "message")
      sink()
      try(close(con), silent = TRUE)
    })
  }

  ctx
}

#' Standard TUSCO theme for consistent plot styling
#' @param base_size Base font size
#' @param base_family Font family
theme_tusco <- function(base_size = 7, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      legend.key.size = unit(0.5, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

#' Read TSV files safely with error handling
#' @param file_path Path to TSV file
#' @param ... Additional arguments passed to readr::read_tsv
read_tsv_safe <- function(file_path, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  tryCatch(
    readr::read_tsv(file_path, show_col_types = FALSE, ...),
    error = function(e) {
      stop("Error reading file ", file_path, ": ", e$message)
    }
  )
}

#' Get figure directory from script path
#' Determines the figure directory based on script location
get_figure_dir <- function() {
  # Try to get script path from command line args
  argv <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", argv[grep("^--file=", argv)])
  
  if (length(script_path) == 1 && nzchar(script_path)) {
    # Script path available - figure dir is parent of code dir
    fig_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
  } else {
    # Fallback - assume we're in a code directory
    if (basename(getwd()) == "code") {
      fig_dir <- normalizePath("..", mustWork = FALSE)
    } else {
      fig_dir <- getwd()
    }
  }
  
  return(fig_dir)
}

# Export commonly used color palettes
TUSCO_COLORS <- list(
  human_tusco = "#a8d5a0",
  mouse_tusco = "#1b9e77", 
  human_gencode = "#fdbf6f",
  mouse_gencode = "#e66101",
  human_mane = "#e41a1c",
  sirvs = "#cab2d6",
  erccs = "#6a3d9a"
)

# Export standard dimensions
NATURE_DIMS <- list(
  single_column = 3.35,    # inches
  double_column = 7.09,    # inches
  full_page_width = 7.87,  # inches
  default_height = 4.0     # inches
)
