## Simple project path helpers (no external deps)
## Usage: source(file.path(project_root(), "R", "paths.R"))

project_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  max_up <- 6
  for (i in seq_len(max_up)) {
    if (file.exists(file.path(cur, ".git")) || file.exists(file.path(cur, "config", "project.yml"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  return(normalizePath(start, winslash = "/", mustWork = TRUE))
}

data_raw <- function(...) file.path(project_root(), "data", "raw", ...)
data_processed <- function(...) file.path(project_root(), "data", "processed", ...)
figs_dir <- function(...) file.path(project_root(), "figs", ...)

ensure_dirs <- function(...) {
  dirs <- list(...)
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  invisible(dirs)
}
