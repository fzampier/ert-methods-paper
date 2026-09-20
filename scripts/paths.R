# Path helpers sourced by every R script. Run all scripts from the repository
# root, for example:
#   Rscript scripts/make_ertb_fair_comparators.R
#
# The ARMA and FACTT scripts need patient-level BioLINCC data, which are not
# part of this repository (see data/README.md). Put the files under
# data/derived/biolincc_validated/, or set the environment variable
# ERT_BIOLINCC_DIR to the folder that holds them.

paper_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    has_layout <- file.exists(file.path(current, "scripts", "paths.R")) &&
      dir.exists(file.path(current, "R")) &&
      dir.exists(file.path(current, "tables"))

    if (has_layout) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find the repository root from ", start, call. = FALSE)
    }
    current <- parent
  }
}

validated_biolincc_dir <- function(root = paper_root()) {
  env_dir <- Sys.getenv("ERT_BIOLINCC_DIR", "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = TRUE))
  }
  file.path(root, "data", "derived", "biolincc_validated")
}

validated_biolincc_file <- function(..., root = paper_root(), must_work = TRUE) {
  path <- file.path(validated_biolincc_dir(root), ...)
  if (must_work && !file.exists(path)) {
    stop("Missing BioLINCC-derived file: ", path, " (see data/README.md)", call. = FALSE)
  }
  path
}

validated_trial_csv <- function(trial = c("arma", "factt"), root = paper_root()) {
  trial <- match.arg(trial)
  validated_biolincc_file("csv_for_analysis", paste0(trial, ".csv"), root = root)
}

paper_file <- function(..., root = paper_root(), must_work = TRUE) {
  path <- file.path(root, ...)
  if (must_work && !file.exists(path)) {
    stop("Missing repository file: ", path, call. = FALSE)
  }
  path
}

local_ert_file <- function(file, root = paper_root()) {
  paper_file("R", file, root = root)
}
