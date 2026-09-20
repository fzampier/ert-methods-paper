# Capture sessionInfo() for the long-term reproducibility record.
#
# Loads the packages used by the manuscript and trial-example scripts, then
# writes sessionInfo() to sessionInfo.txt at the repository root. Run once on a stable build with
#   Rscript R/capture_session_info.R
# or via `make session-info`.

source("scripts/paths.R")

required_pkgs <- c(
  "tidyverse",  # used by R/ertb.R, R/erte.R
  "ggplot2",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "scales",
  "Rcpp",       # used by scripts/make_factt_permutation_figure.R
  "rpact",      # calibrated Lan-DeMets boundaries (make_ertb_fair_comparators.R)
  "mvtnorm",    # independent MVN boundary recursion + power-vs-K integration
  "patchwork",  # figure composition (make_worked_alert_replay.R)
  "ggh4x"       # per-panel axis scales in Figure 5 (guarded requireNamespace,
                # but present in the production environment)
)

for (pkg in required_pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

root <- paper_root()
out_dir <- paper_file(".", root = root, must_work = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "sessionInfo.txt")

con <- file(out_path, open = "wt")
on.exit(close(con), add = TRUE)

writeLines(c(
  paste("# sessionInfo captured", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  ""
), con)

cap <- capture.output(print(sessionInfo()))
cap <- sub("[[:space:]]+$", "", cap)  # R pads package-list lines; strip to keep git diff --check happy
writeLines(cap, con)

cat("Wrote ", out_path, "\n", sep = "")
