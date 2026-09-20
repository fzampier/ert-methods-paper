# Kernel identity check.
#
# The three kernels in R/ are the code behind every result in this repository.
# They must match the supplementary files of the article (supplementary_code/,
# also printed in ESM Sections 1.5, 2.4, and 3.5): byte for byte for e-RTe and
# e-RTc, and for e-RTb apart from the guard on its closing demo block (the
# supplementary file runs that block in any interactive session).
#
# Run via `make check` or `Rscript R/check_kernels.R`.

source("scripts/paths.R")
root <- paper_root()
failures <- 0L

report <- function(label, cond) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) failures <<- failures + 1L
}

read_bytes <- function(...) {
  path <- file.path(root, ...)
  readBin(path, what = "raw", n = file.size(path))
}

for (k in c("e", "c")) {
  report(sprintf("R/ert%s.R is byte-identical to supplementary_code/SuppCode_eRT%s.R", k, k),
         identical(read_bytes("R", sprintf("ert%s.R", k)),
                   read_bytes("supplementary_code", sprintf("SuppCode_eRT%s.R", k))))
}

ours <- readLines(file.path(root, "R", "ertb.R"), warn = FALSE)
theirs <- readLines(file.path(root, "supplementary_code", "SuppCode_eRTb.R"), warn = FALSE)
at <- which(theirs == "if (interactive()) {")
report("SuppCode_eRTb.R has one demo block", length(at) == 1L)
if (length(at) == 1L) {
  guard <- ours[at:(at + 2L)]
  report("R/ertb.R differs from SuppCode_eRTb.R only by the demo guard",
         identical(ours[-(at:(at + 2L))], theirs[-at]) &&
           all(startsWith(guard[1:2], "#")) &&
           identical(guard[3], "if (interactive() && isTRUE(getOption(\"ert.run_demo\", FALSE))) {"))
}

cat(sprintf("\n=== kernels: %d failure(s) ===\n", failures))
quit(status = if (failures) 1L else 0L)
