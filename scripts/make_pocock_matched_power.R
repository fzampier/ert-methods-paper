# Pocock matched-power robustness check for the binary interruption-curve
# comparison. Companion to make_ertb_interruption_curve.R.
#
# Context: the primary comparison in make_ertb_interruption_curve.R uses
# N=2690 (calibrated for OBF 80% power) so that every method runs on the
# same simulated trial streams. At this N, Pocock is somewhat under-powered
# relative to its own matched-power calibration: K=3, two-sided 5%, 80%
# power requires Pocock inflation factor 1.166 over fixed-design vs OBF's
# 1.017 (Jennison & Turnbull 2000 Table 2.3), giving matched-power
# N ≈ 2690 × 1.166/1.017 ≈ 3084.
#
# This script re-runs Pocock and the two e-RTb policies at N=3084 to confirm
# that the qualitative conclusion (design e-RTb earlier than Pocock in the
# overwhelming majority of shared interruptions) holds at fair Pocock sizing.
# Writes pocock_matched_power_summary.csv for manuscript-numbers-check
# gating.

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))

root <- paper_root()
tables_dir <- paper_file("tables", root = root, must_work = FALSE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260511)

alpha <- 0.05
threshold <- 1 / alpha
n_sims <- 5000L
n_patients <- 3084L
p_ctrl <- 0.35
p_trt_alt <- 0.30
p_trt_null <- 0.35
look_fracs <- c(1 / 3, 2 / 3, 1)
look_n <- round(n_patients * look_fracs)
# Pocock K=3, two-sided 5%, equally-spaced looks: constant |Z| >= 2.289
# (Pocock 1977; Jennison & Turnbull 2000 Table 2.1).
pocock_bounds <- rep(2.289, 3)

simulate_trial <- function(n, p_trt, p_ctrl, p_random = 0.5) {
  treatment <- rbinom(n, 1, p_random)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, p_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, p_ctrl)
  list(treatment = treatment, outcome = outcome)
}

pocock_first_crossing <- function(treatment, outcome) {
  for (k in seq_along(look_n)) {
    idx <- seq_len(look_n[[k]])
    trt <- treatment[idx] == 1
    ctrl <- treatment[idx] == 0
    if (sum(trt) == 0 || sum(ctrl) == 0) next
    p_t <- mean(outcome[idx][trt])
    p_c <- mean(outcome[idx][ctrl])
    p_pool <- mean(outcome[idx])
    se <- sqrt(p_pool * (1 - p_pool) * (1 / sum(trt) + 1 / sum(ctrl)))
    if (!is.finite(se) || se <= 0) next
    z <- (p_c - p_t) / se
    if (is.finite(z) && abs(z) >= pocock_bounds[[k]]) {
      return(look_n[[k]])
    }
  }
  NA_integer_
}

ert_first_crossing <- function(wealth) {
  cross <- which(wealth >= threshold)
  if (!length(cross)) NA_integer_ else cross[[1]]
}

run_scenario <- function(p_trt, scenario, seed) {
  set.seed(seed)
  out <- data.frame(
    sim = seq_len(n_sims),
    pocock = NA_integer_,
    ertb_adaptive = NA_integer_,
    ertb_design = NA_integer_
  )
  for (sim in seq_len(n_sims)) {
    trial <- simulate_trial(n_patients, p_trt = p_trt, p_ctrl = p_ctrl)
    out$pocock[sim] <- pocock_first_crossing(trial$treatment, trial$outcome)
    w_adapt <- compute_eRT(
      trial$treatment, trial$outcome,
      burn_in = 50, ramp = 100, wager = "adaptive", kelly_fraction = 0.5
    )
    out$ertb_adaptive[sim] <- ert_first_crossing(w_adapt)
    w_des <- compute_eRT(
      trial$treatment, trial$outcome,
      wager = "design", p_trt_design = p_trt_alt, p_ctrl_design = p_ctrl,
      ramp_fixed = FALSE
    )
    out$ertb_design[sim] <- ert_first_crossing(w_des)
  }
  out$scenario <- scenario
  out
}

cat(sprintf("Running Pocock-N alternative (N=%d, looks at %s) ...\n",
            n_patients, paste(look_n, collapse = "/")))
alt <- run_scenario(p_trt_alt, "Alternative", seed = 20260544)
cat("Running Pocock-N null ...\n")
nul <- run_scenario(p_trt_null, "Null", seed = 20260545)

shared_earlier <- function(df, comparator, ert) {
  a <- df[[comparator]]
  b <- df[[ert]]
  both <- !is.na(a) & !is.na(b)
  if (!any(both)) return(NA_real_)
  mean(b[both] < a[both])
}

summary_row <- function(metric, value) data.frame(metric = metric, value = value)

summary_df <- rbind(
  summary_row("pocock_alt_interruption",            mean(!is.na(alt$pocock))),
  summary_row("pocock_null_interruption",           mean(!is.na(nul$pocock))),
  summary_row("ertb_design_alt_interruption",       mean(!is.na(alt$ertb_design))),
  summary_row("ertb_adaptive_alt_interruption",     mean(!is.na(alt$ertb_adaptive))),
  summary_row("ertb_design_earlier_than_pocock",    shared_earlier(alt, "pocock", "ertb_design")),
  summary_row("ertb_adaptive_earlier_than_pocock",  shared_earlier(alt, "pocock", "ertb_adaptive"))
)
summary_df$n_patients <- n_patients
summary_df$n_sims <- n_sims
summary_df <- summary_df[, c("metric", "value", "n_patients", "n_sims")]

csv_path <- file.path(tables_dir, "pocock_matched_power_summary.csv")
write.csv(summary_df, csv_path, row.names = FALSE)

cat("\n=== Pocock matched-power summary (N=", n_patients, ", ", n_sims,
    " sims) ===\n", sep = "")
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %-40s %.4f\n", summary_df$metric[i], summary_df$value[i]))
}
cat("\nWrote ", csv_path, "\n", sep = "")
