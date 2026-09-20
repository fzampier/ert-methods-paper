# Freeze amendment 6 (author request 2026-08-10): e-RTc detection power under
# normal-shift alternatives, closing the ESM Section 3 gap (which was
# null-only).
#
# UNIFORMITY DIRECTIVE (author): power demonstrations are standalone and
# uniform across flavours. Each trial is the size a traditional two-sample
# t-test would use for the same shift, two-sided alpha 0.05 (power.t.test),
# at 80% and 90% target power:
#   delta 0.2 SD: N = 788 (80%), 1,054 (90%)
#   delta 0.3 SD: N = 352 (80%),   470 (90%)
#
# Scenarios: control N(0,1), treatment N(delta,1); 1:1 simple randomization;
# 5,000 trials per cell. Bettors:
#   - adaptive default      (sign direction; burn 20, ramp 50, c_max 0.6)
#   - adaptive conservative (Cohen's d;      burn 50, ramp 100, c_max 0.6)
#   - design matched        (normal-shift model mu_ctrl 0, mu_trt delta, sd 1)
#
# NEW standalone script under FRESH seed 20260551; no existing random stream
# is touched, so all frozen tables reproduce unchanged.
#
# Output: tables/ertc_power_check.{csv,md}

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertc.R")))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

set.seed(20260551)
N_SIMS <- as.integer(Sys.getenv("ERTCPOW_NSIMS", "5000"))
THRESHOLD <- 20

grid <- list()
gi <- 1L
for (delta in c(0.2, 0.3)) {
  for (tp in c(0.80, 0.90)) {
    n_arm <- ceiling(power.t.test(delta = delta, sd = 1, power = tp,
                                  sig.level = 0.05)$n)
    grid[[gi]] <- list(delta = delta, target_power = tp, N = 2L * as.integer(n_arm))
    gi <- gi + 1L
  }
}

bettors <- list(
  list(label = "Adaptive default (sign; burn 20, ramp 50)",
       fn = function(tr, y, delta) compute_eRTc(tr, y, burn_in = 20L, ramp = 50L,
                                                c_max = 0.6, wager = "adaptive",
                                                adaptive_direction = "sign")),
  list(label = "Adaptive conservative (Cohen's d; burn 50, ramp 100)",
       fn = function(tr, y, delta) compute_eRTc(tr, y, burn_in = 50L, ramp = 100L,
                                                c_max = 0.6, wager = "adaptive",
                                                adaptive_direction = "cohens_d")),
  list(label = "Design matched (normal shift)",
       fn = function(tr, y, delta) compute_eRTc(tr, y, wager = "design",
                                                mu_ctrl_design = 0,
                                                mu_trt_design = delta,
                                                sd_design = 1))
)

first_cross <- function(w) { i <- which(w >= THRESHOLD); if (length(i)) i[1] else NA_integer_ }

results <- list()
ri <- 1L
t0 <- Sys.time()
for (g in grid) {
  for (b in bettors) {
    crossed <- logical(N_SIMS)
    fc <- rep(NA_integer_, N_SIMS)
    for (sim in seq_len(N_SIMS)) {
      treatment <- rbinom(g$N, 1, 0.5)
      outcome <- rnorm(g$N, mean = g$delta * treatment, sd = 1)
      wealth <- b$fn(treatment, outcome, g$delta)
      crossed[sim] <- any(wealth >= THRESHOLD)
      fc[sim] <- first_cross(wealth)
    }
    r <- mean(crossed)
    cat(sprintf("  delta=%.1f %d%% N=%d | %-50s | power=%.3f\n",
                g$delta, round(100 * g$target_power), g$N, b$label, r))
    results[[ri]] <- data.frame(
      delta = g$delta,
      target_power = g$target_power,
      N = g$N,
      bettor = b$label,
      power = r,
      se = sqrt(r * (1 - r) / N_SIMS),
      median_first_crossing_patient = median(fc, na.rm = TRUE),
      pct_of_N = median(fc, na.rm = TRUE) / g$N,
      n_sims = N_SIMS
    )
    ri <- ri + 1L
  }
}
out <- do.call(rbind, results)
write.csv(out, file.path(table_dir, "ertc_power_check.csv"), row.names = FALSE)

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
md <- data.frame(
  `Shift (SD)` = sprintf("%.1f", out$delta),
  `Target power` = sprintf("%d%%", round(100 * out$target_power)),
  `N (t-test)` = format(out$N, big.mark = ","),
  Bettor = out$bettor,
  `Power (W>=20)` = sprintf("%s (%.2fpp)", fmt_pct(out$power), 100 * out$se),
  `Median crossing (patient)` = ifelse(is.na(out$median_first_crossing_patient), "--",
                                       sprintf("%d", as.integer(out$median_first_crossing_patient))),
  `% N` = fmt_pct(out$pct_of_N),
  check.names = FALSE
)
writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  sprintf("Control N(0,1), treatment N(delta,1). Each trial uses the traditional two-sample t-test sample size for its shift and target power, two-sided alpha 0.05 (power.t.test). %s simulated trials per cell, 1:1 simple randomization, seed 20260551. Median crossing among crossed trials.",
          format(N_SIMS, big.mark = ","))
), file.path(table_dir, "ertc_power_check.md"))

cat(sprintf("Done in %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(md, row.names = FALSE)
