# Empirical Type I error for e-RTc under several null data-generating processes.
#
# The simulations cover Type I error exhaustively for e-RTb. This script does
# the same for e-RTc, which the FACTT illustration (98.9% permutation crossings
# on fluid VFD-28) uses with the Cohen's-d wager.
#
# DGPs cover the regimes relevant to the manuscript:
#   - Normal(0,1)             : canonical continuous baseline
#   - t_3                     : heavy-tailed robustness check
#   - Lognormal(0, 1)         : skewed robustness check
#   - Zero-inflated integer   : integer outcome with ~20% zero mass and a
#                               truncated positive integer tail, the shape
#                               of days-free-of-event count endpoints that
#                               motivate the Cohen's-d wager
#
# Settings cover the manuscript-stated default and a more conservative
# variant suited to discrete / zero-inflated endpoints:
#   - default      : burn_in=20, ramp=50, c_max=0.6, adaptive_direction="sign"
#   - conservative : burn_in=50, ramp=100, c_max=0.6, adaptive_direction="cohens_d"
#
# Under the null, treatment assignment is a fair coin and outcome is drawn
# from the same marginal regardless of arm, so by the Section 2.1 validity
# argument the wealth process is a nonnegative martingale and Ville's
# inequality bounds the crossing rate at 1/threshold = 0.05.

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertc.R")))

root <- paper_root()
tables_dir <- paper_file("tables", root = root, must_work = FALSE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260518)

N_PATIENTS <- 500L
N_SIMS <- 5000L
THRESHOLD <- 20

dgp_normal <- function(n) rnorm(n, mean = 0, sd = 1)
dgp_t3 <- function(n) rt(n, df = 3)
dgp_lognormal <- function(n) rlnorm(n, meanlog = 0, sdlog = 1)
dgp_vfd_like <- function(n) {
  zero_mask <- rbinom(n, 1, 0.20) == 1L
  out <- integer(n)
  out[!zero_mask] <- pmin(28L, rpois(sum(!zero_mask), lambda = 15))
  out
}

dgps <- list(
  Normal = dgp_normal,
  `Student's t (df=3)` = dgp_t3,
  Lognormal = dgp_lognormal,
  `Zero-inflated integer` = dgp_vfd_like
)

settings <- list(
  list(label = "Default (sign; burn 20, ramp 50)",
       burn_in = 20L, ramp = 50L, c_max = 0.6, adaptive_direction = "sign"),
  list(label = "Conservative (Cohen's d; burn 50, ramp 100)",
       burn_in = 50L, ramp = 100L, c_max = 0.6, adaptive_direction = "cohens_d")
)

results <- list()
row_i <- 1L

t_start <- Sys.time()
for (dgp_name in names(dgps)) {
  dgp <- dgps[[dgp_name]]
  for (s in settings) {
    crossed <- logical(N_SIMS)
    final_w <- numeric(N_SIMS)
    peak_w <- numeric(N_SIMS)

    for (sim in seq_len(N_SIMS)) {
      treatment <- rbinom(N_PATIENTS, 1, 0.5)
      outcome <- dgp(N_PATIENTS)
      wealth <- compute_eRTc(
        treatment = treatment,
        outcome = outcome,
        burn_in = s$burn_in,
        ramp = s$ramp,
        c_max = s$c_max,
        wager = "adaptive",
        adaptive_direction = s$adaptive_direction
      )
      crossed[sim] <- any(wealth >= THRESHOLD)
      final_w[sim] <- wealth[N_PATIENTS]
      peak_w[sim] <- max(wealth)
    }

    rate <- mean(crossed)
    se <- sqrt(rate * (1 - rate) / N_SIMS)

    cat(sprintf("  %-30s | %-26s | type1=%.3f (SE=%.3f) | med peak=%.3f | q95 peak=%.3f\n",
                dgp_name, s$label, rate, se,
                median(peak_w), quantile(peak_w, 0.95)))

    results[[row_i]] <- data.frame(
      dgp = dgp_name,
      setting = s$label,
      n_patients = N_PATIENTS,
      n_sims = N_SIMS,
      type1 = rate,
      se = se,
      median_peak_wealth = median(peak_w),
      q95_peak_wealth = quantile(peak_w, 0.95),
      median_final_wealth = median(final_w),
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
}

# ---------------------------------------------------------------------------
# Freeze amendment 2 (R4.m9, BMC revision): c_max sensitivity rows.
# Two additional settings at the default policy (sign; burn 20, ramp 50) with
# c_max = 0.4 and 0.8, bracketing the shipped default 0.6. APPENDED AFTER the
# original 4x2 grid under a FRESH seed -- the original rows above run first on
# the untouched RNG stream and must reproduce exactly.
# ---------------------------------------------------------------------------

set.seed(20260547)

cmax_settings <- list(
  list(label = "Default policy, c_max=0.4 (sign; burn 20, ramp 50)",
       burn_in = 20L, ramp = 50L, c_max = 0.4, adaptive_direction = "sign"),
  list(label = "Default policy, c_max=0.8 (sign; burn 20, ramp 50)",
       burn_in = 20L, ramp = 50L, c_max = 0.8, adaptive_direction = "sign")
)

for (dgp_name in names(dgps)) {
  dgp <- dgps[[dgp_name]]
  for (s in cmax_settings) {
    crossed <- logical(N_SIMS)
    final_w <- numeric(N_SIMS)
    peak_w <- numeric(N_SIMS)

    for (sim in seq_len(N_SIMS)) {
      treatment <- rbinom(N_PATIENTS, 1, 0.5)
      outcome <- dgp(N_PATIENTS)
      wealth <- compute_eRTc(
        treatment = treatment,
        outcome = outcome,
        burn_in = s$burn_in,
        ramp = s$ramp,
        c_max = s$c_max,
        wager = "adaptive",
        adaptive_direction = s$adaptive_direction
      )
      crossed[sim] <- any(wealth >= THRESHOLD)
      final_w[sim] <- wealth[N_PATIENTS]
      peak_w[sim] <- max(wealth)
    }

    rate <- mean(crossed)
    se <- sqrt(rate * (1 - rate) / N_SIMS)

    cat(sprintf("  %-30s | %-26s | type1=%.3f (SE=%.3f) | med peak=%.3f | q95 peak=%.3f\n",
                dgp_name, s$label, rate, se,
                median(peak_w), quantile(peak_w, 0.95)))

    results[[row_i]] <- data.frame(
      dgp = dgp_name,
      setting = s$label,
      n_patients = N_PATIENTS,
      n_sims = N_SIMS,
      type1 = rate,
      se = se,
      median_peak_wealth = median(peak_w),
      q95_peak_wealth = quantile(peak_w, 0.95),
      median_final_wealth = median(final_w),
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
}

out <- do.call(rbind, results)
rownames(out) <- NULL

csv_path <- file.path(tables_dir, "ertc_type1_simulation.csv")
md_path <- file.path(tables_dir, "ertc_type1_simulation.md")

write.csv(out, csv_path, row.names = FALSE)

fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
fmt_num <- function(x, d = 2) sprintf(paste0("%.", d, "f"), x)

md <- data.frame(
  DGP = out$dgp,
  Setting = out$setting,
  N = out$n_patients,
  `n_sims` = out$n_sims,
  `Empirical Type I` = paste0(fmt_pct(out$type1), " (SE ", fmt_pct(out$se, 2), ")"),
  `Median peak W` = fmt_num(out$median_peak_wealth),
  `Q95 peak W` = fmt_num(out$q95_peak_wealth),
  check.names = FALSE
)

header <- paste0("| ", paste(names(md), collapse = " | "), " |")
sep <- paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |")
rows <- apply(md, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
writeLines(c(header, sep, rows), md_path)

cat("\nWrote ", csv_path, "\n")
cat("Wrote ", md_path, "\n")
cat("Total runtime: ", format(Sys.time() - t_start), "\n", sep = "")
