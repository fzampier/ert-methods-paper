# Round 2, Reviewer 1 point 2: where does the ~11-point e-RTb design premium come from?
#
# Decomposes, on the FROZEN comparator streams (alt seed 20260542, null seed 20260543;
# regenerated here exactly, kernels consume no RNG), the gap between a fixed-N one-sided
# z-test and the anytime e-RTb design crossing into
#   (1) representation: fixed-N z-test  ->  fixed-N test on the FINAL design-wager wealth
#       (the randomization-model likelihood ratio), calibrated to exact level alpha on an
#       independent null sample (fresh seed 20260554, closed-form cell-count draws);
#   (2) boundary/horizon: fixed-N final-wealth test -> sup_i W_i >= 1/alpha (Ville).
# Reproduction gate: the regenerated streams must reproduce the pinned design rows
# (alt 0.7486 / E[min] 1693.04 / median 1282; null 0.0406 / E[min] 2636.962) or the
# script stops. Run from the repository root.

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))
root <- paper_root()
out_dir <- paper_file("tables", root = root, must_work = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05; threshold <- 20; n_sims <- 5000L; N <- 2690L
p_ctrl <- 0.35; p_trt_alt <- 0.30; p_rand <- 0.5

simulate_trial <- function(n, p_trt, p_ctrl, p_random = 0.5) {   # verbatim from the frozen script
  treatment <- rbinom(n, 1, p_random)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, p_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, p_ctrl)
  list(treatment = treatment, outcome = outcome)
}
z_final <- function(treatment, outcome) {                         # z_at_look at n_look = N
  trt <- treatment == 1; ctrl <- !trt
  p_pool <- mean(outcome)
  se <- sqrt(p_pool * (1 - p_pool) * (1 / sum(trt) + 1 / sum(ctrl)))
  (mean(outcome[ctrl]) - mean(outcome[trt])) / se
}
lam <- binary_design_lambdas(p_trt_alt, p_ctrl, p_rand)
a1 <- log(lam$event / p_rand);    b1 <- log((1 - lam$event) / (1 - p_rand))
a0 <- log(lam$nonevent / p_rand); b0 <- log((1 - lam$nonevent) / (1 - p_rand))
closed_logW <- function(t, y) sum(t * y) * a1 + sum((1 - t) * y) * b1 +
  sum(t * (1 - y)) * a0 + sum((1 - t) * (1 - y)) * b0

run_streams <- function(p_trt, seed) {
  set.seed(seed)
  out <- data.frame(crossed = logical(n_sims), first = NA_integer_, logWN = NA_real_,
                    z = NA_real_, cf_all = NA_real_, cf_skip1 = NA_real_)
  for (s in seq_len(n_sims)) {
    tr <- simulate_trial(N, p_trt = p_trt, p_ctrl = p_ctrl)
    res <- compute_eRT(tr$treatment, tr$outcome, wager = "design",
                       p_trt_design = p_trt_alt, p_ctrl_design = p_ctrl, ramp_fixed = FALSE)
    path <- if (is.list(res)) res$wealth else res
    hit <- which(path >= threshold)
    out$crossed[s] <- length(hit) > 0
    out$first[s] <- if (length(hit)) hit[1] else NA_integer_
    out$logWN[s] <- log(path[N])
    out$z[s] <- z_final(tr$treatment, tr$outcome)
    out$cf_all[s] <- closed_logW(tr$treatment, tr$outcome)
    out$cf_skip1[s] <- closed_logW(tr$treatment[-1], tr$outcome[-1])
  }
  out
}
cat("Regenerating frozen streams (2 x 5000 trials)...\n")
alt <- run_streams(p_trt_alt, 20260542); nul <- run_streams(p_ctrl, 20260543)

# ---- reproduction gate -------------------------------------------------------------
emin <- function(d) mean(ifelse(is.na(d$first), N, pmin(d$first, N)))
chk <- c(alt_cross = mean(alt$crossed), alt_emin = emin(alt), alt_med = median(alt$first, na.rm = TRUE),
         nul_cross = mean(nul$crossed), nul_emin = emin(nul))
pin <- c(alt_cross = 0.7486, alt_emin = 1693.04, alt_med = 1282, nul_cross = 0.0406, nul_emin = 2636.962)
print(rbind(regenerated = chk, pinned = pin))
stopifnot(all(abs(chk - pin) < c(1e-9, 5e-3, 1e-9, 1e-9, 5e-4)))
cat("Reproduction gate PASSED: regenerated streams are the frozen streams.\n")

# ---- closed form vs kernel ----------------------------------------------------------
d_all <- max(abs(c(alt$cf_all, nul$cf_all) - c(alt$logWN, nul$logWN)))
d_sk1 <- max(abs(c(alt$cf_skip1, nul$cf_skip1) - c(alt$logWN, nul$logWN)))
cat(sprintf("closed-form check: max|diff| all patients %.3g ; skipping patient 1 %.3g\n", d_all, d_sk1))
n_bet <- if (d_all < 1e-8) N else if (d_sk1 < 1e-8) N - 1L else stop("closed form matches neither")
cat("kernel bets on", n_bet, "patients\n")

# ---- exact-level calibration of the fixed-N final-wealth test -----------------------
set.seed(20260554); n_cal <- 1e6L
n1 <- rbinom(n_cal, n_bet, p_ctrl)                 # events among betting patients (null)
nT1 <- rbinom(n_cal, n1, p_rand); nT0 <- rbinom(n_cal, n_bet - n1, p_rand)
logW_cal <- nT1 * a1 + (n1 - nT1) * b1 + nT0 * a0 + (n_bet - n1 - nT0) * b0
c_log <- sort(logW_cal)[ceiling((1 - alpha) * n_cal) + 1L]          # smallest c with P0(logW >= c) <= alpha
size_cal <- mean(logW_cal >= c_log)
cat(sprintf("critical value: W_N >= %.4f (log %.4f); calibration size %.5f; mean null W_N %.4f\n",
            exp(c_log), c_log, size_cal, mean(exp(logW_cal))))

# ---- the chain, all on the same frozen streams ---------------------------------------
zc <- qnorm(1 - alpha)
ind <- function(d) data.frame(z = d$z >= zc, lr = d$logWN >= c_log, sup20 = d$crossed,
                              term20 = d$logWN >= log(threshold))
A <- ind(alt); Z <- ind(nul)
rate <- function(x) c(est = mean(x), se = sqrt(mean(x) * (1 - mean(x)) / length(x)))
pdiff <- function(x, y) { d <- x - y; c(est = mean(d), se = sd(d) / sqrt(length(d))) }
tab <- rbind(
  "fixed-N one-sided z-test (alpha 0.05)"              = c(rate(A$z),      null = mean(Z$z)),
  "fixed-N final-wealth test, exact level (W_N >= c)"  = c(rate(A$lr),     null = mean(Z$lr)),
  "anytime crossing sup W >= 20 (e-RTb design)"        = c(rate(A$sup20),  null = mean(Z$sup20)),
  "final wealth W_N >= 20 (no sup; reference)"         = c(rate(A$term20), null = mean(Z$term20)))
gaps <- rbind("representation: z-test minus final-wealth test" = pdiff(A$z, A$lr),
              "boundary/horizon: final-wealth test minus anytime crossing" = pdiff(A$lr, A$sup20),
              "total: z-test minus anytime crossing" = pdiff(A$z, A$sup20))
cat("\n=== power (alt) with MC SE, and null rejection, 5000 frozen trials each ===\n"); print(round(tab, 4))
cat("\n=== paired gaps on the alternative streams (percentage points = x100) ===\n"); print(round(gaps, 4))

res <- data.frame(quantity = c(rownames(tab), rownames(gaps)),
                  estimate = c(tab[, "est"], gaps[, "est"]), mc_se = c(tab[, "se"], gaps[, "se"]),
                  null_rate = c(tab[, "null"], rep(NA, nrow(gaps))))
write.csv(res, file.path(out_dir, "premium_decomposition.csv"), row.names = FALSE)
writeLines(c(sprintf("critical_value_W,%.6f", exp(c_log)), sprintf("calibration_size,%.6f", size_cal),
             sprintf("n_bet,%d", n_bet), sprintf("n_cal,%d", n_cal), "calibration_seed,20260554",
             "stream_seeds,20260542 (alt) / 20260543 (null)"),
           file.path(out_dir, "premium_decomposition_meta.csv"))
cat("\nwritten:", file.path(out_dir, "premium_decomposition.csv"), "\n")
