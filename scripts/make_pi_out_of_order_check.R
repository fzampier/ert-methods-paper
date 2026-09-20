# Round 2, Reviewer 4 comment 1: which information set must p_i condition on when updates
# run in OUTCOME-CONFIRMATION order instead of enrollment order?
#
# Null (35% events in both arms), N = 2,690, permuted blocks of 4. Enrollment 5 patients/day.
# Events are confirmed Uniform(1, 28) days after enrollment, non-events at day 28, ties by
# enrollment index: the update order depends on outcomes and timing, never on assignments.
# About 140 patients are in flight at any time, so blocks interleave heavily.
#
# Bettors, all monitored at W >= 20 (mirrors make_blocked_randomization_sensitivity.R):
#   revealed_count      adaptive e-RTb wager priced at the conditional probability given the
#                       assignments ALREADY REVEALED TO THE MONITOR in the patient's block
#                       (the rule the article states); forced slots skipped.
#   logged_design       same wager priced at the design's own enrollment-order probability
#                       (what the randomization system would log); forced slots skipped.
#   adversary_vs_logged history-only bettor: bets the revealed-count probability against
#                       the logged design price. Uses no outcome. Invalidity demonstration.
#   (An adversary against the revealed-count price bets the price itself: W = 1 forever.)
# Fresh seed 20260555. Writes tables/pi_out_of_order_check.csv.

args <- commandArgs(trailingOnly = FALSE)
script <- gsub("~+~", " ", sub("--file=", "", args[grep("--file=", args)]), fixed = TRUE)
out_dir <- file.path(normalizePath(file.path(dirname(script), "..")), "tables")
dir.create(out_dir, showWarnings = FALSE)
set.seed(20260555)
n_sims <- as.integer(Sys.getenv("PISIM_NSIMS", "5000"))
N <- 2690L; B <- 4L; p_event <- 0.35; threshold <- 20; rate <- 5; burn_in <- 50; ramp <- 100; kf <- 0.5
clamp <- function(x) max(0.001, min(0.999, x))

one_trial <- function() {
  nb <- ceiling(N / B)
  trt_e <- unlist(lapply(seq_len(nb), function(b) sample(rep(c(0, 1), B / 2))))[seq_len(N)]
  y_e <- rbinom(N, 1, p_event)
  t_conf <- seq_len(N) / rate + ifelse(y_e == 1, runif(N, 1, 28), 28)
  ord <- order(t_conf, seq_len(N))                      # processing order (no assignment input)
  # design's own (enrollment-order) probability, as the randomization system would log it
  p_design_e <- numeric(N)
  for (e in seq_len(N)) { pos <- (e - 1L) %% B
    p_design_e[e] <- if (pos == 0L) 0.5 else (B / 2 - sum(trt_e[(e - pos):(e - 1L)])) / (B - pos) }
  blk <- ceiling(ord / B); trt <- trt_e[ord]; y <- y_e[ord]; p_design <- p_design_e[ord]
  rev_t <- integer(nb); rev_n <- integer(nb)
  w <- c(revealed_count = 1, logged_design = 1, adversary_vs_logged = 1)
  fc <- c(revealed_count = NA_integer_, logged_design = NA_integer_, adversary_vs_logged = NA_integer_)
  n_t <- 0L; e_t <- 0L; n_c <- 0L; e_c <- 0L; forced_rev <- 0L; forced_des <- 0L; differ <- 0L
  for (i in seq_len(N)) {
    b <- blk[i]
    p_rev <- (B / 2 - rev_t[b]) / (B - rev_n[b])        # given what the monitor has seen
    # adaptive wager from the processed history (running sums = the kernel's running means)
    lam0 <- 0.5
    if (i > burn_in && n_t > 0L && n_c > 0L) {
      dhat <- e_t / n_t - e_c / n_c; ci <- min(1, (i - burn_in) / ramp)
      lam0 <- clamp(if (y[i] == 1) 0.5 + kf * ci * dhat else 0.5 - kf * ci * dhat)
    }
    mult <- function(lam, p) if (trt[i] == 1) lam / p else (1 - lam) / (1 - p)
    if (p_rev > 0 && p_rev < 1) w[1] <- w[1] * mult(clamp(p_rev + (lam0 - 0.5)), p_rev) else forced_rev <- forced_rev + 1L
    if (p_design[i] > 0 && p_design[i] < 1) {
      w[2] <- w[2] * mult(clamp(p_design[i] + (lam0 - 0.5)), p_design[i])
      w[3] <- w[3] * mult(clamp(p_rev), p_design[i])
    } else forced_des <- forced_des + 1L
    if (abs(p_rev - p_design[i]) > 1e-12) differ <- differ + 1L
    for (k in 1:3) if (is.na(fc[k]) && w[k] >= threshold) fc[k] <- i
    rev_t[b] <- rev_t[b] + trt[i]; rev_n[b] <- rev_n[b] + 1L
    if (trt[i] == 1) { n_t <- n_t + 1L; e_t <- e_t + y[i] } else { n_c <- n_c + 1L; e_c <- e_c + y[i] }
  }
  c(fc, forced_rev = forced_rev, forced_des = forced_des, differ = differ,
    displaced = mean(abs(ord - seq_len(N))))
}
t0 <- Sys.time(); R <- t(replicate(n_sims, one_trial())); cat("elapsed:", format(Sys.time() - t0), "\n")
summ <- do.call(rbind, lapply(c("revealed_count", "logged_design", "adversary_vs_logged"), function(m) {
  x <- !is.na(R[, m]); data.frame(bettor = m, null_crossing = mean(x), mc_se = sqrt(mean(x) * (1 - mean(x)) / n_sims),
                                 median_first_crossing = if (any(x)) median(R[x, m]) else NA) }))
print(summ, digits = 4)
cat(sprintf("\nper trial: forced slots under revealed-count pricing %.0f, under logged-design pricing %.0f;\nupdates where the two prices differ %.0f of %d; mean displacement from enrollment order %.0f patients\n",
            mean(R[, "forced_rev"]), mean(R[, "forced_des"]), mean(R[, "differ"]), N, mean(R[, "displaced"])))
write.csv(summ, file.path(out_dir, "pi_out_of_order_check.csv"), row.names = FALSE)
