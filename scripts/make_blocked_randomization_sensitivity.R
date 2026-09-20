# T9: Blocked-randomization sensitivity (A6) -- answers R4.7 / R1.10.
#
# Under permuted-block randomization the conditional assignment probability
# p_i is NOT constant at 0.5: within a block of size 4, once one arm's slots
# fill, the remaining assignments are increasingly predictable (and the last
# slot is deterministic).  Three bettors, all monitored at W >= 20, all under
# the NULL (no treatment effect; outcomes pure noise):
#
#   (1) history-only bettor, naive p = 0.5:  bets on the under-allocated arm
#       using ONLY the assignment history (no outcome information at all).
#       With naive p = 0.5 in the wealth update this is a money pump: wealth
#       grows deterministically inside every block tail.  Pure invalidity
#       demonstration -- the strongest form of R4.7's "spurious wealth growth".
#   (2) practical outcome wagers (the manuscript's adaptive e-RTb), naive
#       p = 0.5: quantifies the real-world exposure of the shipped method if
#       someone runs it unmodified under permuted blocks.
#   (3) correct conditional p_i (block-aware) with deterministic slots skipped
#       (p_i in {0,1} -> no update), same practical wagers: validity restored.
#
# Self-contained wealth loop (the shipped kernels take scalar p by design and
# are not modified mid-review).  Fresh seed; no shared streams with the
# fair-comparator run.
#
# Output: tables/ertb_blocked_randomization_sensitivity.md / .csv

source("scripts/paths.R")

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

set.seed(20260545)
n_sims <- as.integer(Sys.getenv("BLOCKSIM_NSIMS", "5000"))
n_patients <- 2690
block_size <- 4
p_event <- 0.35          # null: same event rate in both arms
threshold <- 20

permuted_block_assign <- function(n, block = 4) {
  n_blocks <- ceiling(n / block)
  out <- unlist(lapply(seq_len(n_blocks), function(b) sample(rep(c(0, 1), block / 2))))
  out[seq_len(n)]
}

# conditional P(T_i = 1) inside permuted blocks given assignments so far
block_conditional_p <- function(treatment, block = 4) {
  n <- length(treatment)
  p <- numeric(n)
  for (i in seq_len(n)) {
    pos <- (i - 1) %% block          # 0-based position within the current block
    if (pos == 0) {
      p[i] <- 0.5
    } else {
      start <- i - pos
      t_used <- sum(treatment[start:(i - 1)])
      remaining_t <- block / 2 - t_used
      remaining_total <- block - pos
      p[i] <- remaining_t / remaining_total
    }
  }
  p
}

# adaptive e-RTb wager (mirrors the shipped kernel's defaults: burn-in 50,
# ramp 100, half-Kelly on the running risk difference)
adaptive_lambda <- function(i, treatment, outcome, burn_in = 50, ramp = 100,
                            kelly_fraction = 0.5) {
  if (i <= burn_in) return(0.5)
  hist_idx <- seq_len(i - 1)
  trt <- treatment[hist_idx] == 1
  ctrl <- !trt
  if (!any(trt) || !any(ctrl)) return(0.5)
  delta_hat <- mean(outcome[hist_idx][trt]) - mean(outcome[hist_idx][ctrl])
  c_i <- min(1, (i - burn_in) / ramp)
  if (outcome[i] == 1) {
    lam <- 0.5 + kelly_fraction * c_i * delta_hat
  } else {
    lam <- 0.5 - kelly_fraction * c_i * delta_hat
  }
  max(0.001, min(0.999, lam))
}

run_bettor <- function(treatment, outcome, mode) {
  n <- length(treatment)
  p_cond <- block_conditional_p(treatment, block_size)
  w <- 1
  w_max <- 1
  fc <- NA_integer_
  for (i in seq_len(n)) {
    if (mode == "history_naive") {
      # bet on the under-allocated arm within the block; naive p = 0.5
      lam <- max(0.001, min(0.999, p_cond[i]))
      p_used <- 0.5
    } else if (mode == "practical_naive") {
      lam <- adaptive_lambda(i, treatment, outcome)
      p_used <- 0.5
    } else {  # correct_conditional
      if (p_cond[i] %in% c(0, 1)) next   # deterministic slot: skip, no update
      lam <- adaptive_lambda(i, treatment, outcome)
      # recenter the practical wager on the true conditional null probability
      lam <- max(0.001, min(0.999, p_cond[i] + (lam - 0.5)))
      p_used <- p_cond[i]
    }
    if (treatment[i] == 1) {
      w <- w * lam / p_used
    } else {
      w <- w * (1 - lam) / (1 - p_used)
    }
    if (w > w_max) w_max <- w
    if (is.na(fc) && w >= threshold) fc <- i
  }
  list(first_crossing = fc, peak = w_max, final = w)
}

modes <- c("history_naive", "practical_naive", "correct_conditional")
res <- vector("list", n_sims * length(modes))
ri <- 1L
t0 <- Sys.time()
for (sim in seq_len(n_sims)) {
  treatment <- permuted_block_assign(n_patients, block_size)
  outcome <- rbinom(n_patients, 1, p_event)   # null: outcome independent of arm
  for (m in modes) {
    r <- run_bettor(treatment, outcome, m)
    res[[ri]] <- data.frame(sim = sim, mode = m,
                            first_crossing = r$first_crossing,
                            peak = r$peak, final = r$final)
    ri <- ri + 1L
  }
}
res <- do.call(rbind, res)

summ <- do.call(rbind, lapply(split(res, res$mode), function(d) {
  crossed <- mean(!is.na(d$first_crossing))
  data.frame(
    mode = d$mode[1],
    crossed = crossed,
    se = sqrt(crossed * (1 - crossed) / nrow(d)),
    median_crossing = median(d$first_crossing, na.rm = TRUE),
    median_peak = median(d$peak),
    q90_peak = quantile(d$peak, 0.9, names = FALSE)
  )
}))

label <- c(
  history_naive = "History-only bettor, naive p=0.5 (invalidity demonstration)",
  practical_naive = "Adaptive e-RTb wagers, naive p=0.5 (practical exposure)",
  correct_conditional = "Adaptive e-RTb wagers, correct conditional p_i + skip rule"
)
summ$Bettor <- label[summ$mode]

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
md <- data.frame(
  Bettor = summ$Bettor,
  `Null crossing rate (W>=20)` = sprintf("%s (MC SE %.2fpp)", fmt_pct(summ$crossed), 100 * summ$se),
  `Median first crossing` = ifelse(is.na(summ$median_crossing), "--",
                                   format(round(summ$median_crossing), big.mark = ",")),
  `Median peak W` = ifelse(summ$median_peak >= 1e6, sprintf("%.2e", summ$median_peak),
                           sprintf("%.2f", summ$median_peak)),
  `90th pct peak W` = ifelse(summ$q90_peak >= 1e6, sprintf("%.2e", summ$q90_peak),
                             sprintf("%.2f", summ$q90_peak)),
  check.names = FALSE
)

write.csv(res, file.path(table_dir, "ertb_blocked_randomization_sensitivity.csv"), row.names = FALSE)
writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  sprintf("Permuted blocks of %d, N=%d, null event rate %.0f%%, %s simulated trials, seed 20260545.",
          block_size, n_patients, 100 * p_event, format(n_sims, big.mark = ","))
), file.path(table_dir, "ertb_blocked_randomization_sensitivity.md"))

cat(sprintf("Done in %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(md, row.names = FALSE)
