# T6: The price of anytime validity (A9) -- analytic power(K) plateau.
#
# For the manuscript's binary scenario (N=2690, control 35% vs treatment 30%,
# 1:1), compute the power of exactly calibrated Lan-DeMets OBF-type
# error-spending group-sequential designs as the number of equally spaced
# looks K grows (K = 1, 2, 3, 5, 10, 20), one-sided alpha=0.05 (D1 primary
# convention; the comparator for the directional design wager) and two-sided
# alpha=0.05 (comparator for adaptive/mixture e-RT).  Under the canonical
# group-sequential Brownian approximation, Z_k ~ N(theta*sqrt(t_k), 1) with
# cov(Z_j, Z_k) = sqrt(t_j/t_k); power and expected stopping information have
# closed multivariate-normal forms (mvtnorm), so these numbers carry no Monte
# Carlo error.  This quantifies R3.1's K=5-suffices point in the reviewer's
# own currency and, next to the e-RT crossing rates from the fair-comparator
# simulation, states the price of anytime validity as a number (R1.3).
#
# Outputs: tables/power_vs_K.csv, tables/power_vs_K.md
# Figure (W4) is generated separately once the empirical overlay is frozen.

source("scripts/paths.R")
suppressPackageStartupMessages(library(mvtnorm))
stopifnot(requireNamespace("rpact", quietly = TRUE))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

n_patients <- 2690
p_ctrl <- 0.35
p_trt <- 0.30

# drift of the Z-statistic at full information (pooled-variance convention,
# matching z_at_look in the simulation scripts)
n_arm <- n_patients / 2
p_bar <- (p_ctrl + p_trt) / 2
se_full <- sqrt(p_bar * (1 - p_bar) * (2 / n_arm))
theta <- (p_ctrl - p_trt) / se_full
cat(sprintf("Drift at full information: theta = %.3f (fixed-sample one-sided 5%% power %.1f%%)\n",
            theta, 100 * pnorm(theta - qnorm(0.95))))

Ks <- c(1, 2, 3, 5, 10, 20)

design_row <- function(K, sided, alpha = 0.05) {
  t_k <- (1:K) / K
  if (K == 1) {
    crit <- if (sided == 2) qnorm(1 - alpha / 2) else qnorm(1 - alpha)
  } else {
    crit <- rpact::getDesignGroupSequential(
      kMax = K, alpha = alpha, sided = sided, typeOfDesign = "asOF"
    )$criticalValues
  }
  corr <- outer(t_k, t_k, function(a, b) sqrt(pmin(a, b) / pmax(a, b)))
  mu <- theta * sqrt(t_k)

  no_cross <- function(shift) {
    if (K == 1) {
      if (sided == 2) pnorm(crit - shift[1]) - pnorm(-crit - shift[1])
      else pnorm(crit - shift[1])
    } else {
      lo <- if (sided == 2) -crit else rep(-Inf, K)
      pmvnorm(lower = lo, upper = crit, mean = shift, corr = corr,
              algorithm = GenzBretz(abseps = 1e-9, maxpts = 500000))[1]
    }
  }
  power <- 1 - no_cross(mu)
  size <- 1 - no_cross(rep(0, K))

  # expected stopping fraction E[min(tau, N)]/N under the alternative:
  # P(stop at look k) via successive no-crossing probabilities
  p_no_by_k <- vapply(seq_len(K), function(k) {
    if (k == 1) {
      if (sided == 2) pnorm(crit[1] - mu[1]) - pnorm(-crit[1] - mu[1])
      else pnorm(crit[1] - mu[1])
    } else {
      lo <- if (sided == 2) -crit[seq_len(k)] else rep(-Inf, k)
      pmvnorm(lower = lo, upper = crit[seq_len(k)], mean = mu[seq_len(k)],
              corr = corr[seq_len(k), seq_len(k), drop = FALSE],
              algorithm = GenzBretz(abseps = 1e-9, maxpts = 500000))[1]
    }
  }, numeric(1))
  p_stop_at <- c(1 - p_no_by_k[1], -diff(p_no_by_k))
  e_frac <- sum(p_stop_at * t_k) + p_no_by_k[K] * 1.0

  data.frame(
    K = K, sided = sided, alpha = alpha,
    final_boundary = crit[K],
    empirical_size = size,
    power = power,
    expected_n = e_frac * n_patients
  )
}

grid <- do.call(rbind, c(
  lapply(Ks, design_row, sided = 1),
  lapply(Ks, design_row, sided = 2)
))

write.csv(grid, file.path(table_dir, "power_vs_K.csv"), row.names = FALSE)

fmt <- function(x, d) sprintf(paste0("%.", d, "f"), x)
md <- data.frame(
  K = grid$K,
  Sided = ifelse(grid$sided == 2, "two-sided", "one-sided"),
  `Final boundary z` = fmt(grid$final_boundary, 3),
  `Size` = paste0(fmt(100 * grid$empirical_size, 2), "%"),
  `Power` = paste0(fmt(100 * grid$power, 1), "%"),
  `E[min(tau,N)] under alt` = format(round(grid$expected_n), big.mark = ","),
  check.names = FALSE
)
write_md <- function(x, path) {
  x[] <- lapply(x, as.character)
  writeLines(c(
    paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  ), path)
}
write_md(md, file.path(table_dir, "power_vs_K.md"))

cat("\nPower vs K (LD-OBF error spending, exactly calibrated):\n")
print(md, row.names = FALSE)
cat("\nWrote tables/power_vs_K.{csv,md}\n")
