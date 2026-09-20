# Freeze amendment 5 (author request 2026-08-10): e-RTe detection power under
# clean design-anticipated benefits, closing the ESM Section 2 gap (which was
# null-only). No censoring: lemma condition (b) holds; the event coin tilts
# because the effect is real.
#
# UNIFORMITY DIRECTIVE (author): power demonstrations are standalone and
# uniform across flavours. Scenarios and trial sizes MATCH the e-RTb tables
# (Supplementary Tables S1-S2): control event rate 40%; design ARR 5 pp
# (treatment 35%) and 10 pp (treatment 30%); N is the traditional
# two-proportion (chi-square) sample size at two-sided alpha 0.05,
# power.prop.test, for 80% and 90% target power:
#   5 pp: N = 2,942 (80%), 3,938 (90%)   [identical to Tables S1-S2]
#   10 pp: N = 712 (80%), 954 (90%)
# e-RTe sees only the events of those trials: per stream the event count is
# D ~ Binomial(N, pooled rate) and event arms are iid Bernoulli(q_true),
#   q_true = p_trt / (p_trt + p_ctrl)  (0.467 at 5 pp; 0.429 at 10 pp).
# Wealth is priced at the null coin 0.5 (kernel default).
#
# Bettors: adaptive e-RTe (kernel defaults, burn 30 / ramp 50) and the
# matched design event wager (the Table S2 "Design matched" analogue).
#
# NEW standalone script under FRESH seed 20260550; no existing random stream
# is touched, so all frozen tables reproduce unchanged.
#
# Output: tables/erte_power_check.{csv,md}

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("erte.R")))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

set.seed(20260550)
n_sims <- as.integer(Sys.getenv("ERTEPOW_NSIMS", "5000"))
threshold <- 20
p_ctrl <- 0.40

scenarios <- list()
si <- 1L
for (p_trt in c(0.35, 0.30)) {
  for (tp in c(0.80, 0.90)) {
    n_arm <- ceiling(power.prop.test(p1 = p_ctrl, p2 = p_trt, power = tp,
                                     sig.level = 0.05)$n)
    scenarios[[si]] <- list(
      arr = p_ctrl - p_trt, p_trt = p_trt, target_power = tp,
      N = 2L * as.integer(n_arm),
      pooled = (p_ctrl + p_trt) / 2,
      q_true = p_trt / (p_trt + p_ctrl)
    )
    si <- si + 1L
  }
}

first_cross <- function(w) { i <- which(w >= threshold); if (length(i)) i[1] else NA_integer_ }

res <- list(); ri <- 1L
t0 <- Sys.time()
for (sc in scenarios) {
  cross <- list(adaptive = logical(n_sims), design = logical(n_sims))
  fc <- list(adaptive = rep(NA_integer_, n_sims), design = rep(NA_integer_, n_sims))
  n_ev <- integer(n_sims)
  for (sim in seq_len(n_sims)) {
    d <- rbinom(1, sc$N, sc$pooled)
    arms <- rbinom(d, 1, sc$q_true)
    w_ad <- compute_eRTe(arms, wager = "adaptive")
    w_de <- compute_eRTe(arms, wager = "design",
                         p_trt_design = sc$p_trt, p_ctrl_design = p_ctrl)
    cross$adaptive[sim] <- any(w_ad$wealth >= threshold)
    cross$design[sim] <- any(w_de$wealth >= threshold)
    fc$adaptive[sim] <- first_cross(w_ad$wealth)
    fc$design[sim] <- first_cross(w_de$wealth)
    n_ev[sim] <- d
  }
  for (b in c("adaptive", "design")) {
    r <- mean(cross[[b]])
    res[[ri]] <- data.frame(
      design_arr_pp = 100 * sc$arr,
      target_power = sc$target_power,
      N = sc$N,
      median_events = median(n_ev),
      bettor = if (b == "adaptive") "Adaptive" else "Design matched",
      power = r,
      se = sqrt(r * (1 - r) / n_sims),
      median_first_crossing_event = median(fc[[b]], na.rm = TRUE),
      pct_of_events = median(fc[[b]], na.rm = TRUE) / median(n_ev),
      n_sims = n_sims
    )
    ri <- ri + 1L
  }
  cat(sprintf("ARR %.0fpp %du%% done\n", 100 * sc$arr, round(100 * sc$target_power)))
}
out <- do.call(rbind, res)
write.csv(out, file.path(table_dir, "erte_power_check.csv"), row.names = FALSE)

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
md <- data.frame(
  `Design ARR` = sprintf("%.1f pp", out$design_arr_pp),
  `Target power` = sprintf("%d%%", round(100 * out$target_power)),
  `N (chi-square)` = format(out$N, big.mark = ","),
  `Median events` = sprintf("%d", as.integer(out$median_events)),
  Bettor = out$bettor,
  `Power (W>=20)` = sprintf("%s (%.2fpp)", fmt_pct(out$power), 100 * out$se),
  `Median crossing (event)` = sprintf("%d", as.integer(out$median_first_crossing_event)),
  `% events` = fmt_pct(out$pct_of_events),
  check.names = FALSE
)
writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  sprintf("True benefit, no censoring (lemma condition (b) holds). Scenarios and trial sizes match Supplementary Tables S1-S2: control event rate 40%%; N is the two-proportion (chi-square) sample size at two-sided alpha 0.05 for the stated target power. Events per stream D~Binomial(N, pooled rate); event-arm coin q=p_trt/(p_trt+p_ctrl) (0.467 at 5 pp, 0.429 at 10 pp), priced at the null coin 0.5. %s simulated streams per scenario, seed 20260550. Median crossing among crossed streams.",
          format(n_sims, big.mark = ","))
), file.path(table_dir, "erte_power_check.md"))

cat(sprintf("Done in %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(md, row.names = FALSE)
