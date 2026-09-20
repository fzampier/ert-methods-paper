# Freeze amendment 3 (T10; R4.8 / R4.10): e-RTe under treatment-dependent
# event ascertainment, with NO treatment effect on the outcome itself.
#
# The Appendix lemma's condition (b) requires the observation process to be
# assignment-independent. This script demonstrates what happens when it is
# not: both arms have identical true event rates (the outcome null holds),
# but a fraction c of TREATMENT-arm events is never ascertained (censored
# before observation). Among observed events the arm-of-event coin is then
#   q = (1-c)/(2-c)  < 0.5,
# constant across events, so the event stream mimics a treatment benefit and
# e-RTe accumulates wealth: the alert fires on an observation-process
# difference, not a treatment effect -- exactly the "difference in the
# observation process" entry of the causes-of-crossing list (R4.3), and the
# reason condition (b) appears in the lemma rather than as a footnote.
#
# Bettors: adaptive e-RTe (kernel defaults, burn 30 / ramp 50) and the
# design event-only wager calibrated to a 35%-vs-30% anticipated benefit
# (event coin q* = 0.4615), which also fires because censoring mimics benefit.
#
# Event count: 430 observed events per simulated stream (ARMA-scale, cf. the
# 173/429 manuscript invariant). c in {0, 0.10, 0.25, 0.50}; c = 0 is the
# clean null and must sit at or below 5%.
#
# Output: tables/erte_censoring_null.{csv,md}

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("erte.R")))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

set.seed(20260548)
n_sims <- as.integer(Sys.getenv("CENSNULL_NSIMS", "5000"))
n_events <- 430L
threshold <- 20
cens <- c(0, 0.10, 0.25, 0.50)

q_obs <- function(c) (1 - c) / (2 - c)   # P(observed event is from treatment)

res <- list()
ri <- 1L
t0 <- Sys.time()
for (c_frac in cens) {
  q <- q_obs(c_frac)
  crossed_ad <- logical(n_sims)
  crossed_de <- logical(n_sims)
  for (sim in seq_len(n_sims)) {
    arms <- rbinom(n_events, 1, q)   # observed event-arm indicators X_j
    w_ad <- compute_eRTe(arms, wager = "adaptive")
    w_de <- compute_eRTe(arms, wager = "design",
                         p_trt_design = 0.30, p_ctrl_design = 0.35)
    crossed_ad[sim] <- any(w_ad$wealth >= threshold)
    crossed_de[sim] <- any(w_de$wealth >= threshold)
  }
  for (bet in c("adaptive", "design (35% vs 30% benefit)")) {
    cr <- if (bet == "adaptive") crossed_ad else crossed_de
    res[[ri]] <- data.frame(
      censoring_fraction = c_frac,
      q_observed = q,
      bettor = bet,
      crossing_rate = mean(cr),
      se = sqrt(mean(cr) * (1 - mean(cr)) / n_sims),
      n_sims = n_sims
    )
    ri <- ri + 1L
  }
}
out <- do.call(rbind, res)
write.csv(out, file.path(table_dir, "erte_censoring_null.csv"), row.names = FALSE)

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
md <- data.frame(
  `Treatment-event censoring c` = fmt_pct(out$censoring_fraction),
  `Observed event coin q` = sprintf("%.3f", out$q_observed),
  Bettor = out$bettor,
  `Crossing rate (W>=20)` = sprintf("%s (%.2fpp)", fmt_pct(out$crossing_rate), 100 * out$se),
  check.names = FALSE
)
writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  sprintf("Both arms share identical true event rates (outcome null); a fraction c of treatment-arm events is never ascertained, giving observed event coin q=(1-c)/(2-c). %s observed events per stream, %s simulated streams, seed 20260548. c=0 is the clean null. Crossings under c>0 are detections of the observation-process difference, not of a treatment effect (lemma condition (b); causes-of-crossing list).",
          format(n_events, big.mark = ","), format(n_sims, big.mark = ","))
), file.path(table_dir, "erte_censoring_null.md"))

cat(sprintf("Done in %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(md, row.names = FALSE)
