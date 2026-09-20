# Freeze amendment 1 (R2.2, supports R4.9): batched updating and
# outcome-ascertainment delay.
#
# R2.2 asks how e-RT copes with two practical hurdles: (a) data are cleaned and
# analyzed in blocks, not after every patient; (b) there can be a long delay
# between an alert and the DSMB meeting / between enrollment and outcome
# ascertainment.
#
# (a) BATCHED INSPECTION. The wealth path is unchanged; the monitor simply
#     inspects W only at batch boundaries (every B patients). Inspecting a
#     subset of times can only miss crossings, so validity is preserved and the
#     procedure is conservative under the null; the cost is a detection delay
#     bounded by the batch size. Quantified here for B in {1, 50, 100, 250} on
#     the SAME simulated streams as the fair-comparator run (seeds 20260542/43;
#     the B=1 rows must reproduce the pinned crossing rates exactly, which
#     doubles as a stream-replication check).
#
# (b) ASCERTAINMENT DELAY. Under a constant lag of L patients (outcome of
#     patient i ascertained when patient i+L enrolls), the UPDATE ORDER is
#     unchanged, so the wealth path is identical patient-for-patient; every
#     alert simply occurs L enrollments later on the calendar clock. This needs
#     no new simulation: it is post-processing of the frozen first-crossing
#     distribution. Reported: share of alternative-scenario alerts that still
#     occur before enrollment ends (fc <= N - L) and the median calendar alert
#     time min(fc + L, end), for L in {0, 100, 300}.
#
# Outputs: tables/batched_updating_sensitivity.{csv,md}
#          tables/ascertainment_delay.{csv,md}

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))
suppressPackageStartupMessages(library(dplyr))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

n_sims <- as.integer(Sys.getenv("BATCH_NSIMS", "5000"))
n_patients <- 2690
p_ctrl <- 0.35
p_trt_alt <- 0.30
threshold <- 20
batches <- c(1, 50, 100, 250)

simulate_trial <- function(n, p_trt, p_ctrl, p_random = 0.5) {
  treatment <- rbinom(n, 1, p_random)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, p_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, p_ctrl)
  list(treatment = treatment, outcome = outcome)
}

batch_first_crossing <- function(wealth, B) {
  looks <- seq(B, n_patients, by = B)
  if (tail(looks, 1) != n_patients) looks <- c(looks, n_patients)
  w_at <- wealth[looks]
  hit <- which(w_at >= threshold)
  if (!length(hit)) NA_integer_ else looks[hit[1]]
}

run_scenario <- function(p_trt, scenario, seed) {
  set.seed(seed)
  out <- vector("list", n_sims)
  for (sim in seq_len(n_sims)) {
    trial <- simulate_trial(n_patients, p_trt = p_trt, p_ctrl = p_ctrl)
    w_adaptive <- compute_eRT(trial$treatment, trial$outcome,
                              burn_in = 50, ramp = 100,
                              wager = "adaptive", kelly_fraction = 0.5)
    w_up <- compute_eRT(trial$treatment, trial$outcome, wager = "design",
                        p_trt_design = p_trt_alt, p_ctrl_design = p_ctrl,
                        ramp_fixed = FALSE)
    w_down <- compute_eRT(trial$treatment, trial$outcome, wager = "design",
                          p_trt_design = p_ctrl, p_ctrl_design = p_trt_alt,
                          ramp_fixed = FALSE)
    w_mix <- 0.5 * w_up + 0.5 * w_down

    rows <- list()
    for (B in batches) {
      rows[[length(rows) + 1L]] <- data.frame(
        sim = sim, scenario = scenario, batch = B,
        method = c("e-RTb design directional", "e-RTb adaptive", "e-RTb mixture two-sided"),
        first_crossing = c(batch_first_crossing(w_up, B),
                           batch_first_crossing(w_adaptive, B),
                           batch_first_crossing(w_mix, B))
      )
    }
    out[[sim]] <- bind_rows(rows)
  }
  bind_rows(out)
}

cat(sprintf("Batched inspection: %d sims x 2 scenarios...\n", n_sims))
t0 <- Sys.time()
res <- bind_rows(
  run_scenario(p_trt_alt, "Alternative: 35% vs 30%", seed = 20260542),
  run_scenario(p_ctrl,    "Null: 35% vs 35%",        seed = 20260543)
)
cat(sprintf("  done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summ <- res %>%
  group_by(scenario, method, batch) %>%
  summarise(
    crossed = mean(!is.na(first_crossing)),
    se = sqrt(crossed * (1 - crossed) / n()),
    median_crossing = median(first_crossing, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summ, file.path(table_dir, "batched_updating_sensitivity.csv"), row.names = FALSE)

# stream-replication check: B=1 must equal the pinned fair-comparator rates
pin <- read.csv(file.path(table_dir, "fair_comparators_method_summary.csv")) %>%
  filter(method %in% unique(summ$method), scenario %in% unique(summ$scenario)) %>%
  select(scenario, method, pinned = crossed_any)
chk <- summ %>% filter(batch == 1) %>% inner_join(pin, by = c("scenario", "method"))
max_dev <- max(abs(chk$crossed - chk$pinned))
cat(sprintf("B=1 vs pinned fair-comparator rates: max dev = %.6f\n", max_dev))
if (n_sims == 5000 && max_dev > 1e-12) {
  cat("WARNING: B=1 rows did NOT reproduce the pinned rates -- investigate.\n")
} else if (n_sims == 5000) {
  cat("B=1 rows reproduce the pinned rates exactly (stream replication verified).\n")
}

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_n <- function(x) ifelse(is.na(x), "--", format(round(x), big.mark = ","))
md <- summ %>%
  arrange(scenario, method, batch) %>%
  transmute(
    Scenario = scenario, Method = method, `Inspect every` = batch,
    `Crossing rate` = sprintf("%s (%.2fpp)", fmt_pct(crossed), 100 * se),
    `Median crossing` = fmt_n(median_crossing)
  )
writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  "Same simulated streams as the fair-comparator run (seeds 20260542/43); wealth computed patient-by-patient, then inspected only at multiples of the batch size (final look always included) -- batched INSPECTION of a sequentially updated process, at the point-in-time wealth. Inspecting a subset of times preserves anytime validity and is conservative under the null. Detection delay is NOT bounded by the batch size under this rule: an excursion that crosses 20 and falls back between inspections is missed permanently, not delayed (design wager: 74.9% patient-level vs 71.1% at every-100 inspection). If per-update wealth is logged and the inspector instead checks the running maximum since the last inspection, every patient-level crossing is caught at the next boundary (delay then bounded by the batch size, crossing set identical, validity unchanged -- Ville bounds the supremum). [Caption corrected 2026-08-13, audit; simulated numbers untouched.]"
), file.path(table_dir, "batched_updating_sensitivity.md"))

# ---- (b) ascertainment delay: post-processing of the frozen crossings ----
fc <- read.csv(file.path(table_dir, "fair_comparators_first_crossings.csv")) %>%
  filter(method %in% c("e-RTb design directional", "e-RTb adaptive", "e-RTb mixture two-sided"),
         scenario == "Alternative: 35% vs 30%")
lags <- c(0, 100, 300)
delay <- bind_rows(lapply(lags, function(L) {
  fc %>% group_by(method) %>% summarise(
    lag = L,
    crossing_rate = mean(!is.na(first_crossing)),
    alert_before_enrollment_end = mean(!is.na(first_crossing) & first_crossing <= n_patients - L),
    median_calendar_alert = median(ifelse(is.na(first_crossing), NA, pmin(first_crossing + L, n_patients)), na.rm = TRUE),
    .groups = "drop")
}))
write.csv(delay, file.path(table_dir, "ascertainment_delay.csv"), row.names = FALSE)
md2 <- delay %>% transmute(
  Method = method, `Lag L (patients)` = lag,
  `Crossing rate (unchanged)` = fmt_pct(crossing_rate),
  `Alert before enrollment ends` = fmt_pct(alert_before_enrollment_end),
  `Median calendar alert` = fmt_n(median_calendar_alert)
)
writeLines(c(
  paste0("| ", paste(names(md2), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md2)), collapse = " | "), " |"),
  apply(md2, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  "Under a constant ascertainment lag of L patients the update order is unchanged, so the wealth path -- and therefore the crossing rate -- is identical; every alert occurs L enrollments later on the calendar clock (capped at enrollment end). Post-processing of the frozen alternative-scenario first crossings; no new simulation."
), file.path(table_dir, "ascertainment_delay.md"))

cat("Wrote batched_updating_sensitivity.{csv,md} and ascertainment_delay.{csv,md}\n")
