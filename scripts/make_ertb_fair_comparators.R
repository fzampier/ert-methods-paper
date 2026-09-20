# Fair-comparator and union-doctrine simulation for the BMC revision (T1-T5).
#
# Extends the submitted interruption-curve simulation (make_ertb_interruption_curve.R)
# with, on the SAME simulated trial streams (identical seeds; kernels consume no RNG):
#   T1  exactly calibrated error-spending boundaries (rpact canonical, mvtnorm
#       cross-check): Lan-DeMets OBF-type and Kim-DeMets rho-family, K=3 and K=5
#   T2  one-sided comparators matched to the directional design wager
#       (one-sided alpha=0.05 primary, 0.025 sensitivity; D1 decision)
#   T3  mixture two-sided e-RT: W_mix = 0.5*W_up + 0.5*W_down (convex combination
#       of test martingales is a test martingale; threshold 20 unchanged)
#   T4  harm-direction scenario (35% vs 40%) under a fresh seed (appended stream,
#       never inserted mid-loop) -> detectable-directions table
#   T5  union-rule doctrine (D2 decision): worst-case always-stop at W>=20;
#       Repair A "DSMB confirms" (one-sided z>=3.0 primary, z>=2.5 sensitivity);
#       Repair B "partitioned alpha" (trigger W>=100 with GS at alpha=0.04;
#       W>=200 with GS at alpha=0.045)
#
# The five submitted method rows (OBF/Pocock/HP approximate two-sided, adaptive,
# design) are recomputed on the identical streams and must reproduce the pinned
# numbers; this script writes NEW output files only and never touches the
# submitted artifacts.
#
# Outputs: tables/fair_comparators_first_crossings.csv
#          tables/fair_comparators_method_summary.csv   (with MC standard errors)
#          tables/fair_comparators_table2.md
#          tables/fair_comparators_boundaries.md
#          tables/union_doctrine.csv / union_doctrine.md
#          tables/detectable_directions.md

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(mvtnorm))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05
threshold <- 20            # 1/alpha, submitted e-RT threshold
n_sims <- as.integer(Sys.getenv("FAIRCOMP_NSIMS", "5000"))
n_patients <- 2690
p_ctrl <- 0.35
p_trt_alt <- 0.30
p_trt_null <- 0.35
p_trt_harm <- 0.40         # T4: reversal scenario, same magnitude, opposite direction

look3 <- round(n_patients * (1:3) / 3)
look5 <- round(n_patients * (1:5) / 5)
t3 <- (1:3) / 3
t5 <- (1:5) / 5

# ---------------------------------------------------------------------------
# Boundaries (T1/T2). rpact is canonical; mvtnorm recursion is the cross-check.
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("rpact", quietly = TRUE))

rpact_bounds <- function(kMax, alpha, sided, type, gammaA = NULL) {
  args <- list(kMax = kMax, alpha = alpha, sided = sided, typeOfDesign = type)
  if (!is.null(gammaA)) args$gammaA <- gammaA
  do.call(rpact::getDesignGroupSequential, args)$criticalValues
}

# mvtnorm spending-boundary solver (cross-check; one-sided convention)
solve_spending_bounds <- function(t, spend_cum, sided = 2) {
  K <- length(t)
  corr <- outer(t, t, function(a, b) sqrt(pmin(a, b) / pmax(a, b)))
  bounds <- numeric(K)
  no_cross <- function(b, k) {
    if (k == 1) {
      if (sided == 2) pnorm(b[1]) - pnorm(-b[1]) else pnorm(b[1])
    } else {
      lo <- if (sided == 2) -b else rep(-Inf, k)
      pmvnorm(lower = lo, upper = b, corr = corr[seq_len(k), seq_len(k), drop = FALSE],
              algorithm = GenzBretz(abseps = 1e-9, maxpts = 250000))[1]
    }
  }
  for (k in seq_len(K)) {
    f <- function(ck) (1 - no_cross(c(bounds[seq_len(k - 1)], ck), k)) - spend_cum[k]
    bounds[k] <- uniroot(f, c(0.3, 8.5), tol = 1e-7)$root
  }
  bounds
}
obf_spend_1s <- function(t, a) 2 * (1 - pnorm(qnorm(1 - a / 2) / sqrt(t)))

# canonical boundary sets
B <- list(
  ldobf_k3_2s_05  = rpact_bounds(3, 0.05,  2, "asOF"),
  ldobf_k5_2s_05  = rpact_bounds(5, 0.05,  2, "asOF"),
  ldobf_k3_1s_05  = rpact_bounds(3, 0.05,  1, "asOF"),
  ldobf_k5_1s_05  = rpact_bounds(5, 0.05,  1, "asOF"),
  ldobf_k3_1s_025 = rpact_bounds(3, 0.025, 1, "asOF"),
  ldobf_k5_1s_025 = rpact_bounds(5, 0.025, 1, "asOF"),
  kd1_k5_1s_05    = rpact_bounds(5, 0.05,  1, "asKD", gammaA = 1),
  kd3_k5_1s_05    = rpact_bounds(5, 0.05,  1, "asKD", gammaA = 3),
  # Repair B recalibrated primary monitors (LD-OBF K=3 two-sided)
  ldobf_k3_2s_04  = rpact_bounds(3, 0.04,  2, "asOF"),
  ldobf_k3_2s_045 = rpact_bounds(3, 0.045, 2, "asOF")
)

# cross-check the one-sided sets against the mvtnorm recursion
xchk <- list(
  ldobf_k3_1s_05  = solve_spending_bounds(t3, obf_spend_1s(t3, 0.05),  1),
  ldobf_k5_1s_05  = solve_spending_bounds(t5, obf_spend_1s(t5, 0.05),  1),
  ldobf_k3_1s_025 = solve_spending_bounds(t3, obf_spend_1s(t3, 0.025), 1)
)
for (nm in names(xchk)) {
  dev <- max(abs(B[[nm]] - xchk[[nm]]))
  cat(sprintf("boundary cross-check %s: rpact vs mvtnorm max dev %.4f\n", nm, dev))
  stopifnot(dev < 0.005)
}

# legacy approximate boundaries (submitted rows; unchanged)
obf_approx <- qnorm(1 - alpha / 2) / sqrt(t3)
hp_bounds <- c(rep(qnorm(1 - 0.001 / 2), 2), qnorm(1 - alpha / 2))
pocock_bounds <- rep(2.289, 3)

boundary_table <- bind_rows(lapply(names(B), function(nm) {
  data.frame(set = nm, look = seq_along(B[[nm]]), bound = B[[nm]])
}))

# ---------------------------------------------------------------------------
# Trial machinery (identical to the submitted script where shared)
# ---------------------------------------------------------------------------

simulate_trial <- function(n, p_trt, p_ctrl, p_random = 0.5) {
  treatment <- rbinom(n, 1, p_random)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, p_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, p_ctrl)
  list(treatment = treatment, outcome = outcome)
}

z_at_look <- function(treatment, outcome, n_look) {
  idx <- seq_len(n_look)
  trt <- treatment[idx] == 1
  ctrl <- treatment[idx] == 0
  n_t <- sum(trt)
  n_c <- sum(ctrl)
  if (n_t == 0 || n_c == 0) return(NA_real_)
  p_t <- mean(outcome[idx][trt])
  p_c <- mean(outcome[idx][ctrl])
  p_pool <- mean(outcome[idx])
  se <- sqrt(p_pool * (1 - p_pool) * (1 / n_t + 1 / n_c))
  if (!is.finite(se) || se <= 0) return(NA_real_)
  (p_c - p_t) / se   # positive Z = benefit (fewer events in treatment)
}

# generic GS first crossing on precomputed look z-values
gs_first_crossing <- function(z_by_look, look_n, bounds, sided) {
  for (k in seq_along(look_n)) {
    z <- z_by_look[[k]]
    if (!is.finite(z)) next
    hit <- if (sided == 2) abs(z) >= bounds[[k]] else z >= bounds[[k]]
    if (hit) {
      return(list(first_crossing = look_n[[k]],
                  direction = ifelse(z > 0, "benefit", "harm"),
                  statistic = z))
    }
  }
  list(first_crossing = NA_integer_, direction = NA_character_, statistic = NA_real_)
}

ert_first_crossing_at <- function(wealth, thr) {
  cross <- which(wealth >= thr)
  if (!length(cross)) NA_integer_ else cross[[1]]
}

arr_at <- function(treatment, outcome, n_look) {
  idx <- seq_len(n_look)
  trt <- treatment[idx] == 1
  ctrl <- treatment[idx] == 0
  if (!any(trt) || !any(ctrl)) return(NA_real_)
  mean(outcome[idx][ctrl]) - mean(outcome[idx][trt])
}

# method spec: name, kind (gs/ert), and parameters
GS_METHODS <- list(
  list(name = "OBF two-sided",                      looks = look3, bounds = obf_approx,          sided = 2),
  list(name = "Pocock two-sided",                   looks = look3, bounds = pocock_bounds,       sided = 2),
  list(name = "Haybittle-Peto",                     looks = look3, bounds = hp_bounds,           sided = 2),
  list(name = "LD-OBF K3 two-sided a=0.05",         looks = look3, bounds = B$ldobf_k3_2s_05,    sided = 2),
  list(name = "LD-OBF K5 two-sided a=0.05",         looks = look5, bounds = B$ldobf_k5_2s_05,    sided = 2),
  list(name = "LD-OBF K3 one-sided a=0.05",         looks = look3, bounds = B$ldobf_k3_1s_05,    sided = 1),
  list(name = "LD-OBF K5 one-sided a=0.05",         looks = look5, bounds = B$ldobf_k5_1s_05,    sided = 1),
  list(name = "LD-OBF K3 one-sided a=0.025",        looks = look3, bounds = B$ldobf_k3_1s_025,   sided = 1),
  list(name = "LD-OBF K5 one-sided a=0.025",        looks = look5, bounds = B$ldobf_k5_1s_025,   sided = 1),
  list(name = "KD rho=1 K5 one-sided a=0.05",       looks = look5, bounds = B$kd1_k5_1s_05,      sided = 1),
  list(name = "KD rho=3 K5 one-sided a=0.05",       looks = look5, bounds = B$kd3_k5_1s_05,      sided = 1)
)

run_scenario <- function(p_trt, scenario, seed) {
  set.seed(seed)
  method_rows <- vector("list", n_sims)
  union_rows <- vector("list", n_sims)

  for (sim in seq_len(n_sims)) {
    trial <- simulate_trial(n_patients, p_trt = p_trt, p_ctrl = p_ctrl)

    # z at every K3 and K5 look, computed once
    z3 <- lapply(look3, function(nl) z_at_look(trial$treatment, trial$outcome, nl))
    z5 <- lapply(look5, function(nl) z_at_look(trial$treatment, trial$outcome, nl))
    z_for <- function(looks) if (identical(looks, look3)) z3 else z5

    # wealth processes: adaptive, design-up (anticipated benefit), design-down (mirror)
    w_adaptive <- compute_eRT(trial$treatment, trial$outcome,
                              burn_in = 50, ramp = 100,
                              wager = "adaptive", kelly_fraction = 0.5)
    w_up <- compute_eRT(trial$treatment, trial$outcome,
                        wager = "design",
                        p_trt_design = p_trt_alt, p_ctrl_design = p_ctrl,
                        ramp_fixed = FALSE)
    w_down <- compute_eRT(trial$treatment, trial$outcome,
                          wager = "design",
                          p_trt_design = p_ctrl, p_ctrl_design = p_trt_alt,
                          ramp_fixed = FALSE)
    w_mix <- 0.5 * w_up + 0.5 * w_down

    res <- list()
    for (m in GS_METHODS) {
      res[[m$name]] <- gs_first_crossing(z_for(m$looks), m$looks, m$bounds, m$sided)
    }
    ert_res <- function(wealth, thr = threshold) {
      fc <- ert_first_crossing_at(wealth, thr)
      if (is.na(fc)) {
        list(first_crossing = NA_integer_, direction = NA_character_, statistic = NA_real_)
      } else {
        arr <- arr_at(trial$treatment, trial$outcome, fc)
        list(first_crossing = fc,
             direction = ifelse(is.finite(arr) && arr >= 0, "benefit", "harm"),
             statistic = arr)
      }
    }
    res[["e-RTb adaptive"]] <- ert_res(w_adaptive)
    res[["e-RTb design directional"]] <- ert_res(w_up)
    res[["e-RTb mixture two-sided"]] <- ert_res(w_mix)

    method_rows[[sim]] <- bind_rows(lapply(names(res), function(nm) {
      r <- res[[nm]]
      data.frame(sim = sim, scenario = scenario, method = nm,
                 first_crossing = r$first_crossing,
                 direction = r$direction,
                 statistic_at_crossing = r$statistic,
                 stringsAsFactors = FALSE)
    }))

    # ---- T5 union doctrine (driver = design e-RT; base monitor = LD-OBF K3 2s) ----
    gs_time <- res[["LD-OBF K3 two-sided a=0.05"]]$first_crossing
    fc20 <- ert_first_crossing_at(w_up, 20)
    fc100 <- ert_first_crossing_at(w_up, 100)
    fc200 <- ert_first_crossing_at(w_up, 200)
    z_alert <- if (!is.na(fc20)) z_at_look(trial$treatment, trial$outcome, fc20) else NA_real_

    # recalibrated monitors for Repair B
    gs04 <- gs_first_crossing(z3, look3, B$ldobf_k3_2s_04, 2)$first_crossing
    gs045 <- gs_first_crossing(z3, look3, B$ldobf_k3_2s_045, 2)$first_crossing

    stop_time <- function(a, b) suppressWarnings(min(c(a, b), na.rm = TRUE))
    fin <- function(x) if (is.finite(x)) x else NA_integer_

    union_rows[[sim]] <- data.frame(
      sim = sim, scenario = scenario,
      gs_time = gs_time, fc20 = fc20, fc100 = fc100, fc200 = fc200,
      z_at_alert = z_alert,
      worst_case = fin(stop_time(gs_time, fc20)),
      repair_a30 = fin(stop_time(gs_time, if (!is.na(fc20) && is.finite(z_alert) && z_alert >= 3.0) fc20 else NA_integer_)),
      repair_a25 = fin(stop_time(gs_time, if (!is.na(fc20) && is.finite(z_alert) && z_alert >= 2.5) fc20 else NA_integer_)),
      repair_b100 = fin(stop_time(gs04, fc100)),
      repair_b200 = fin(stop_time(gs045, fc200))
    )
  }

  list(methods = bind_rows(method_rows), union = bind_rows(union_rows))
}

cat(sprintf("Running %d sims x 3 scenarios (n=%d)...\n", n_sims, n_patients))
t0 <- Sys.time()
alt  <- run_scenario(p_trt_alt,  "Alternative: 35% vs 30%", seed = 20260542)
cat(sprintf("  alt done  (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
null <- run_scenario(p_trt_null, "Null: 35% vs 35%",        seed = 20260543)
cat(sprintf("  null done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
harm <- run_scenario(p_trt_harm, "Harm: 35% vs 40%",        seed = 20260544)
cat(sprintf("  harm done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

crossings <- bind_rows(alt$methods, null$methods, harm$methods)
unions <- bind_rows(alt$union, null$union, harm$union)

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

method_order <- c(vapply(GS_METHODS, `[[`, "", "name"),
                  "e-RTb adaptive", "e-RTb design directional", "e-RTb mixture two-sided")

method_summary <- crossings %>%
  mutate(method = factor(method, levels = method_order)) %>%
  group_by(scenario, method) %>%
  summarise(
    crossed_any = mean(!is.na(first_crossing)),
    crossed_benefit = mean(!is.na(first_crossing) & direction == "benefit"),
    crossed_harm = mean(!is.na(first_crossing) & direction == "harm"),
    se_any = sqrt(crossed_any * (1 - crossed_any) / n()),
    median_crossing = median(first_crossing, na.rm = TRUE),
    q25_crossing = quantile(first_crossing, 0.25, na.rm = TRUE, names = FALSE),
    q75_crossing = quantile(first_crossing, 0.75, na.rm = TRUE, names = FALSE),
    expected_stop = mean(ifelse(is.na(first_crossing), n_patients, pmin(first_crossing, n_patients))),
    n_sims = n(),
    .groups = "drop"
  )

fmt_pct <- function(x, d = 1) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f%%"), 100 * x))
fmt_se <- function(x) ifelse(is.na(x), "--", sprintf("(%.2fpp)", 100 * x))
fmt_num <- function(x, d = 0) ifelse(is.na(x), "--", format(round(x, d), big.mark = ",", trim = TRUE))

write_md_table <- function(x, path) {
  x[] <- lapply(x, as.character)
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), path)
  invisible(path)
}

# main Table 2 rebuild source
tab2 <- method_summary %>%
  select(scenario, method, crossed_any, se_any, median_crossing, expected_stop) %>%
  pivot_wider(names_from = scenario,
              values_from = c(crossed_any, se_any, median_crossing, expected_stop))
alt_nm <- "Alternative: 35% vs 30%"; null_nm <- "Null: 35% vs 35%"; harm_nm <- "Harm: 35% vs 40%"
tab2_out <- tab2 %>%
  transmute(
    Method = as.character(method),
    `Alt crossing` = paste(fmt_pct(.data[[paste0("crossed_any_", alt_nm)]]),
                           fmt_se(.data[[paste0("se_any_", alt_nm)]])),
    `Alt median n` = fmt_num(.data[[paste0("median_crossing_", alt_nm)]]),
    `Alt E[min(tau,N)]` = fmt_num(.data[[paste0("expected_stop_", alt_nm)]]),
    `Null crossing` = paste(fmt_pct(.data[[paste0("crossed_any_", null_nm)]]),
                            fmt_se(.data[[paste0("se_any_", null_nm)]])),
    `Harm crossing` = paste(fmt_pct(.data[[paste0("crossed_any_", harm_nm)]]),
                            fmt_se(.data[[paste0("se_any_", harm_nm)]])),
    `Harm median n` = fmt_num(.data[[paste0("median_crossing_", harm_nm)]])
  )

# detectable-directions table (T4/R4.m8): benefit- and harm-direction crossings
detect <- method_summary %>%
  filter(method %in% c("e-RTb adaptive", "e-RTb design directional", "e-RTb mixture two-sided",
                       "LD-OBF K3 two-sided a=0.05", "LD-OBF K3 one-sided a=0.05")) %>%
  select(scenario, method, crossed_any, crossed_benefit, crossed_harm) %>%
  mutate(across(c(crossed_any, crossed_benefit, crossed_harm), fmt_pct))

# union doctrine summary (T5)
union_summary <- unions %>%
  group_by(scenario) %>%
  summarise(
    `GS alone (LD-OBF K3 2s 5%)` = mean(!is.na(gs_time)),
    `Worst-case union (W>=20 always-stop)` = mean(!is.na(worst_case)),
    `Repair A confirm z>=3.0` = mean(!is.na(repair_a30)),
    `Repair A confirm z>=2.5` = mean(!is.na(repair_a25)),
    `Repair B W>=100 + GS a=0.04` = mean(!is.na(repair_b100)),
    `Repair B W>=200 + GS a=0.045` = mean(!is.na(repair_b200)),
    med_worst = median(worst_case, na.rm = TRUE),
    med_a30 = median(repair_a30, na.rm = TRUE),
    med_b100 = median(repair_b100, na.rm = TRUE),
    med_gs = median(gs_time, na.rm = TRUE),
    n_sims = n(),
    .groups = "drop"
  )

union_long <- union_summary %>%
  pivot_longer(cols = c(`GS alone (LD-OBF K3 2s 5%)`,
                        `Worst-case union (W>=20 always-stop)`,
                        `Repair A confirm z>=3.0`, `Repair A confirm z>=2.5`,
                        `Repair B W>=100 + GS a=0.04`, `Repair B W>=200 + GS a=0.045`),
               names_to = "Rule", values_to = "rate") %>%
  mutate(se = sqrt(rate * (1 - rate) / n_sims))

union_md <- union_long %>%
  select(scenario, Rule, rate, se) %>%
  mutate(`Stop rate` = paste(fmt_pct(rate), fmt_se(se))) %>%
  select(-rate, -se) %>%
  pivot_wider(names_from = scenario, values_from = `Stop rate`)

median_md <- union_summary %>%
  transmute(Scenario = scenario,
            `Median stop, GS alone` = fmt_num(med_gs),
            `Median stop, worst-case union` = fmt_num(med_worst),
            `Median stop, Repair A z3.0` = fmt_num(med_a30),
            `Median stop, Repair B W100` = fmt_num(med_b100))

# ---------------------------------------------------------------------------
# Write outputs (new files only)
# ---------------------------------------------------------------------------

write.csv(crossings, file.path(table_dir, "fair_comparators_first_crossings.csv"), row.names = FALSE)
write.csv(method_summary, file.path(table_dir, "fair_comparators_method_summary.csv"), row.names = FALSE)
write.csv(unions, file.path(table_dir, "union_doctrine_first_crossings.csv"), row.names = FALSE)
write.csv(union_summary, file.path(table_dir, "union_doctrine.csv"), row.names = FALSE)
write_md_table(tab2_out, file.path(table_dir, "fair_comparators_table2.md"))
write_md_table(detect, file.path(table_dir, "detectable_directions.md"))
write_md_table(union_md, file.path(table_dir, "union_doctrine.md"))
write_md_table(median_md, file.path(table_dir, "union_doctrine_medians.md"))
write_md_table(
  boundary_table %>% mutate(bound = sprintf("%.3f", bound)),
  file.path(table_dir, "fair_comparators_boundaries.md")
)

# pin-reproduction check against the submitted summary (same streams, same rows)
legacy <- read.csv(file.path(table_dir, "ertb_interruption_method_summary.csv"))
check_rows <- crossings %>%
  filter(method %in% c("OBF two-sided", "Pocock two-sided", "Haybittle-Peto",
                       "e-RTb adaptive", "e-RTb design directional"),
         scenario != harm_nm) %>%
  group_by(scenario, method) %>%
  summarise(crossed_any = mean(!is.na(first_crossing)), .groups = "drop")
merged <- merge(check_rows, legacy[, c("scenario", "method", "interrupted_any")],
                by = c("scenario", "method"))
max_dev <- max(abs(merged$crossed_any - merged$interrupted_any))
cat(sprintf("\nPin-reproduction check (5 legacy rows x 2 scenarios): max dev = %.5f\n", max_dev))
if (n_sims == 5000 && max_dev > 1e-12) {
  cat("WARNING: legacy rows did NOT reproduce exactly -- investigate before using outputs.\n")
} else if (n_sims == 5000) {
  cat("Legacy rows reproduce exactly. Additive extension is stream-consistent.\n")
}

cat(sprintf("\nDone in %.1f min. Outputs in %s\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), table_dir))
