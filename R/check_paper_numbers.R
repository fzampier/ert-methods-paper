# Cross-check the numbers quoted in the article and its ESM against the
# generated CSV / JSON / markdown table outputs.
#
# Purpose: future seed/rerun drift in the generated tables will cause this
# script to fail loudly rather than silently desynchronising the article
# prose from the supporting tables. The hardcoded `expected` values mirror
# the article prose (gate rule: no number without a pin); the actual
# values come from the frozen artefacts under tables/ and, for the ARMA and
# FACTT pins, from the BioLINCC-derived files described in data/README.md.
# Those files are not part of this repository: without them the ARMA/FACTT
# pins that need them are reported as SKIPPED and do not fail the gate.
#
# NOTE: the `expected` values below are pinned to the current Monte-Carlo
# seeds. A deliberate reseed, a change in `n_sims`, or a change in the
# simulation design will move the table values and require re-pinning the
# expected list here. This is by design — a failure is a signal that prose
# and tables must be re-synced.
#
# Layout follows the article: construction and the e-RT family, practical
# deployment (Table 1, Table 2, worked alert, lead time), completed trials,
# discussion bridges, then ESM anchor pins and carried submission-freeze
# pins that back ESM Section 5's continuity tables, then the round-2
# premium decomposition (ESM Table S20) and the out-of-order p_i check.
#
# Run via `make check` or `Rscript R/check_paper_numbers.R`.

source("scripts/paths.R")

root <- paper_root()
tables_dir <- file.path(root, "tables")
results_dir <- file.path(validated_biolincc_dir(root), "results")
csv_dir <- file.path(validated_biolincc_dir(root), "csv_for_analysis")
have_trial_data <- all(file.exists(
  file.path(results_dir, c("arma_ert_summary_R.json", "factt_ert_summary_R.json")),
  file.path(csv_dir, "arma.csv")
))
if (have_trial_data) suppressPackageStartupMessages(library(jsonlite))

failures <- character(0)
skips <- character(0)

require_eq <- function(label, actual, expected, tol) {
  if (length(actual) != 1L || is.na(actual)) {
    failures[[length(failures) + 1L]] <<- sprintf("%s: actual is NA/absent", label)
    cat(sprintf("  [FAIL] %s: actual=NA\n", label))
    return(invisible(FALSE))
  }
  diff <- abs(actual - expected)
  if (diff <= tol) {
    cat(sprintf("  [PASS] %s: %.4f vs %.4f (tol %.4f)\n", label, actual, expected, tol))
    invisible(TRUE)
  } else {
    msg <- sprintf("%s: actual %.4f, expected %.4f, diff %.4f > tol %.4f",
                   label, actual, expected, diff, tol)
    failures[[length(failures) + 1L]] <<- msg
    cat(sprintf("  [FAIL] %s\n", msg))
    invisible(FALSE)
  }
}

require_true <- function(label, cond) {
  if (isTRUE(cond)) {
    cat(sprintf("  [PASS] %s\n", label))
  } else {
    failures[[length(failures) + 1L]] <<- label
    cat(sprintf("  [FAIL] %s\n", label))
  }
}

note_skip <- function(label, why) {
  skips[[length(skips) + 1L]] <<- sprintf("%s — %s", label, why)
  cat(sprintf("  [SKIP] %s — %s\n", label, why))
}

md_cells <- function(md_path, row_pattern) {
  lines <- readLines(md_path, warn = FALSE)
  hit <- grep(row_pattern, lines, value = TRUE)
  if (length(hit) == 0L) return(NULL)
  cells <- strsplit(hit[[1]], "\\|")[[1]]
  cells <- trimws(cells)
  cells[nchar(cells) > 0L]
}

md_num <- function(md_path, row_pattern, col_index, pct = FALSE) {
  cells <- md_cells(md_path, row_pattern)
  if (is.null(cells) || length(cells) < col_index) return(NA_real_)
  raw <- gsub(",", "", cells[[col_index]])
  m <- regmatches(raw, regexpr("-?[0-9]+\\.?[0-9]*", raw))
  if (length(m) == 0L) return(NA_real_)
  val <- as.numeric(m)
  if (pct) val / 100 else val
}

md_range <- function(md_path, row_pattern, col_index) {
  cells <- md_cells(md_path, row_pattern)
  if (is.null(cells) || length(cells) < col_index) return(c(NA_real_, NA_real_))
  raw <- gsub(",", "", cells[[col_index]])
  m <- regmatches(raw, gregexpr("[0-9]+\\.?[0-9]*", raw))[[1]]
  if (length(m) < 2L) return(c(NA_real_, NA_real_))
  as.numeric(m[1:2])
}

cat("=== Construction, the e-RT family: overall behaviour (main text; ESM Sections 1-3) ===\n")

oc <- read.csv(file.path(tables_dir, "ertb_section3_fresh_operating_characteristics.csv"),
               stringsAsFactors = FALSE)
oc_null <- oc[oc$scenario_type == "null", ]
oc_alt <- oc[oc$scenario_type == "alternative" & oc$policy == "matched", ]
# "0.5–4.6% for e-RTb across wager policies"
require_eq("Construction: e-RTb null min (prose 0.5%)", min(oc_null$rejection_rate), 0.0052, 0.0005)
require_eq("Construction: e-RTb null max (prose 4.6%)", max(oc_null$rejection_rate), 0.0464, 0.0005)
# "matched design wagers cross in 75–86% (e-RTb)"
require_eq("Construction: e-RTb matched power min (prose 75%)", min(oc_alt$rejection_rate), 0.7494, 0.003)
require_eq("Construction: e-RTb matched power max (prose 86%)", max(oc_alt$rejection_rate), 0.8610, 0.003)

cens <- read.csv(file.path(tables_dir, "erte_censoring_null.csv"), stringsAsFactors = FALSE)
# "2.4% for e-RTe" (uncensored null, adaptive)
require_eq("Construction: e-RTe null (prose 2.4%)",
           cens$crossing_rate[cens$censoring_fraction == 0 & cens$bettor == "adaptive"],
           0.0238, 0.0015)
# "null crossing rate rises smoothly to certainty" under one-arm non-ascertainment
require_true("Construction: e-RTe censoring detection: 50% censoring null rate >= 99.9%",
             all(cens$crossing_rate[cens$censoring_fraction == 0.5] >= 0.999))

epw <- read.csv(file.path(tables_dir, "erte_power_check.csv"), stringsAsFactors = FALSE)
epw_design <- epw[epw$bettor == "Design matched", ]
# "43–67% (e-RTe, which discards the non-event stream)"
require_eq("Construction: e-RTe design power min (prose 43%)", min(epw_design$power), 0.434, 0.003)
require_eq("Construction: e-RTe design power max (prose 67%)", max(epw_design$power), 0.670, 0.003)

ertc_t1 <- read.csv(file.path(tables_dir, "ertc_type1_simulation.csv"),
                    stringsAsFactors = FALSE)
t1_default <- ertc_t1[ertc_t1$setting == "Default (sign; burn 20, ramp 50)", ]
t1_conserv <- ertc_t1[grepl("^Conservative", ertc_t1$setting), ]
t1_caps <- ertc_t1[grepl("c_max=0\\.[48]", ertc_t1$setting), ]
# "4.1–4.5% (default) or 2.4–2.6% (conservative) for e-RTc" (cap-sensitivity
# rows are a separate ESM S10 block, pinned below)
require_eq("Construction: e-RTc null default min (prose 4.1%)", min(t1_default$type1), 0.041, 0.002)
require_eq("Construction: e-RTc null default max (prose 4.5%)", max(t1_default$type1), 0.045, 0.002)
require_eq("Construction: e-RTc null conservative min (prose 2.4%)", min(t1_conserv$type1), 0.024, 0.002)
require_eq("Construction: e-RTc null conservative max (prose 2.6%)", max(t1_conserv$type1), 0.026, 0.002)
# ESM S10 wager-cap sensitivity (R4.m9 freeze amendment): c_max 0.4/0.8 bracket
require_true("ESM S10 cap-sensitivity rows present (8)", nrow(t1_caps) == 8L)
require_eq("ESM S10 cap-sensitivity null min (3.7%)", min(t1_caps$type1), 0.037, 0.002)
require_eq("ESM S10 cap-sensitivity null max (4.8%)", max(t1_caps$type1), 0.048, 0.002)

cpw <- read.csv(file.path(tables_dir, "ertc_power_check.csv"), stringsAsFactors = FALSE)
cpw_design <- cpw[grepl("^Design matched", cpw$bettor), ]
# "71–85% (e-RTc)"
require_eq("Construction: e-RTc design power min (prose 71%)", min(cpw_design$power), 0.7128, 0.003)
require_eq("Construction: e-RTc design power max (prose 85%)", max(cpw_design$power), 0.8486, 0.003)

blocked_md <- file.path(tables_dir, "ertb_blocked_randomization_sensitivity.md")
# "permuted blocks with the conditional p_i (null 2.2%)"
require_eq("Construction: blocked conditional p_i null (prose 2.2%)",
           md_num(blocked_md, "correct conditional", 2, pct = TRUE), 0.022, 0.0015)
# "a bettor exploiting block bookkeeping against a naive 0.5 is a certain false
# signal" — ESM 1.3 worked numbers 100.0% / median 14; adaptive naive exposure 2.6%
require_eq("Construction/ESM 1.3 history-only naive bettor null (100.0%)",
           md_num(blocked_md, "History-only bettor", 2, pct = TRUE), 1.000, 0.0005)
require_eq("Construction/ESM 1.3 history-only naive bettor median crossing (14)",
           md_num(blocked_md, "History-only bettor", 3), 14, 0.5)
require_eq("Completed trials/ESM 1.3 adaptive naive-0.5 exposure (2.6%)",
           md_num(blocked_md, "practical exposure", 2, pct = TRUE), 0.026, 0.0015)

uneq <- read.csv(file.path(tables_dir, "ertb_unequal_allocation_check.csv"),
                 stringsAsFactors = FALSE)
# "under 2:1 allocation (design power 69.8%)"
require_eq("Construction: 2:1 allocation design power (prose 69.8%)",
           uneq$crossing_rate[grepl("^Alternative", uneq$scenario) &
                                grepl("design", uneq$method)],
           0.698, 0.003)

cat("\n=== Deployment: worked alert (Figure 3; ESM Section 4.1) ===\n")

walert_md <- file.path(tables_dir, "worked_alert_replay.md")
# "Wealth crosses 20 at patient 1,310 ... z is 2.56 ... 483 patients (18% of enrollment)"
require_eq("Deployment: worked alert crossing patient (1,310)",
           md_num(walert_md, "ALERT_N", 2), 1310, 0.5)
require_eq("Deployment: worked alert companion z (2.56)",
           md_num(walert_md, "ALERT_Z", 2), 2.56, 0.005)
require_eq("Deployment: worked alert lead (483 patients)",
           md_num(walert_md, "LEAD_N", 2), 483, 0.5)
require_eq("Deployment: worked alert lead (18% of enrollment)",
           md_num(walert_md, "LEAD_PCT", 2, pct = TRUE), 0.18, 0.005)

cat("\n=== Deployment, Table 1: fair comparators (fresh seeds) ===\n")

fair <- read.csv(file.path(tables_dir, "fair_comparators_method_summary.csv"),
                 stringsAsFactors = FALSE)
alt <- "Alternative: 35% vs 30%"
nul <- "Null: 35% vs 35%"
hrm <- "Harm: 35% vs 40%"
get_fair <- function(scen, meth) {
  r <- fair[fair$scenario == scen & fair$method == meth, ]
  if (nrow(r) == 0L) stop("Missing fair-comparators row: ", scen, " / ", meth)
  r
}
pin_t1 <- function(meth, alt_p, med, estop, null_p, harm_p, med_tol = 1) {
  ra <- get_fair(alt, meth); rn <- get_fair(nul, meth); rh <- get_fair(hrm, meth)
  require_eq(sprintf("Table 1 %s alt%%", meth), ra$crossed_any, alt_p, 0.003)
  require_eq(sprintf("Table 1 %s median crossing", meth), ra$median_crossing, med, med_tol)
  require_eq(sprintf("Table 1 %s E[min(tau,N)]", meth), ra$expected_stop, estop, 1)
  require_eq(sprintf("Table 1 %s null%%", meth), rn$crossed_any, null_p, 0.003)
  require_eq(sprintf("Table 1 %s harm-side%%", meth), rh$crossed_any, harm_p, 0.003)
}
pin_t1("LD-OBF K3 one-sided a=0.05",  0.865, 1793, 2156, 0.053, 0.000)
pin_t1("LD-OBF K5 one-sided a=0.05",  0.860, 1614, 2003, 0.052, 0.000)
pin_t1("KD rho=1 K5 one-sided a=0.05", 0.822, 1614, 1760, 0.052, 0.000)
pin_t1("LD-OBF K3 two-sided a=0.05",  0.783, 1793, 2316, 0.052, 0.766)
pin_t1("e-RTb design directional",    0.749, 1282, 1693, 0.041, 0.000)
pin_t1("e-RTb mixture two-sided",     0.644, 1506, 1952, 0.033, 0.614)
pin_t1("e-RTb adaptive",              0.459, 1343, 2081, 0.031, 0.444)
# Continuity rows (submission's approximate textbook boundaries, fresh streams)
pin_t1("OBF two-sided",               0.795, 1793, 2259, 0.059, 0.777)
pin_t1("Pocock two-sided",            0.716, 1793, 2012, 0.051, 0.694)
pin_t1("Haybittle-Peto",              0.789, 2690, 2505, 0.054, 0.772)

# Derived headline (deployment prose + Abstract): "beats design e-RTb on power by 11.6
# percentage points"
gap <- get_fair(alt, "LD-OBF K3 one-sided a=0.05")$crossed_any -
  get_fair(alt, "e-RTb design directional")$crossed_any
require_eq("Deployment/Abstract derived 11.6-pt power gap", gap, 0.116, 0.002)
# "The mixture ... pays about ten points relative to the design wager"
require_eq("Deployment: mixture-vs-design ~10-pt cost",
           get_fair(alt, "e-RTb design directional")$crossed_any -
             get_fair(alt, "e-RTb mixture two-sided")$crossed_any,
           0.105, 0.005)

cat("\n=== Deployment: time-resolved trade (cumulative at looks; ESM S15) ===\n")

cum <- read.csv(file.path(tables_dir, "fair_comparators_cumulative_at_looks.csv"),
                stringsAsFactors = FALSE)
get_cum <- function(scen, meth) {
  r <- cum[cum$scenario == scen & cum$method == meth, ]
  if (nrow(r) == 0L) stop("Missing cumulative row: ", scen, " / ", meth)
  r
}
# "by the first scheduled look design e-RTb has crossed in 20.5% ... against 5.2%
# for its calibrated one-sided rival, the curves are nearly level by two-thirds
# information (55.5 vs 54.3), and the scheduled designs have won by the final look"
require_eq("Deployment: design e-RTb cum at look 1 (20.5%)",
           get_cum(alt, "e-RTb design directional")$cum_897, 0.205, 0.003)
require_eq("Deployment: LD-OBF K3 1s cum at look 1 (5.2%)",
           get_cum(alt, "LD-OBF K3 one-sided a=0.05")$cum_897, 0.052, 0.003)
require_eq("ESM S15 KD rho=1 cum at look 1 (13.0%)",
           get_cum(alt, "KD rho=1 K5 one-sided a=0.05")$cum_897, 0.130, 0.003)
require_eq("Deployment: design e-RTb cum at look 2 (55.5%)",
           get_cum(alt, "e-RTb design directional")$cum_1793, 0.555, 0.003)
require_eq("Deployment: LD-OBF K3 1s cum at look 2 (54.3%)",
           get_cum(alt, "LD-OBF K3 one-sided a=0.05")$cum_1793, 0.543, 0.003)
require_eq("ESM S15 KD rho=1 cum at look 2 (54.7%)",
           get_cum(alt, "KD rho=1 K5 one-sided a=0.05")$cum_1793, 0.547, 0.003)

cat("\n=== Deployment, Figure 5: price of anytime validity (power_vs_K) ===\n")

pvk <- read.csv(file.path(tables_dir, "power_vs_K.csv"), stringsAsFactors = FALSE)
pvk1 <- pvk[pvk$sided == 1 & pvk$alpha == 0.05, ]
row_k <- function(k) pvk1[pvk1$K == k, ]
require_eq("Fig 5 K=1 power (86.9%)", row_k(1)$power, 0.8694, 0.001)
require_eq("Fig 5 K=1 boundary (1.645)", row_k(1)$final_boundary, 1.645, 0.0015)
require_eq("Fig 5 K=1 expected N (2,690)", row_k(1)$expected_n, 2690, 1)
require_eq("Fig 5 K=5 expected N (2,002)", row_k(5)$expected_n, 2002, 1)
require_eq("Fig 5 K=20 power (85.2%)", row_k(20)$power, 0.8524, 0.001)
require_eq("Fig 5 K=20 boundary (1.839)", row_k(20)$final_boundary, 1.839, 0.0015)
require_eq("Fig 5 K=20 expected N (1,848)", row_k(20)$expected_n, 1848, 1)
# Derived: "one look to twenty costs under two points of power"
require_true("Deployment: K1->K20 power cost < 2pp",
             (row_k(1)$power - row_k(20)$power) < 0.02)
# Derived: "cutting the expected stopping time by more than 800 patients"
require_true("Deployment: K1->K20 expected-N saving > 800",
             (row_k(1)$expected_n - row_k(20)$expected_n) > 800)
# Derived: "five looks already capture about four-fifths of that saving"
frac5 <- (row_k(1)$expected_n - row_k(5)$expected_n) /
  (row_k(1)$expected_n - row_k(20)$expected_n)
require_true("Deployment: K=5 captures ~4/5 of the saving (0.75–0.88)",
             frac5 > 0.75 && frac5 < 0.88)
# Derived: "the anytime premium is about 11 points" (e-RT horizontal below plateau)
prem <- mean(c(row_k(1)$power, row_k(20)$power)) -
  get_fair(alt, "e-RTb design directional")$crossed_any
require_true("Deployment/Discussion anytime premium ~11 points (0.095–0.125)",
             prem > 0.095 && prem < 0.125)

cat("\n=== Deployment, Table 2: union rule and repairs ===\n")

ud <- read.csv(file.path(tables_dir, "union_doctrine.csv"),
               stringsAsFactors = FALSE, check.names = FALSE)
ud_alt <- ud[grepl("^Alternative", ud$scenario), ]
ud_nul <- ud[grepl("^Null", ud$scenario), ]
pin_t2 <- function(col, null_p, alt_p) {
  require_eq(sprintf("Table 2 %s null stop rate", col), ud_nul[[col]], null_p, 0.003)
  require_eq(sprintf("Table 2 %s power", col), ud_alt[[col]], alt_p, 0.003)
}
pin_t2("GS alone (LD-OBF K3 2s 5%)",              0.052, 0.783)
pin_t2("Worst-case union (W>=20 always-stop)",    0.072, 0.812)
pin_t2("Repair A confirm z>=3.0",                 0.056, 0.784)
pin_t2("Repair A confirm z>=2.5",                 0.067, 0.800)
pin_t2("Repair B W>=100 + GS a=0.04",             0.043, 0.758)
pin_t2("Repair B W>=200 + GS a=0.045",            0.047, 0.770)
# "its false stops arrive early (median null stop at patient 1,953, versus 2,690)"
udm_md <- file.path(tables_dir, "union_doctrine_medians.md")
require_eq("Deployment: worst-case union median null stop (1,953)",
           md_num(udm_md, "^\\|\\s*Null", 3), 1953, 0.5)
require_eq("Deployment: schedule-alone median null stop (2,690)",
           md_num(udm_md, "^\\|\\s*Null", 2), 2690, 0.5)

cat("\n=== Deployment: union path decomposition (ESM S16) ===\n")

upd <- read.csv(file.path(tables_dir, "union_path_decomposition.csv"),
                stringsAsFactors = FALSE)
upd_alt <- upd[grepl("^Alternative", upd$scenario), ]
get_upd <- function(rule) {
  r <- upd_alt[upd_alt$rule == rule, ]
  if (nrow(r) == 0L) stop("Missing union-path row: ", rule)
  r
}
# "When the alert path does win — 7.1% of alternative trials — it wins by a
# median of 1,317 patients"
require_eq("Deployment: Repair A z3.0 alert-path share (7.1%)",
           get_upd("Repair A z>=3.0")$via_confirmed_alert, 0.071, 0.002)
require_eq("Deployment: Repair A z3.0 alert-path median lead (1,317)",
           get_upd("Repair A z>=3.0")$alert_won_median_lead, 1317, 0.5)
# "the confirmation gate blocks roughly nine of ten alerts at the alert instant"
blk <- get_upd("Repair A z>=3.0")$alert_unconfirmed / get_upd("Repair A z>=3.0")$alert_fired
require_true("Deployment: confirmation gate blocks ~9/10 alerts (0.85–0.95)",
             blk > 0.85 && blk < 0.95)
# "The gentler settings give the alert path a substantial share of stops" (ESM S16:
# 46.2% under z>=2.5, 41.2% under Repair B W>=100)
require_eq("ESM S16 Repair A z2.5 alert-path share (46.2%)",
           get_upd("Repair A z>=2.5")$via_confirmed_alert, 0.462, 0.003)
require_eq("ESM S16 Repair B W100 alert-path share (41.2%)",
           get_upd("Repair B W>=100 (+GS a=0.04)")$via_confirmed_alert, 0.412, 0.003)

cat("\n=== Deployment: alert lead time (Abstract; ESM S17) ===\n")

alt_md <- file.path(tables_dir, "alert_lead_time.md")
# "Among alternative trials where the scheduled design stopped, 91.8% had a
# preceding e-RT alert, at a median of 776 patients before the stop (IQR 445–1,105)"
require_eq("Deployment/Abstract alert precedes stop (91.8%)",
           md_num(alt_md, "^\\|\\s*Alternative", 3, pct = TRUE), 0.918, 0.003)
require_eq("Deployment/Abstract median lead (776 patients)",
           md_num(alt_md, "^\\|\\s*Alternative", 4), 776, 0.5)
iqr <- md_range(alt_md, "^\\|\\s*Alternative", 5)
require_eq("Deployment: lead IQR lower (445)", iqr[1], 445, 0.5)
require_eq("Deployment: lead IQR upper (1,105)", iqr[2], 1105, 0.5)
require_eq("ESM S17 median z at alert (2.54)",
           md_num(alt_md, "^\\|\\s*Alternative", 6), 2.54, 0.005)
# "Under the null, preceding alerts accompany only 39.1% of the rare scheduled
# stops (5.2%); under harm, none."
require_eq("Deployment: null GS-stop rate context (5.2%)",
           md_num(alt_md, "^\\|\\s*Null", 2, pct = TRUE), 0.052, 0.003)
require_eq("Deployment: null preceded share (39.1%)",
           md_num(alt_md, "^\\|\\s*Null", 3, pct = TRUE), 0.391, 0.003)
require_eq("Deployment: harm preceded share (0.0%)",
           md_num(alt_md, "^\\|\\s*Harm", 3, pct = TRUE), 0.000, 0.0005)

cat("\n=== Deployment: design-versus-actual grid (ESM S18) ===\n")

dap <- read.csv(file.path(tables_dir, "ertb_design_actual_power.csv"),
                stringsAsFactors = FALSE)
dap_rev <- dap[dap$direction == "reversed", ]
dap_des <- dap[dap$direction == "as-designed", ]
require_true("ESM S18 grid has 16 as-designed + 16 reversed cells",
             nrow(dap_rev) == 16L && nrow(dap_des) == 16L)
# "design e-RTb sits a few points below the fixed-sample power that set the
# sample size (75.8% against 80.0%)"
match_row <- dap_des[dap_des$p_ctrl == 0.4 & dap_des$design_arr == 0.05 &
                       dap_des$actual_arr == 0.05 & dap_des$target_power == 0.8, ]
require_eq("Deployment: as-designed design e-RTb power (75.8%)",
           match_row$ertb_design_power, 0.758, 0.003)
require_eq("Deployment: as-designed frequentist power (80.0%)",
           match_row$actual_frequentist_power, 0.800, 0.003)
# "design e-RTb crossed in 0.0–0.5% of trials in every one of the 16 reversal cells"
require_true("Deployment: reversal: design e-RTb <= 0.5% in all 16 cells",
             max(dap_rev$ertb_design_power) <= 0.005 + 1e-9)
require_true("Deployment: reversal: design e-RTb >= 0.0% (sanity)",
             min(dap_rev$ertb_design_power) >= 0)
# Artifact truth for the adaptive column in the reversal cells (ESM S18):
require_eq("ESM S18 reversal adaptive min (6.7%)", min(dap_rev$ertb_adaptive_power), 0.067, 0.002)
require_eq("ESM S18 reversal adaptive max (100%)", max(dap_rev$ertb_adaptive_power), 1.000, 0.002)
# Main-text claim (corrected 2026-08-13, was "48–100%"): "the adaptive wager
# retained 91–112% of its mirrored-benefit power in those same cells".
# Definition: per reversal cell, adaptive reversal power / adaptive power at
# the mirrored as-designed cell (matched on p_ctrl, design_arr, |actual_arr|,
# target_power, n_patients).
dap_asd <- dap[dap$direction == "as-designed", ]
mir_key <- function(p, da, aa, tp, n) paste(p, da, aa, tp, n)
mir_idx <- match(mir_key(dap_rev$p_ctrl, dap_rev$design_arr, abs(dap_rev$actual_arr),
                         dap_rev$target_power, dap_rev$n_patients),
                 mir_key(dap_asd$p_ctrl, dap_asd$design_arr, dap_asd$actual_arr,
                         dap_asd$target_power, dap_asd$n_patients))
require_true("Deployment: mirrored-benefit retention: all 16 reversal cells matched",
             !anyNA(mir_idx) && length(mir_idx) == 16)
mir_ratio <- dap_rev$ertb_adaptive_power / dap_asd$ertb_adaptive_power[mir_idx]
require_eq("Deployment: adaptive mirrored-benefit retention min (91%)", min(mir_ratio), 0.913, 0.005)
require_eq("Deployment: adaptive mirrored-benefit retention max (112%)", max(mir_ratio), 1.121, 0.005)

cat("\n=== Deployment: Bayesian sensitivity (Figure 6) ===\n")

bayes_method <- read.csv(file.path(tables_dir, "bayesian_monitor_method_summary.csv"),
                         stringsAsFactors = FALSE)
get_bayes <- function(scen, thr) {
  r <- bayes_method[bayes_method$scenario == scen & abs(bayes_method$threshold - thr) < 1e-9, ]
  if (nrow(r) == 0L) stop("Missing bayes row: scenario=", scen, " threshold=", thr)
  r
}
# "stopping at posterior probability 0.98 signals in 78.9% (alternative) and 4.7%
# (null) by itself"
require_eq("Deployment: Bayes 0.98 standalone alt (78.9%)",
           get_bayes(alt, 0.98)$interrupted_any, 0.789, 0.005)
require_eq("Deployment: Bayes 0.98 standalone null (4.7%)",
           get_bayes(nul, 0.98)$interrupted_any, 0.047, 0.003)

union <- read.csv(file.path(tables_dir, "bayesian_ert_union_summary.csv"),
                  stringsAsFactors = FALSE)
union$threshold <- as.numeric(union$threshold)
union$combined_signal <- as.numeric(union$combined_signal)
get_union <- function(scen, thr, ert) {
  r <- union[union$scenario == scen & abs(union$threshold - thr) < 1e-9 &
               union$ert_method == ert, ]
  if (nrow(r) == 0L) stop("Missing union row: scenario=", scen,
                          " threshold=", thr, " ert=", ert)
  r
}
# "as an automatic union with design e-RTb, raising the threshold from 0.98 to
# 0.99 returns the null rate from 5.5% to 4.2%" (Figure 6 also plots 0.985)
for (spec in list(
  list(thr = 0.98,  alt_exp = 0.800, nul_exp = 0.055),
  list(thr = 0.985, alt_exp = 0.777, nul_exp = 0.048),
  list(thr = 0.99,  alt_exp = 0.752, nul_exp = 0.042)
)) {
  thr <- spec$thr
  require_eq(sprintf("Fig 6 Bayes %.3f + design e-RTb alt", thr),
             get_union(alt, thr, "e-RTb design directional")$combined_signal,
             spec$alt_exp, 0.008)
  require_eq(sprintf("Fig 6 Bayes %.3f + design e-RTb null", thr),
             get_union(nul, thr, "e-RTb design directional")$combined_signal,
             spec$nul_exp, 0.005)
}

cat("\n=== Completed trials: ARMA (validated BioLINCC) ===\n")

# "55.6% of random-order permutations crossed for e-RTb" (+ ESM Section 6 BWA rows)
arma_md <- file.path(tables_dir, "arma_adaptive_permutation_sensitivity.md")
require_eq("Completed trials: ARMA mortality e-RTb perm crossings (55.6%)",
           md_num(arma_md, "^\\|\\s*Mortality e-RTb\\b", 3, pct = TRUE), 0.556, 0.005)
require_eq("ESM Section 6 BWA-28 e-RTb perm crossings (94.2%)",
           md_num(arma_md, "^\\|\\s*BWA-28 failure e-RTb\\b", 3, pct = TRUE), 0.942, 0.005)
require_eq("ESM Section 6 BWA-28 e-RTb observed crossing index (751)",
           md_num(arma_md, "^\\|\\s*BWA-28 failure e-RTb\\b", 2), 751, 0.5)

if (have_trial_data) {
  arma <- fromJSON(file.path(results_dir, "arma_ert_summary_R.json"))
  arma_csv <- read.csv(file.path(csv_dir, "arma.csv"), stringsAsFactors = FALSE)
  # "the analysis cohort has 861 patients (429 control, 432 low tidal volume)"
  require_eq("Completed trials: ARMA cohort n (861)", arma$n_total, 861, 0.5)
  require_eq("Completed trials: ARMA control arm (429)", sum(arma_csv$arm == 0), 429, 0.5)
  require_eq("Completed trials: ARMA low-tidal-volume arm (432)", sum(arma_csv$arm == 1), 432, 0.5)
  # Counts quoted in the article: 173 control vs 134 treatment deaths
  require_eq("Completed trials: ARMA control deaths (173)", arma$ertd$n_deaths_control, 173, 0.5)
  require_eq("Completed trials: ARMA treatment deaths (134)", arma$ertd$n_deaths_treatment, 134, 0.5)
  # "mortality crossed with neither adaptive e-RTb nor e-RTe"
  require_true("Completed trials: ARMA mortality e-RTb did not cross (masked order)",
               !isTRUE(arma$ert_binary$crossed_threshold))
  require_true("Completed trials: ARMA mortality e-RTe did not cross (masked order)",
               !isTRUE(arma$ertd$crossed_threshold))
  # "the event-only stream ... a peak W ... of about 16.8"
  require_eq("Completed trials: ARMA e-RTe mortality peak W (prose ~16.8)",
             arma$ertd$peak_wealth, 16.77, 0.05)
  # "BWA-28 crossed with both variants"
  require_true("Completed trials: ARMA BWA-28 crossed (patient-level e-RTb)",
               isTRUE(arma$bwa28_binary$crossed_threshold))
  require_true("Completed trials: ARMA BWA-28 crossed (event-only e-RTe)",
               isTRUE(arma$bwa28_ertd$crossed_threshold))
  # ESM Section 6 carried: ARR on mortality (9.3 pp; JSON rd is treatment-minus-control)
  require_eq("ESM Section 6 ARMA mortality ARR (9.3 pp)", -100 * arma$ert_binary$rd, 9.3, 0.1)
} else {
  note_skip("Completed trials: ARMA cohort, deaths, observed-order crossings, peak W, ARR (11 pins)",
            "needs the BioLINCC-derived files (data/README.md)")
}

cat("\n=== Completed trials: FACTT (validated BioLINCC) ===\n")

# "the fluid VFD endpoint crossed in the observed masked order and in 98.9% of
# permutations"
factt_md <- file.path(tables_dir, "factt_adaptive_permutation_sensitivity.md")
require_eq("Completed trials: FACTT fluid VFD-28 perm crossings (98.9%)",
           md_num(factt_md, "^\\|\\s*Fluid VFD-28 e-RTc\\b", 3, pct = TRUE), 0.989, 0.003)
# "a same-tuning null crossing rate of about 2.5%" (conservative e-RTc, zero-
# inflated integer DGP; descriptive contrast)
require_eq("Completed trials: same-tuning null rate (~2.5%; ZII conservative)",
           ertc_t1$type1[ertc_t1$dgp == "Zero-inflated integer" &
                           grepl("^Conservative", ertc_t1$setting)],
           0.025, 0.002)
# "The catheter VFD trajectory ... 2.7% of permutations crossed"
require_eq("Completed trials: FACTT catheter perm crossings (2.7%)",
           md_num(factt_md, "^\\|\\s*Catheter VFD-28 e-RTc\\b", 3, pct = TRUE), 0.027, 0.002)

if (have_trial_data) {
  factt <- fromJSON(file.path(results_dir, "factt_ert_summary_R.json"))
  # "the validated cohort is the full randomized sample of 1,000 patients (497
  # liberal versus 503 conservative; 487 central venous versus 513 pulmonary artery)"
  require_eq("Completed trials: FACTT n_total (1,000)", factt$n_total, 1000, 0.5)
  require_eq("Completed trials: FACTT n_liberal (497)", factt$n_liberal, 497, 0.5)
  require_eq("Completed trials: FACTT n_conservative (503)", factt$n_conservative, 503, 0.5)
  require_eq("Completed trials: FACTT n_cvc (487)", factt$n_cvc, 487, 0.5)
  require_eq("Completed trials: FACTT n_pac (513)", factt$n_pac, 513, 0.5)
  require_true("Completed trials: FACTT fluid VFD-28 crossed (masked order)",
               isTRUE(factt$ertc_vfd28$fluid_comparison$crossed_threshold))
  require_eq("ESM Section 6 fluid VFD-28 observed crossing index (131)",
             factt$ertc_vfd28$fluid_comparison$crossing_patient, 131, 0.5)
  # "The catheter VFD trajectory ... peak W ≈ 2.42 ... drifted to W ≈ 0.05 by trial end"
  require_true("Completed trials: FACTT catheter VFD-28 did not cross",
               !isTRUE(factt$ertc_vfd28$catheter_comparison$crossed_threshold))
  require_eq("Completed trials: FACTT catheter peak W (2.42)",
             factt$ertc_vfd28$catheter_comparison$peak_evalue, 2.42, 0.005)
  # "W ≈ 0.05" in the prose matches the frozen 0.0543.
  require_eq("Completed trials: FACTT catheter final W (0.054; prose '~0.05')",
             factt$ertc_vfd28$catheter_comparison$final_evalue, 0.054, 0.005)
} else {
  note_skip("Completed trials: FACTT cohort, observed-order crossings, catheter peak and final W (10 pins)",
            "needs the BioLINCC-derived files (data/README.md)")
}

cat("\n=== Discussion bridges: Type M, batching, delay ===\n")

# Type M numbers are carried from the submission freeze (per-patient policies vs
# comparators; ESM S19 continuity family)
method_old <- read.csv(file.path(tables_dir, "ertb_interruption_method_summary.csv"),
                       stringsAsFactors = FALSE)
get_old <- function(scen, meth) {
  r <- method_old[method_old$scenario == scen & method_old$method == meth, ]
  if (nrow(r) == 0L) stop("Missing submission-freeze row: ", scen, " / ", meth)
  r
}
require_eq("Discussion: adaptive e-RTb Type M (1.26)",
           get_old(alt, "e-RTb adaptive")$median_type_m, 1.26, 0.02)
require_eq("Discussion: design e-RTb Type M (1.19)",
           get_old(alt, "e-RTb design directional")$median_type_m, 1.19, 0.02)
require_eq("Discussion: Pocock Type M (1.07)",
           get_old(alt, "Pocock two-sided")$median_type_m, 1.07, 0.02)
require_eq("Discussion: OBF Type M (~1.00, fires at final look)",
           get_old(alt, "OBF two-sided")$median_type_m, 1.00, 0.02)
require_eq("Discussion: Haybittle-Peto Type M (~1.00, fires at final look)",
           get_old(alt, "Haybittle-Peto")$median_type_m, 1.00, 0.02)

bat <- read.csv(file.path(tables_dir, "batched_updating_sensitivity.csv"),
                stringsAsFactors = FALSE)
bat_des <- bat[grepl("^Alternative", bat$scenario) &
                 bat$method == "e-RTb design directional", ]
bat_des_null <- bat[grepl("^Null", bat$scenario) &
                      bat$method == "e-RTb design directional", ]
# "every-100-patient inspection lowered the design wager crossing rate from 74.9%
# to 71.1% (median alert 1,282 to 1,400 ...) while null rates fell" (ESM S11)
require_eq("Discussion: batching design per-patient rate (74.9%)",
           bat_des$crossed[bat_des$batch == 1], 0.749, 0.003)
require_eq("Discussion: batching design batch-100 rate (71.1%)",
           bat_des$crossed[bat_des$batch == 100], 0.711, 0.003)
require_eq("Discussion: batching design per-patient median (1,282)",
           bat_des$median_crossing[bat_des$batch == 1], 1282, 1)
require_eq("Discussion: batching design batch-100 median (1,400)",
           bat_des$median_crossing[bat_des$batch == 100], 1400, 1)
require_true("Discussion: batching: null rate falls under batch-100 inspection",
             bat_des_null$crossed[bat_des_null$batch == 100] <
               bat_des_null$crossed[bat_des_null$batch == 1])
# "even with a 300-patient lag, a design alert preceded the end of enrollment
# in 69.6% of all simulated trials — 93.0% of the trials that alerted at all"
# (all simulated trials as denominator; ESM S12)
lag <- read.csv(file.path(tables_dir, "ascertainment_delay.csv"),
                stringsAsFactors = FALSE)
lag_des <- lag[lag$method == "e-RTb design directional" & lag$lag == 300, ]
require_eq("Discussion: 300-patient lag: alert before enrollment end, all trials (69.6%)",
           lag_des$alert_before_enrollment_end, 0.696, 0.003)
require_eq("Discussion: 300-patient lag: among alerting trials (derived 93.0%)",
           lag_des$alert_before_enrollment_end / lag_des$crossing_rate,
           0.930, 0.003)

cat("\n=== ESM S13: calibrated boundary values ===\n")

bnd_md <- file.path(tables_dir, "fair_comparators_boundaries.md")
bnd <- readLines(bnd_md, warn = FALSE)
bnd_val <- function(set, look) {
  hit <- grep(sprintf("^\\|\\s*%s\\s*\\|\\s*%d\\s*\\|", set, look), bnd, value = TRUE)
  if (length(hit) == 0L) return(NA_real_)
  as.numeric(trimws(strsplit(hit[[1]], "\\|")[[1]][4]))
}
# One-sided K=3 at alpha=0.05: 3.200 / 2.141 / 1.695 (quoted in ESM S13 text)
require_eq("ESM S13 1s K3 look 1 (3.200)", bnd_val("ldobf_k3_1s_05", 1), 3.200, 0.0015)
require_eq("ESM S13 1s K3 look 2 (2.141)", bnd_val("ldobf_k3_1s_05", 2), 2.141, 0.0015)
require_eq("ESM S13 1s K3 look 3 (1.695)", bnd_val("ldobf_k3_1s_05", 3), 1.695, 0.0015)
# Two-sided K=3 at alpha=0.05: 3.710 / 2.511 / 1.993
require_eq("ESM S13 2s K3 look 1 (3.710)", bnd_val("ldobf_k3_2s_05", 1), 3.710, 0.0015)
require_eq("ESM S13 2s K3 look 2 (2.511)", bnd_val("ldobf_k3_2s_05", 2), 2.511, 0.0015)
require_eq("ESM S13 2s K3 look 3 (1.993)", bnd_val("ldobf_k3_2s_05", 3), 1.993, 0.0015)
# Repair-B recalibrated schedules exist (anchor one value each)
require_eq("ESM S13 2s K3 a=0.04 look 1 (3.863)", bnd_val("ldobf_k3_2s_04", 1), 3.863, 0.0015)
require_eq("ESM S13 2s K3 a=0.045 look 1 (3.783)", bnd_val("ldobf_k3_2s_045", 1), 3.783, 0.0015)

cat("\n=== ESM S19: submission-freeze continuity (conditional comparisons) ===\n")

joint <- read.csv(file.path(tables_dir, "ertb_interruption_joint_summary.csv"),
                  stringsAsFactors = FALSE)
get_joint <- function(scen, comp) {
  r <- joint[joint$scenario == scen & joint$comparison == comp, ]
  if (nrow(r) == 0L) stop("Missing joint row: scenario=", scen, " comparison=", comp)
  r
}
# Carried shared-crossing conditionals (submission Table 2; now ESM S19)
require_eq("ESM S19 design earlier | both (OBF) (99.5%)",
           get_joint(alt, "OBF + design e-RTb")$ert_before_base_given_both, 0.995, 0.005)
require_eq("ESM S19 design earlier | both (Pocock) (90.5%)",
           get_joint(alt, "Pocock + design e-RTb")$ert_before_base_given_both, 0.905, 0.01)
require_eq("ESM S19 design earlier | both (HP) (99.9%)",
           get_joint(alt, "HP + design e-RTb")$ert_before_base_given_both, 0.999, 0.005)
require_eq("ESM S19 adaptive earlier | both (OBF) (80.6%)",
           get_joint(alt, "OBF + adaptive e-RTb")$ert_before_base_given_both, 0.806, 0.01)
require_eq("ESM S19 adaptive earlier | both (Pocock) (62.3%)",
           get_joint(alt, "Pocock + adaptive e-RTb")$ert_before_base_given_both, 0.623, 0.01)
require_eq("ESM S19 adaptive earlier | both (HP) (98.0%)",
           get_joint(alt, "HP + adaptive e-RTb")$ert_before_base_given_both, 0.980, 0.01)
# Carried either-trigger null unions (Section 3.2 of the original submission; now ESM S19 context)
require_eq("ESM S19 OBF+design null union (7.7%)",
           get_joint(nul, "OBF + design e-RTb")$combined, 0.077, 0.003)
require_eq("ESM S19 Pocock+design null union (6.7%)",
           get_joint(nul, "Pocock + design e-RTb")$combined, 0.067, 0.003)
require_eq("ESM S19 HP+design null union (7.4%)",
           get_joint(nul, "HP + design e-RTb")$combined, 0.074, 0.003)

# Carried Pocock matched-power robustness (submission; now ESM S19)
pmp <- read.csv(file.path(tables_dir, "pocock_matched_power_summary.csv"),
                stringsAsFactors = FALSE)
get_pmp <- function(metric) {
  r <- pmp[pmp$metric == metric, ]
  if (nrow(r) == 0L) stop("Missing matched-power row: ", metric)
  r$value[[1]]
}
require_eq("ESM S19 Pocock alt% @ N=3084 (78.3%)",
           get_pmp("pocock_alt_interruption"), 0.783, 0.008)
require_eq("ESM S19 Pocock null% @ N=3084 (4.9%)",
           get_pmp("pocock_null_interruption"), 0.049, 0.005)
require_eq("ESM S19 design e-RTb alt% @ N=3084 (80.0%)",
           get_pmp("ertb_design_alt_interruption"), 0.800, 0.008)
require_eq("ESM S19 design earlier than Pocock @ N=3084 (92.8%)",
           get_pmp("ertb_design_earlier_than_pocock"), 0.928, 0.01)

cat("\n=== ESM S20: the design-wager premium decomposed (round 2) ===\n")

prem_tab <- read.csv(file.path(tables_dir, "premium_decomposition.csv"), stringsAsFactors = FALSE)
prem_meta <- readLines(file.path(tables_dir, "premium_decomposition_meta.csv"))
get_prem <- function(pattern, col = "estimate") {
  r <- prem_tab[[col]][grepl(pattern, prem_tab$quantity, fixed = TRUE)]
  if (length(r) != 1L) stop("Missing or duplicated premium row: ", pattern)
  r
}
get_meta <- function(key) {
  as.numeric(sub(paste0(key, ","), "", prem_meta[grepl(paste0("^", key, ","), prem_meta)]))
}
half <- 0.00051
require_eq("ESM S20 fixed-N z-test power (87.2%)", get_prem("fixed-N one-sided z-test"), 0.872, half)
require_eq("ESM S20 fixed-N z-test MC SE (0.5)", get_prem("fixed-N one-sided z-test", "mc_se"), 0.005, half)
require_eq("ESM S20 fixed-N final-wealth test power (87.6%)", get_prem("fixed-N final-wealth test"), 0.876, half)
require_eq("ESM S20 final-wealth MC SE (0.5)", get_prem("fixed-N final-wealth test", "mc_se"), 0.005, half)
require_eq("ESM S20 anytime crossing (74.9%)", get_prem("anytime crossing sup W"), 0.749, half)
require_eq("ESM S20 anytime MC SE (0.6)", get_prem("anytime crossing sup W", "mc_se"), 0.006, half)
require_eq("ESM S20 terminal W_N >= 20 (61.8%)", get_prem("final wealth W_N >= 20"), 0.618, half)
require_eq("ESM S20 terminal MC SE (0.7)", get_prem("final wealth W_N >= 20", "mc_se"), 0.007, half)
require_eq("ESM S20 null z-test (5.3%)", get_prem("fixed-N one-sided z-test", "null_rate"), 0.053, half)
require_eq("ESM S20 null final-wealth test (5.3%)", get_prem("fixed-N final-wealth test", "null_rate"), 0.053, half)
require_eq("ESM S20 null anytime crossing (4.1%)", get_prem("anytime crossing sup W", "null_rate"), 0.041, half)
require_eq("ESM S20 null terminal W_N >= 20 (0.8%)", get_prem("final wealth W_N >= 20", "null_rate"), 0.008, half)
require_eq("ESM S20 representation gap (-0.4 pt)", get_prem("representation:"), -0.004, half)
require_eq("ESM S20 representation SE (0.1)", get_prem("representation:", "mc_se"), 0.001, half)
require_eq("ESM S20 boundary gap (12.7 pt)", get_prem("boundary/horizon:"), 0.127, half)
require_eq("ESM S20 boundary SE (0.5)", get_prem("boundary/horizon:", "mc_se"), 0.005, half)
# "a test calibrated for a fixed N rejects at a final wealth of about 1.9"
require_eq("Discussion: fixed-N critical final wealth (1.93)", get_meta("critical_value_W"), 1.93, 0.0051)
require_eq("Discussion: fixed-N critical final wealth ('about 1.9')", get_meta("critical_value_W"), 1.9, 0.051)
require_eq("ESM S20 betting patients (2,689)", get_meta("n_bet"), 2689, 0)
# "the supremum recovers about half" of the gap between the terminal and the fixed-N test
require_eq("ESM S20 supremum recovers about half",
           (get_prem("anytime crossing sup W") - get_prem("final wealth W_N >= 20")) /
             (get_prem("fixed-N final-wealth test") - get_prem("final wealth W_N >= 20")), 0.5, 0.05)
require_true("ESM S20 calibration seed 20260554 recorded",
             any(prem_meta == "calibration_seed,20260554"))
# "defaults n_0 = 30 and n_r = 50 events" (e-RTe kernel)
require_true("Construction: e-RTe kernel defaults burn_in = 30, ramp = 50",
             any(grepl("burn_in = 30, ramp = 50", readLines(local_ert_file("erte.R"), warn = FALSE), fixed = TRUE)))

cat("\n=== Out-of-order updates: which p_i prices the wager (round 2; not tabulated in the article) ===\n")

pi_chk <- read.csv(file.path(tables_dir, "pi_out_of_order_check.csv"), stringsAsFactors = FALSE)
require_eq("p_i check: revealed-count price, null crossing (1.6%)",
           pi_chk$null_crossing[pi_chk$bettor == "revealed_count"], 0.016, half)
require_eq("p_i check: logged-design price, null crossing (2.3%)",
           pi_chk$null_crossing[pi_chk$bettor == "logged_design"], 0.023, half)
require_eq("p_i check: adversary against the logged price (100%)",
           pi_chk$null_crossing[pi_chk$bettor == "adversary_vs_logged"], 1, 0)
require_eq("p_i check: adversary median first crossing (61 updates)",
           pi_chk$median_first_crossing[pi_chk$bettor == "adversary_vs_logged"], 61, 0)

cat("\n=== Summary ===\n")
if (length(skips) > 0L) {
  cat(sprintf("SKIPPED (tracked, not failing): %d\n", length(skips)))
  for (msg in skips) cat("  ~ ", msg, "\n", sep = "")
}
if (length(failures) == 0L) {
  cat("All pinned numbers reproduce within tolerance.\n")
  quit(status = 0)
} else {
  cat(sprintf("FAILED: %d assertions did not match.\n", length(failures)))
  for (msg in failures) cat("  - ", msg, "\n", sep = "")
  quit(status = 1)
}
