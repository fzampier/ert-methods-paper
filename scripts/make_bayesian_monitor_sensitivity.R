# Bayesian scheduled-monitoring sensitivity (Figure 6 inputs). This is a comparator table only; it does not alter the e-RT
# interruption-curve figure or Table 1.

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260514)

n_sims <- 5000
n_patients <- 2690
p_ctrl <- 0.35
p_trt_alt <- 0.30
p_trt_null <- 0.35
look_fracs <- c(1 / 3, 2 / 3, 1)
look_n <- round(n_patients * look_fracs)
posterior_thresholds <- c(0.975, 0.98, 0.985, 0.99, 0.995, 0.9954, 0.999)

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f%%"), 100 * x))
}

fmt_num <- function(x, digits = 0) {
  ifelse(is.na(x), "--", format(round(x, digits), big.mark = ",", nsmall = digits, trim = TRUE))
}

write_md_table <- function(x, path) {
  x[] <- lapply(x, as.character)
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), path)
  invisible(path)
}

simulate_trial <- function(n, p_trt, p_ctrl, p_random = 0.5) {
  treatment <- rbinom(n, 1, p_random)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, p_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, p_ctrl)
  list(treatment = treatment, outcome = outcome)
}

arr_at <- function(treatment, outcome, n_look) {
  idx <- seq_len(n_look)
  trt <- treatment[idx] == 1
  ctrl <- treatment[idx] == 0
  if (!any(trt) || !any(ctrl)) return(NA_real_)
  mean(outcome[idx][ctrl]) - mean(outcome[idx][trt])
}

beta_gt <- function(a_x, b_x, a_y, b_y) {
  # Exact P[X > Y] for independent beta variables with integer shape a_x.
  i <- 0:(a_x - 1)
  sum(exp(
    lbeta(a_y + i, b_x + b_y) -
      log(b_x + i) -
      lbeta(1 + i, b_x) -
      lbeta(a_y, b_y)
  ))
}

posterior_prob_benefit <- function(treatment, outcome, n_look) {
  idx <- seq_len(n_look)
  trt <- treatment[idx] == 1
  ctrl <- treatment[idx] == 0

  y_trt <- sum(outcome[idx][trt])
  y_ctrl <- sum(outcome[idx][ctrl])
  n_trt <- sum(trt)
  n_ctrl <- sum(ctrl)

  # Independent beta(1, 1) priors. Benefit is a lower event risk in treatment.
  a_trt <- 1L + y_trt
  b_trt <- 1L + n_trt - y_trt
  a_ctrl <- 1L + y_ctrl
  b_ctrl <- 1L + n_ctrl - y_ctrl

  beta_gt(a_ctrl, b_ctrl, a_trt, b_trt)
}

bayes_first_crossing <- function(treatment, outcome, threshold) {
  for (n_look in look_n) {
    prob_benefit <- posterior_prob_benefit(treatment, outcome, n_look)
    if (is.finite(prob_benefit) && prob_benefit >= threshold) {
      return(list(
        first_crossing = n_look,
        statistic = prob_benefit
      ))
    }
  }
  list(first_crossing = NA_integer_, statistic = NA_real_)
}

bayes_look_probabilities <- function(treatment, outcome) {
  data.frame(
    look = look_n,
    posterior_benefit = vapply(
      look_n,
      function(n_look) posterior_prob_benefit(treatment, outcome, n_look),
      numeric(1)
    )
  )
}

bayes_first_crossing_from_looks <- function(look_probs, threshold) {
  hit <- which(look_probs$posterior_benefit >= threshold)
  if (!length(hit)) {
    return(list(first_crossing = NA_integer_, statistic = NA_real_))
  }
  k <- hit[[1]]
  list(
    first_crossing = look_probs$look[[k]],
    statistic = look_probs$posterior_benefit[[k]]
  )
}

ert_first_crossing <- function(wealth, treatment, outcome) {
  cross <- which(wealth >= 20)
  if (!length(cross)) {
    return(list(first_crossing = NA_integer_, statistic = NA_real_))
  }
  fc <- cross[[1]]
  list(first_crossing = fc, statistic = arr_at(treatment, outcome, fc))
}

run_bayes <- function(p_trt, scenario, seed) {
  set.seed(seed)
  rows <- vector("list", n_sims * length(posterior_thresholds))
  combo_rows <- vector("list", n_sims * length(posterior_thresholds) * 2L)
  row_i <- 1L
  combo_i <- 1L

  for (sim in seq_len(n_sims)) {
    trial <- simulate_trial(n_patients, p_trt = p_trt, p_ctrl = p_ctrl)
    final_effect <- arr_at(trial$treatment, trial$outcome, n_patients)
    bayes_looks <- bayes_look_probabilities(trial$treatment, trial$outcome)

    wealth_adaptive <- compute_eRT(
      trial$treatment,
      trial$outcome,
      burn_in = 50,
      ramp = 100,
      wager = "adaptive",
      kelly_fraction = 0.5
    )
    ert_adaptive <- ert_first_crossing(wealth_adaptive, trial$treatment, trial$outcome)

    wealth_design <- compute_eRT(
      trial$treatment,
      trial$outcome,
      wager = "design",
      p_trt_design = p_trt_alt,
      p_ctrl_design = p_ctrl,
      ramp_fixed = FALSE
    )
    ert_design <- ert_first_crossing(wealth_design, trial$treatment, trial$outcome)

    for (threshold in posterior_thresholds) {
      res <- bayes_first_crossing_from_looks(bayes_looks, threshold)
      effect_at_crossing <- if (!is.na(res$first_crossing)) {
        arr_at(trial$treatment, trial$outcome, res$first_crossing)
      } else {
        NA_real_
      }

      rows[[row_i]] <- data.frame(
        sim = sim,
        scenario = scenario,
        threshold = threshold,
        method = sprintf("Bayes scheduled Pr(benefit)>%.4f", threshold),
        first_crossing = res$first_crossing,
        posterior_at_crossing = res$statistic,
        effect_at_crossing = effect_at_crossing,
        final_effect = final_effect,
        type_m = abs(effect_at_crossing) / abs(final_effect)
      )
      row_i <- row_i + 1L

      for (ert_name in c("e-RTb adaptive", "e-RTb design directional")) {
        ert_res <- if (ert_name == "e-RTb adaptive") ert_adaptive else ert_design
        combo_first <- min(res$first_crossing, ert_res$first_crossing, na.rm = TRUE)
        if (!is.finite(combo_first)) combo_first <- NA_integer_

        combo_rows[[combo_i]] <- data.frame(
          sim = sim,
          scenario = scenario,
          threshold = threshold,
          bayes_method = sprintf("Bayes scheduled Pr(benefit)>%.4f", threshold),
          ert_method = ert_name,
          bayes_first_crossing = res$first_crossing,
          ert_first_crossing = ert_res$first_crossing,
          combined_first_crossing = combo_first,
          bayes_only = !is.na(res$first_crossing) && is.na(ert_res$first_crossing),
          ert_only = is.na(res$first_crossing) && !is.na(ert_res$first_crossing),
          both = !is.na(res$first_crossing) && !is.na(ert_res$first_crossing),
          neither = is.na(res$first_crossing) && is.na(ert_res$first_crossing)
        )
        combo_i <- combo_i + 1L
      }
    }
  }

  list(
    bayes = do.call(rbind, rows),
    combo = do.call(rbind, combo_rows)
  )
}

alternative <- run_bayes(p_trt_alt, "Alternative: 35% vs 30%", seed = 20260552)
null <- run_bayes(p_trt_null, "Null: 35% vs 35%", seed = 20260553)
bayes_rows <- rbind(alternative$bayes, null$bayes)
combo_rows <- rbind(alternative$combo, null$combo)

write.csv(
  bayes_rows,
  file.path(table_dir, "bayesian_monitor_first_crossings.csv"),
  row.names = FALSE
)

write.csv(
  combo_rows,
  file.path(table_dir, "bayesian_ert_union_first_crossings.csv"),
  row.names = FALSE
)

summarise_group <- function(rows) {
  crossed <- !is.na(rows$first_crossing)
  data.frame(
    interrupted_any = mean(crossed),
    median_interruption = median(rows$first_crossing[crossed], na.rm = TRUE),
    median_type_m = median(rows$type_m[crossed], na.rm = TRUE)
  )
}

split_rows <- split(bayes_rows, list(bayes_rows$scenario, bayes_rows$threshold), drop = TRUE)
summary_list <- lapply(names(split_rows), function(key) {
  rows <- split_rows[[key]]
  cbind(
    scenario = rows$scenario[[1]],
    threshold = rows$threshold[[1]],
    method = rows$method[[1]],
    summarise_group(rows)
  )
})
summary_df <- do.call(rbind, summary_list)
summary_df$threshold <- as.numeric(summary_df$threshold)
summary_df$interrupted_any <- as.numeric(summary_df$interrupted_any)
summary_df$median_interruption <- as.numeric(summary_df$median_interruption)
summary_df$median_type_m <- as.numeric(summary_df$median_type_m)

write.csv(
  summary_df,
  file.path(table_dir, "bayesian_monitor_method_summary.csv"),
  row.names = FALSE
)

alt_summary <- summary_df[summary_df$scenario == "Alternative: 35% vs 30%", ]
null_summary <- summary_df[summary_df$scenario == "Null: 35% vs 35%", ]
names(null_summary)[names(null_summary) == "interrupted_any"] <- "null_interruption"

table_df <- merge(
  alt_summary[, c("threshold", "method", "interrupted_any", "median_interruption", "median_type_m")],
  null_summary[, c("threshold", "null_interruption")],
  by = "threshold",
  all.x = TRUE
)
table_df <- table_df[order(table_df$threshold), ]

md_table <- data.frame(
  `Bayesian rule` = table_df$method,
  `Looks` = paste(format(look_n, big.mark = ","), collapse = ", "),
  `Alternative signal` = fmt_pct(table_df$interrupted_any),
  `Median signal` = fmt_num(table_df$median_interruption),
  `Null signal` = fmt_pct(table_df$null_interruption),
  `Median Type M` = ifelse(is.na(table_df$median_type_m), "--", sprintf("%.2fx", table_df$median_type_m)),
  check.names = FALSE
)

write_md_table(
  md_table,
  file.path(table_dir, "bayesian_monitor_sensitivity.md")
)

combo_split <- split(
  combo_rows,
  list(combo_rows$scenario, combo_rows$threshold, combo_rows$ert_method),
  drop = TRUE
)
combo_summary_list <- lapply(combo_split, function(rows) {
  crossed <- !is.na(rows$combined_first_crossing)
  data.frame(
    scenario = rows$scenario[[1]],
    threshold = rows$threshold[[1]],
    bayes_method = rows$bayes_method[[1]],
    ert_method = rows$ert_method[[1]],
    combined_signal = mean(crossed),
    median_combined_signal = median(rows$combined_first_crossing[crossed], na.rm = TRUE),
    bayes_only = mean(rows$bayes_only),
    ert_only = mean(rows$ert_only),
    both = mean(rows$both)
  )
})
combo_summary <- do.call(rbind, combo_summary_list)
combo_summary$threshold <- as.numeric(combo_summary$threshold)
combo_summary$combined_signal <- as.numeric(combo_summary$combined_signal)
combo_summary$median_combined_signal <- as.numeric(combo_summary$median_combined_signal)
combo_summary$bayes_only <- as.numeric(combo_summary$bayes_only)
combo_summary$ert_only <- as.numeric(combo_summary$ert_only)
combo_summary$both <- as.numeric(combo_summary$both)
combo_summary <- combo_summary[order(combo_summary$ert_method, combo_summary$threshold, combo_summary$scenario), ]

write.csv(
  combo_summary,
  file.path(table_dir, "bayesian_ert_union_summary.csv"),
  row.names = FALSE
)

combo_alt <- combo_summary[combo_summary$scenario == "Alternative: 35% vs 30%", ]
combo_null <- combo_summary[combo_summary$scenario == "Null: 35% vs 35%", ]
names(combo_null)[names(combo_null) == "combined_signal"] <- "null_combined_signal"

combo_table_df <- merge(
  combo_alt[, c("threshold", "bayes_method", "ert_method", "combined_signal", "median_combined_signal")],
  combo_null[, c("threshold", "ert_method", "null_combined_signal")],
  by = c("threshold", "ert_method"),
  all.x = TRUE
)
combo_table_df <- combo_table_df[order(combo_table_df$ert_method, combo_table_df$threshold), ]

combo_md_table <- data.frame(
  `Union rule` = paste(combo_table_df$bayes_method, "+", combo_table_df$ert_method),
  `Alternative signal` = fmt_pct(combo_table_df$combined_signal),
  `Median signal` = fmt_num(combo_table_df$median_combined_signal),
  `Null signal` = fmt_pct(combo_table_df$null_combined_signal),
  check.names = FALSE
)

write_md_table(
  combo_md_table,
  file.path(table_dir, "bayesian_ert_union_sensitivity.md")
)

cat("Wrote Bayesian monitoring sensitivity table.\n")
cat("Table: ", file.path(table_dir, "bayesian_monitor_sensitivity.md"), "\n", sep = "")
cat("Union table: ", file.path(table_dir, "bayesian_ert_union_sensitivity.md"), "\n", sep = "")
