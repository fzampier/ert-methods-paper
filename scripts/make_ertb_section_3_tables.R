# Fresh e-RTb simulations for ESM Section 1 (Tables S1-S3). This script is
# stand-alone inside this repository and uses only the
# local R/ertb.R implementation.

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

set.seed(20260502)

p_ctrl <- 0.40
alpha <- 0.05
threshold <- 1 / alpha
design_arrs <- c(0.05, 0.10)
target_powers <- c(0.80, 0.90)
n_sims_main <- 5000
n_sims_tuning <- 1000
default_burn_in <- 50
default_ramp <- 100
default_kelly <- 0.5

fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
fmt_pp <- function(x, digits = 1) {
  ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f pp"), 100 * x))
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

design_arr_grid <- function(arr) {
  c(
    under = max(arr / 2, arr - 0.05),
    matched = arr,
    over = arr + 0.05
  )
}

policy_label <- function(policy) {
  switch(
    policy,
    adaptive = "Adaptive half-Kelly",
    under = "Design under",
    matched = "Design matched",
    over = "Design over",
    oracle = "Oracle",
    policy
  )
}

sample_size_total <- function(arr, target_power) {
  ss <- power.prop.test(
    p1 = p_ctrl,
    p2 = p_ctrl - arr,
    power = target_power,
    sig.level = alpha
  )
  2 * ceiling(ss$n)
}

simulate_trial_fast <- function(n, rate_trt, rate_ctrl) {
  treatment <- rbinom(n, 1, 0.5)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, rate_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, rate_ctrl)
  list(treatment = treatment, outcome = outcome)
}

compute_design_wealth <- function(treatment, outcome, p_trt_design, p_ctrl_design) {
  lambdas <- binary_design_lambdas(p_trt_design, p_ctrl_design, p_random = 0.5)
  lambda <- ifelse(outcome == 1, lambdas$event, lambdas$nonevent)
  multiplier <- ifelse(treatment == 1, lambda / 0.5, (1 - lambda) / 0.5)
  multiplier[[1L]] <- 1
  cumprod(multiplier)
}

arr_estimate <- function(treatment, outcome, upto = length(outcome)) {
  idx <- seq_len(upto)
  trt <- treatment[idx] == 1
  ctrl <- treatment[idx] == 0
  if (!any(trt) || !any(ctrl)) {
    return(NA_real_)
  }
  mean(outcome[idx][ctrl]) - mean(outcome[idx][trt])
}

summarise_policy <- function(rows, scenario_type, true_arr) {
  rejection_rate <- mean(rows$crossed)
  crossing_effects <- rows$crossing_effect[rows$crossed]
  data.frame(
    scenario_type = scenario_type,
    rejection_rate = rejection_rate,
    se = sqrt(rejection_rate * (1 - rejection_rate) / nrow(rows)),
    median_crossing = median(rows$first_crossing[rows$crossed], na.rm = TRUE),
    median_final_evalue = median(rows$final_evalue),
    median_crossing_effect = median(crossing_effects, na.rm = TRUE),
    median_final_effect = median(rows$final_effect, na.rm = TRUE),
    median_type_m = if (true_arr > 0) median(abs(crossing_effects) / true_arr, na.rm = TRUE) else NA_real_,
    type_m_q75 = if (true_arr > 0) quantile(abs(crossing_effects) / true_arr, 0.75, na.rm = TRUE) else NA_real_,
    type_m_q90 = if (true_arr > 0) quantile(abs(crossing_effects) / true_arr, 0.90, na.rm = TRUE) else NA_real_
  )
}

run_policy_grid <- function(arr, target_power, scenario_type, n_sims) {
  n <- sample_size_total(arr, target_power)
  true_arr <- if (scenario_type == "alternative") arr else 0
  p_trt_true <- p_ctrl - true_arr
  design_arrs_for_policy <- design_arr_grid(arr)

  policies <- c("adaptive", names(design_arrs_for_policy))
  if (scenario_type == "alternative") {
    policies <- c(policies, "oracle")
  }

  row_store <- vector("list", length(policies))
  names(row_store) <- policies
  for (policy in policies) {
    row_store[[policy]] <- data.frame(
      crossed = logical(0),
      first_crossing = integer(0),
      final_evalue = numeric(0),
      crossing_effect = numeric(0),
      final_effect = numeric(0)
    )
  }

  for (sim in seq_len(n_sims)) {
    trial <- simulate_trial_fast(n, p_trt_true, p_ctrl)
    treatment <- trial$treatment
    outcome <- trial$outcome
    final_effect <- arr_estimate(treatment, outcome)

    for (policy in policies) {
      wealth <- switch(
        policy,
        adaptive = compute_eRT(
          treatment,
          outcome,
          burn_in = default_burn_in,
          ramp = default_ramp,
          wager = "adaptive",
          kelly_fraction = default_kelly
        ),
        under = compute_design_wealth(
          treatment,
          outcome,
          p_trt_design = p_ctrl - design_arrs_for_policy[["under"]],
          p_ctrl_design = p_ctrl
        ),
        matched = compute_design_wealth(
          treatment,
          outcome,
          p_trt_design = p_ctrl - design_arrs_for_policy[["matched"]],
          p_ctrl_design = p_ctrl
        ),
        over = compute_design_wealth(
          treatment,
          outcome,
          p_trt_design = p_ctrl - design_arrs_for_policy[["over"]],
          p_ctrl_design = p_ctrl
        ),
        oracle = compute_design_wealth(
          treatment,
          outcome,
          p_trt_design = p_trt_true,
          p_ctrl_design = p_ctrl
        )
      )

      crossing <- which(wealth >= threshold)
      crossed <- length(crossing) > 0
      first_crossing <- if (crossed) crossing[[1L]] else NA_integer_
      crossing_effect <- if (crossed) arr_estimate(treatment, outcome, first_crossing) else NA_real_

      row_store[[policy]] <- rbind(
        row_store[[policy]],
        data.frame(
          crossed = crossed,
          first_crossing = first_crossing,
          final_evalue = wealth[[length(wealth)]],
          crossing_effect = crossing_effect,
          final_effect = final_effect
        )
      )
    }
  }

  summaries <- lapply(policies, function(policy) {
    wager_arr <- switch(
      policy,
      adaptive = NA_real_,
      under = design_arrs_for_policy[["under"]],
      matched = design_arrs_for_policy[["matched"]],
      over = design_arrs_for_policy[["over"]],
      oracle = true_arr
    )
    out <- summarise_policy(row_store[[policy]], scenario_type, true_arr)
    cbind(
      endpoint = "e-RTb",
      arr = arr,
      target_power = target_power,
      n_patients = n,
      policy = policy,
      policy_label = policy_label(policy),
      wager_arr = wager_arr,
      true_arr = true_arr,
      n_sims = n_sims,
      alpha = alpha,
      burn_in = ifelse(policy == "adaptive", default_burn_in, NA_integer_),
      ramp = ifelse(policy == "adaptive", default_ramp, NA_integer_),
      kelly_fraction = ifelse(policy == "adaptive", default_kelly, 1),
      out
    )
  })

  do.call(rbind, summaries)
}

message("Running fresh e-RTb operating-characteristic simulations...")
main_results <- do.call(rbind, lapply(design_arrs, function(arr) {
  do.call(rbind, lapply(target_powers, function(target_power) {
    message(sprintf("  ARR %.0f pp, target %.0f%%: null", 100 * arr, 100 * target_power))
    null_res <- run_policy_grid(arr, target_power, "null", n_sims_main)
    message(sprintf("  ARR %.0f pp, target %.0f%%: alternative", 100 * arr, 100 * target_power))
    alt_res <- run_policy_grid(arr, target_power, "alternative", n_sims_main)
    rbind(null_res, alt_res)
  }))
}))

run_tuning <- function(arr, target_power, n_sims) {
  n <- sample_size_total(arr, target_power)
  p_trt_true <- p_ctrl - arr
  settings <- data.frame(
    setting = c(
      "Default half-Kelly",
      "Short burn-in/ramp",
      "Proportional 5/10%",
      "Quarter-Kelly",
      "Three-quarter Kelly",
      "Full-Kelly"
    ),
    burn_in = c(default_burn_in, 30, round(0.05 * n), default_burn_in, default_burn_in, default_burn_in),
    ramp = c(default_ramp, 50, round(0.10 * n), default_ramp, default_ramp, default_ramp),
    kelly_fraction = c(0.5, 0.5, 0.5, 0.25, 0.75, 1.0),
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(settings))
  for (setting_row in seq_len(nrow(settings))) {
    crossed <- logical(n_sims)
    first_crossing <- rep(NA_integer_, n_sims)

    for (sim in seq_len(n_sims)) {
      trial <- simulate_trial_fast(n, p_trt_true, p_ctrl)
      wealth <- compute_eRT(
        trial$treatment,
        trial$outcome,
        burn_in = settings$burn_in[[setting_row]],
        ramp = settings$ramp[[setting_row]],
        wager = "adaptive",
        kelly_fraction = settings$kelly_fraction[[setting_row]]
      )
      crossing <- which(wealth >= threshold)
      crossed[[sim]] <- length(crossing) > 0
      if (crossed[[sim]]) {
        first_crossing[[sim]] <- crossing[[1L]]
      }
    }

    power <- mean(crossed)
    rows[[setting_row]] <- cbind(
      endpoint = "e-RTb",
      arr = arr,
      target_power = target_power,
      n_patients = n,
      n_sims = n_sims,
      alpha = alpha,
      settings[setting_row, ],
      power = power,
      se = sqrt(power * (1 - power) / n_sims),
      median_crossing = median(first_crossing[crossed], na.rm = TRUE)
    )
  }

  do.call(rbind, rows)
}

message("Running fresh e-RTb tuning sensitivity simulations...")
tuning_results <- do.call(rbind, lapply(design_arrs, function(arr) {
  message(sprintf("  ARR %.0f pp, target 80%%", 100 * arr))
  run_tuning(arr, 0.80, n_sims_tuning)
}))

write_csv(main_results, file.path(table_dir, "ertb_section3_fresh_operating_characteristics.csv"))
write_csv(tuning_results, file.path(table_dir, "ertb_section3_fresh_tuning.csv"))

type1_table <- main_results %>%
  filter(scenario_type == "null") %>%
  transmute(
    `Design ARR` = fmt_pp(arr),
    `Target power` = fmt_pct(target_power, digits = 0),
    `N` = fmt_num(n_patients),
    `Policy` = policy_label,
    `Wager ARR` = fmt_pp(wager_arr),
    `Type I error` = fmt_pct(rejection_rate),
    `SE` = fmt_pct(se, digits = 2)
  )

power_table <- main_results %>%
  filter(scenario_type == "alternative") %>%
  transmute(
    `True ARR` = fmt_pp(arr),
    `Target power` = fmt_pct(target_power, digits = 0),
    `N` = fmt_num(n_patients),
    `Policy` = policy_label,
    `Wager ARR` = fmt_pp(wager_arr),
    `Power` = fmt_pct(rejection_rate),
    `Median crossing` = fmt_num(median_crossing),
    `% N` = fmt_pct(median_crossing / n_patients, digits = 0),
    `ARR at crossing` = fmt_pp(median_crossing_effect),
    `Final ARR` = fmt_pp(median_final_effect),
    `Median Type M` = fmt_num(median_type_m, digits = 2)
  )

tuning_table <- tuning_results %>%
  transmute(
    `True ARR` = fmt_pp(arr),
    `N` = fmt_num(n_patients),
    `Setting` = setting,
    `Burn-in` = fmt_num(burn_in),
    `Ramp` = fmt_num(ramp),
    `Kelly fraction` = fmt_num(kelly_fraction, digits = 2),
    `Power` = fmt_pct(power),
    `Median crossing` = fmt_num(median_crossing)
  )

paths <- c(
  write_md_table(type1_table, file.path(table_dir, "ertb_section3_type1.md")),
  write_md_table(power_table, file.path(table_dir, "ertb_section3_power_type_m.md")),
  write_md_table(tuning_table, file.path(table_dir, "ertb_section3_tuning.md"))
)

message("Wrote:\n", paste(paths, collapse = "\n"))
