# =============================================================================
# e-RT Simulations: Type I Error and Power (Unified)
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# Core e-RT function
# -----------------------------------------------------------------------------

binary_design_lambdas <- function(p_trt, p_ctrl, p_random = 0.5) {
  event_denom <- p_random * p_trt + (1 - p_random) * p_ctrl
  nonevent_denom <- p_random * (1 - p_trt) + (1 - p_random) * (1 - p_ctrl)

  if (event_denom <= 0 || nonevent_denom <= 0) {
    stop("Design rates imply a degenerate event or non-event denominator")
  }

  list(
    event = p_random * p_trt / event_denom,
    nonevent = p_random * (1 - p_trt) / nonevent_denom
  )
}

compute_eRT <- function(treatment, outcome, p = 0.5, burn_in = 50, ramp = 100,
                        wager = c("adaptive", "design", "fixed", "oracle"),
                        p_trt_design = NULL, p_ctrl_design = NULL,
                        p_trt_oracle = NULL, p_ctrl_oracle = NULL,
                        fixed_event_lambda = NULL,
                        fixed_nonevent_lambda = NULL,
                        ramp_fixed = TRUE,
                        kelly_fraction = 0.5) {
  n <- length(treatment)
  wager <- match.arg(wager)

  if (is.factor(treatment)) treatment <- as.numeric(treatment) - 1
  if (is.factor(outcome)) outcome <- as.numeric(outcome) - 1

  if (wager == "design" &&
      (is.null(p_trt_design) || is.null(p_ctrl_design))) {
    stop("Design wager requires p_trt_design and p_ctrl_design")
  }

  if (wager == "oracle" &&
      (is.null(p_trt_oracle) || is.null(p_ctrl_oracle))) {
    stop("Oracle wager requires p_trt_oracle and p_ctrl_oracle")
  }

  if (wager == "fixed" &&
      (is.null(fixed_event_lambda) || is.null(fixed_nonevent_lambda))) {
    stop("Fixed wager requires fixed_event_lambda and fixed_nonevent_lambda")
  }

  if (wager == "design") {
    design_lambdas <- binary_design_lambdas(p_trt_design, p_ctrl_design, p)
  } else if (wager == "oracle") {
    design_lambdas <- binary_design_lambdas(p_trt_oracle, p_ctrl_oracle, p)
  } else if (wager == "fixed") {
    design_lambdas <- list(event = fixed_event_lambda, nonevent = fixed_nonevent_lambda)
  } else {
    design_lambdas <- NULL
  }

  wealth <- numeric(n)
  wealth[1] <- 1

  n_trt <- as.integer(treatment[1] == 1)
  n_ctrl <- as.integer(treatment[1] == 0)
  events_trt <- if (treatment[1] == 1) outcome[1] else 0
  events_ctrl <- if (treatment[1] == 0) outcome[1] else 0

  for (i in 2:n) {
    c_i <- max(0, min(1, (i - burn_in) / ramp))

    if (wager == "adaptive") {
      rate_trt <- if (n_trt > 0) events_trt / n_trt else 0.5
      rate_ctrl <- if (n_ctrl > 0) events_ctrl / n_ctrl else 0.5

      delta_hat <- rate_trt - rate_ctrl

      if (outcome[i] == 1) {
        lambda <- 0.5 + kelly_fraction * c_i * delta_hat
      } else {
        lambda <- 0.5 - kelly_fraction * c_i * delta_hat
      }
    } else {
      lambda_target <- if (outcome[i] == 1) {
        design_lambdas$event
      } else {
        design_lambdas$nonevent
      }

      if (ramp_fixed) {
        lambda <- p + c_i * (lambda_target - p)
      } else {
        lambda <- lambda_target
      }
    }

    lambda <- max(0.001, min(0.999, lambda))

    if (treatment[i] == 1) {
      multiplier <- lambda / p
    } else {
      multiplier <- (1 - lambda) / (1 - p)
    }

    wealth[i] <- wealth[i-1] * multiplier

    if (treatment[i] == 1) {
      n_trt <- n_trt + 1
      events_trt <- events_trt + outcome[i]
    } else {
      n_ctrl <- n_ctrl + 1
      events_ctrl <- events_ctrl + outcome[i]
    }
  }

  return(wealth)
}

# -----------------------------------------------------------------------------
# Simulate a single trial
# -----------------------------------------------------------------------------

simulate_trial <- function(n, p_trt = 0.5, rate_trt, rate_ctrl) {
  treatment <- rbinom(n, 1, p_trt)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, rate_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, rate_ctrl)
  return(data.frame(treatment = treatment, outcome = outcome))
}

# -----------------------------------------------------------------------------
# Unified simulation function
# -----------------------------------------------------------------------------

simulate_eRT <- function(n_sims = 5000,
                         p_ctrl,
                         p_trt = NULL,
                         hypothesized_ARR = NULL,
                         target_power = 0.80,
                         alpha = 0.05,
                         burn_in = 50,
                         ramp = 100,
                         wager = c("adaptive", "design", "fixed", "oracle"),
                         p_trt_design = NULL,
                         p_ctrl_design = NULL,
                         fixed_event_lambda = NULL,
                         fixed_nonevent_lambda = NULL,
                         ramp_fixed = TRUE,
                         kelly_fraction = 0.5) {
  wager <- match.arg(wager)

  # Determine hypothesized ARR for sample size calculation
  if (is.null(hypothesized_ARR) && is.null(p_trt)) {
    stop("Must specify either p_trt or hypothesized_ARR")
  }

  if (is.null(hypothesized_ARR)) {
    hypothesized_ARR <- p_ctrl - p_trt
  }

  # Calculate sample size based on hypothesized effect
  ss <- power.prop.test(p1 = p_ctrl, p2 = p_ctrl - hypothesized_ARR,
                        power = target_power, sig.level = alpha)
  n_per_arm <- ceiling(ss$n)
  n_patients <- 2 * n_per_arm

  # Determine true p_trt (null if not specified)
  if (is.null(p_trt)) {
    p_trt <- p_ctrl  # Null is true
    true_ARR <- 0
    sim_type <- "Type I Error"
  } else {
    true_ARR <- p_ctrl - p_trt
    sim_type <- "Power"
  }

  cat(sprintf("%s Simulation\n", sim_type))
  cat(sprintf("  Design: p_ctrl = %.2f, hypothesized ARR = %.2f, target power = %.0f%%\n",
              p_ctrl, hypothesized_ARR, target_power * 100))
  cat(sprintf("  Sample size: n = %d per arm = %d total\n", n_per_arm, n_patients))
  cat(sprintf("  Truth: p_ctrl = %.2f, p_trt = %.2f, true ARR = %.2f\n",
              p_ctrl, p_trt, true_ARR))
  cat(sprintf("  Wager policy: %s\n", wager))
  cat(sprintf("  n_sims = %d\n", n_sims))

  rejections <- 0
  first_crossing <- numeric(n_sims)
  final_evalues <- numeric(n_sims)

  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

  for (sim in 1:n_sims) {
    trial <- simulate_trial(n_patients, rate_trt = p_trt, rate_ctrl = p_ctrl)
    wealth <- compute_eRT(
      trial$treatment,
      trial$outcome,
      burn_in = burn_in,
      ramp = ramp,
      wager = wager,
      p_trt_design = p_trt_design,
      p_ctrl_design = p_ctrl_design,
      p_trt_oracle = p_trt,
      p_ctrl_oracle = p_ctrl,
      fixed_event_lambda = fixed_event_lambda,
      fixed_nonevent_lambda = fixed_nonevent_lambda,
      ramp_fixed = ramp_fixed,
      kelly_fraction = kelly_fraction
    )

    final_evalues[sim] <- wealth[n_patients]

    crossing <- which(wealth >= 1/alpha)
    if (length(crossing) > 0) {
      rejections <- rejections + 1
      first_crossing[sim] <- crossing[1]
    } else {
      first_crossing[sim] <- NA
    }

    setTxtProgressBar(pb, sim)
  }
  close(pb)

  rejection_rate <- rejections / n_sims
  se <- sqrt(rejection_rate * (1 - rejection_rate) / n_sims)
  median_stop <- median(first_crossing, na.rm = TRUE)

  cat(sprintf("\nResults:\n"))
  cat(sprintf("  Rejection rate: %.3f (SE: %.3f)\n", rejection_rate, se))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n", rejection_rate - 1.96*se, rejection_rate + 1.96*se))
  if (sim_type == "Power") {
    cat(sprintf("  Target power: %.0f%%\n", target_power * 100))
    cat(sprintf("  Median crossing: %.0f patients\n", median_stop))
  } else {
    cat(sprintf("  Nominal alpha: %.3f\n", alpha))
  }

  return(list(
    sim_type = sim_type,
    rejection_rate = rejection_rate,
    se = se,
    n_per_arm = n_per_arm,
    n_patients = n_patients,
    p_ctrl = p_ctrl,
    p_trt = p_trt,
    hypothesized_ARR = hypothesized_ARR,
    true_ARR = true_ARR,
    wager = wager,
    target_power = target_power,
    alpha = alpha,
    median_crossing = median_stop,
    first_crossing = first_crossing,
    final_evalues = final_evalues
  ))
}

# -----------------------------------------------------------------------------
# Plot trajectories
# -----------------------------------------------------------------------------

plot_trajectories <- function(n_trials = 30,
                              p_ctrl,
                              p_trt = NULL,
                              hypothesized_ARR = NULL,
                              target_power = 0.80,
                              alpha = 0.05,
                              burn_in = 50,
                              ramp = 100,
                              wager = c("adaptive", "design", "fixed", "oracle"),
                              p_trt_design = NULL,
                              p_ctrl_design = NULL,
                              fixed_event_lambda = NULL,
                              fixed_nonevent_lambda = NULL,
                              ramp_fixed = TRUE,
                              title = NULL) {
  wager <- match.arg(wager)

  if (is.null(hypothesized_ARR) && is.null(p_trt)) {
    stop("Must specify either p_trt or hypothesized_ARR")
  }

  if (is.null(hypothesized_ARR)) {
    hypothesized_ARR <- p_ctrl - p_trt
  }

  if (is.null(p_trt)) {
    p_trt <- p_ctrl  # Null
  }

  # Calculate sample size
  ss <- power.prop.test(p1 = p_ctrl, p2 = p_ctrl - hypothesized_ARR,
                        power = target_power, sig.level = alpha)
  n_patients <- 2 * ceiling(ss$n)

  true_ARR <- p_ctrl - p_trt

  if (is.null(title)) {
    if (true_ARR == 0) {
      title <- sprintf("e-RTb Under Null (n = %d, designed for %.0f%% ARR)",
                       n_patients, hypothesized_ARR * 100)
    } else {
      title <- sprintf("e-RTb Under Alternative (n = %d, %.0f%% power)",
                       n_patients, target_power * 100)
    }
  }

  trajectories <- list()

  for (i in 1:n_trials) {
    trial <- simulate_trial(n_patients, rate_trt = p_trt, rate_ctrl = p_ctrl)
    wealth <- compute_eRT(
      trial$treatment,
      trial$outcome,
      burn_in = burn_in,
      ramp = ramp,
      wager = wager,
      p_trt_design = p_trt_design,
      p_ctrl_design = p_ctrl_design,
      p_trt_oracle = p_trt,
      p_ctrl_oracle = p_ctrl,
      fixed_event_lambda = fixed_event_lambda,
      fixed_nonevent_lambda = fixed_nonevent_lambda,
      ramp_fixed = ramp_fixed
    )
    trajectories[[i]] <- data.frame(
      patient = 1:n_patients,
      wealth = wealth,
      trial = i
    )
  }

  df <- bind_rows(trajectories)

  p <- ggplot(df, aes(x = patient, y = wealth, group = trial)) +
    geom_line(alpha = 0.4, color = "steelblue") +
    geom_hline(yintercept = 1/alpha, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray50") +
    scale_y_log10() +
    labs(
      title = title,
      subtitle = sprintf("Control = %.0f%%, Treatment = %.0f%%, True ARR = %.0f%%",
                         p_ctrl * 100, p_trt * 100, true_ARR * 100),
      x = "Patient",
      y = "e-value (log scale)"
    ) +
    annotate("text", x = n_patients * 0.95, y = 1/alpha * 1.5,
             label = sprintf("1/alpha = %.0f", 1/alpha), color = "red", hjust = 1) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

# -----------------------------------------------------------------------------
# Run all simulations
# -----------------------------------------------------------------------------

run_all <- function(n_sims = 5000,
                    p_ctrl = 0.40,
                    ARRs = c(0.05, 0.10),
                    target_powers = c(0.80, 0.90),
                    alpha = 0.05,
                    seed = 42) {

  set.seed(seed)
  results <- list()

  # Create grid of scenarios
  scenarios <- expand.grid(
    hypothesized_ARR = ARRs,
    target_power = target_powers,
    stringsAsFactors = FALSE
  )

  # For each scenario, run Type I error and Power
  all_results <- data.frame()

  for (i in 1:nrow(scenarios)) {
    hyp_ARR <- scenarios$hypothesized_ARR[i]
    tgt_pow <- scenarios$target_power[i]

    cat(sprintf("\n========================================\n"))
    cat(sprintf("Scenario: ARR = %.0f%%, Target Power = %.0f%%\n",
                hyp_ARR * 100, tgt_pow * 100))
    cat(sprintf("========================================\n\n"))

    # Type I error (null true, same sample size)
    cat("--- Type I Error ---\n")
    t1 <- simulate_eRT(
      n_sims = n_sims,
      p_ctrl = p_ctrl,
      p_trt = NULL,  # Null is true
      hypothesized_ARR = hyp_ARR,
      target_power = tgt_pow,
      alpha = alpha
    )

    # Power (alternative true)
    cat("\n--- Power ---\n")
    pow <- simulate_eRT(
      n_sims = n_sims,
      p_ctrl = p_ctrl,
      p_trt = p_ctrl - hyp_ARR,  # True effect = hypothesized
      target_power = tgt_pow,
      alpha = alpha
    )

    all_results <- rbind(all_results, data.frame(
      hypothesized_ARR = hyp_ARR,
      target_power = tgt_pow,
      n_patients = t1$n_patients,
      type1_error = t1$rejection_rate,
      type1_se = t1$se,
      eRT_power = pow$rejection_rate,
      power_se = pow$se,
      median_crossing = pow$median_crossing
    ))
  }

  results$summary <- all_results

  # Print summary
  cat("\n========================================\n")
  cat("SUMMARY\n")
  cat("========================================\n\n")

  print(all_results %>%
          mutate(
            hypothesized_ARR = sprintf("%.0f%%", hypothesized_ARR * 100),
            target_power = sprintf("%.0f%%", target_power * 100),
            type1_error = sprintf("%.3f", type1_error),
            eRT_power = sprintf("%.1f%%", eRT_power * 100)
          ) %>%
          select(hypothesized_ARR, target_power, n_patients,
                 type1_error, eRT_power, median_crossing))

  return(results)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------

# Guard added in this repository only. The supplementary file runs this block in any
# interactive session; here, set options(ert.run_demo = TRUE) before sourcing to run it.
if (interactive() && isTRUE(getOption("ert.run_demo", FALSE))) {

  results <- run_all(
    n_sims = 5000,
    p_ctrl = 0.40,
    ARRs = c(0.05, 0.10),
    target_powers = c(0.80, 0.90),
    alpha = 0.05
  )

  # Trajectory plots - Alternative
  p1 <- plot_trajectories(p_ctrl = 0.40, p_trt = 0.30, target_power = 0.80)
  ggsave("traj_alt_10pct_80pow.pdf", p1, width = 8, height = 5)

  p2 <- plot_trajectories(p_ctrl = 0.40, p_trt = 0.30, target_power = 0.90)
  ggsave("traj_alt_10pct_90pow.pdf", p2, width = 8, height = 5)

  # Trajectory plots - Null
  p3 <- plot_trajectories(p_ctrl = 0.40, p_trt = NULL, hypothesized_ARR = 0.10, target_power = 0.80)
  ggsave("traj_null_10pct_80pow.pdf", p3, width = 8, height = 5)

  p4 <- plot_trajectories(p_ctrl = 0.40, p_trt = NULL, hypothesized_ARR = 0.10, target_power = 0.90)
  ggsave("traj_null_10pct_90pow.pdf", p4, width = 8, height = 5)
}
