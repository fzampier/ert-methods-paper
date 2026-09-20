# e-RTe: Event-Only Sequential Monitoring
# Supplementary R Code for e-RT V7
# Date: 2026-03-08
#
# A simplified e-process that monitors only events (e.g., deaths, disease
# progression, myocardial infarction), without requiring non-event data or
# denominator tracking. Based on the insight that under 1:1 randomization,
# P(event from treatment | event occurred) = 0.5 under null.
#
# Key advantage: No non-event ascertainment needed; just log events.
# Active manuscript comparisons use explicit sample sizes, including same-N
# comparisons against e-RTb.

library(tidyverse)

# === CORE E-PROCESS ===

#' Compute the e-RTe wealth process from a stream of events.
#'
#' @param arms Character or numeric vector: 1 = treatment, 0 = control.
#'   Each element represents one event, in the order events occurred.
#' @param burn_in Number of events before betting begins (default 30).
#' @param ramp Number of events over which betting ramps from 0 to full (default 50).
#' @param threshold Rejection threshold, typically 1/alpha (default 20).
#' @param kelly_fraction Fractional adaptive Kelly intensity; 1 preserves the
#'   default full-intensity event-only wager.
#' @return A list with:
#'   - wealth: numeric vector of wealth after each event
#'   - crossed: logical, did wealth ever reach threshold?
#'   - crossed_at: event number at first crossing (NA if never)
#'   - final_p: final P(event from trt | event)
#'   - final_rr: final event rate ratio p/(1-p)
event_coin <- function(p_ctrl, p_trt, p_random = 0.5) {
  denom <- p_random * p_trt + (1 - p_random) * p_ctrl
  if (denom <= 0) stop("Event rates imply a degenerate event-coin denominator")
  p_random * p_trt / denom
}

compute_eRTe <- function(arms, burn_in = 30, ramp = 50, threshold = 20,
                         p_random = 0.5,
                         wager = c("adaptive", "design", "fixed", "oracle"),
                         p_trt_design = NULL, p_ctrl_design = NULL,
                         p_trt_oracle = NULL, p_ctrl_oracle = NULL,
                         fixed_lambda = NULL,
                         ramp_fixed = TRUE,
                         kelly_fraction = 1) {
  n <- length(arms)
  wager <- match.arg(wager)

  if (is.factor(arms)) arms <- as.numeric(arms) - 1
  arms <- as.numeric(arms)

  if (wager == "design" &&
      (is.null(p_trt_design) || is.null(p_ctrl_design))) {
    stop("Design wager requires p_trt_design and p_ctrl_design")
  }

  if (wager == "oracle" &&
      (is.null(p_trt_oracle) || is.null(p_ctrl_oracle))) {
    stop("Oracle wager requires p_trt_oracle and p_ctrl_oracle")
  }

  if (wager == "fixed" && is.null(fixed_lambda)) {
    stop("Fixed wager requires fixed_lambda on the event-coin scale")
  }

  fixed_target <- switch(
    wager,
    design = event_coin(p_ctrl_design, p_trt_design, p_random),
    oracle = event_coin(p_ctrl_oracle, p_trt_oracle, p_random),
    fixed = fixed_lambda,
    adaptive = NA_real_
  )

  wealth <- numeric(n)
  wealth[1] <- 1

  d_trt <- 0
  d_ctrl <- 0
  crossed <- FALSE
  crossed_at <- NA

  for (i in 1:n) {
    # Step 1: Compute p_obs from PREVIOUS events (before updating counts)
    total <- d_trt + d_ctrl
    p_obs <- if (total > 0) d_trt / total else 0.5

    # Step 2: Compute lambda (wager)
    c_i <- min(1, max(0, (i - burn_in) / ramp))
    if (wager == "adaptive") {
      if (i > burn_in && total > 0) {
        lambda <- p_random + kelly_fraction * c_i * (p_obs - p_random)
        lambda <- max(0.001, min(0.999, lambda))
      } else {
        lambda <- p_random  # Neutral bet during burn-in
      }
    } else {
      if (ramp_fixed) {
        lambda <- p_random + c_i * (fixed_target - p_random)
      } else {
        lambda <- fixed_target
      }
      lambda <- max(0.001, min(0.999, lambda))
    }

    # Step 3: Update wealth
    if (arms[i] == 1) {
      # Treatment event
      mult <- lambda / p_random
    } else {
      # Control event
      mult <- (1 - lambda) / (1 - p_random)
    }

    if (i == 1) {
      wealth[i] <- mult
    } else {
      wealth[i] <- wealth[i - 1] * mult
    }

    # Step 4: Update counts AFTER betting
    if (arms[i] == 1) {
      d_trt <- d_trt + 1
    } else {
      d_ctrl <- d_ctrl + 1
    }

    # Check threshold (record first crossing only)
    if (!crossed && wealth[i] >= threshold) {
      crossed <- TRUE
      crossed_at <- i
    }
  }

  # Final estimates
  total <- d_trt + d_ctrl
  final_p <- if (total > 0) d_trt / total else 0.5
  final_rr <- if (final_p > 0.001 && final_p < 0.999) {
    final_p / (1 - final_p)
  } else if (final_p <= 0.001) {
    0
  } else {
    Inf
  }

  list(
    wealth = wealth,
    crossed = crossed,
    crossed_at = crossed_at,
    final_p = final_p,
    final_rr = final_rr,
    d_trt = d_trt,
    d_ctrl = d_ctrl
  )
}

# === SIMULATION HELPERS ===

#' Generate a stream of events for simulation.
#'
#' @param n_events Number of events to simulate.
#' @param p_trt_given_event P(event from treatment | event occurred).
#'   Under null: 0.5. Under alternative with treatment benefit: < 0.5.
#' @return Numeric vector of arm labels (1 = treatment, 0 = control).
simulate_events <- function(n_events, p_trt_given_event) {
  rbinom(n_events, 1, p_trt_given_event)
}

#' Calculate expected number of events from event rates and sample size.
#'
#' @param n_patients Total number of patients (both arms).
#' @param p_ctrl Control arm event rate.
#' @param p_trt Treatment arm event rate.
#' @return Expected number of events.
expected_events <- function(n_patients, p_ctrl, p_trt) {
  ceiling(n_patients / 2 * (p_ctrl + p_trt))
}

# === FULL SIMULATION ===

#' Run e-RTe simulations under null and alternative hypotheses.
#'
#' @param n_sims Number of simulations (default 1000).
#' @param p_ctrl Control arm event rate.
#' @param p_trt Treatment arm event rate (under alternative).
#'   For null simulation, both arms use p_ctrl.
#' @param n_patients Total number of patients. If NULL, uses a legacy
#'   convenience default equal to the frequentist sample size with 2.5x
#'   inflation. Active manuscript scripts pass explicit sample sizes.
#' @param target_power Target power for sample size calculation (default 0.80).
#' @param alpha Significance level (default 0.05).
#' @param burn_in Burn-in period in events (default 30).
#' @param ramp Ramp period in events (default 50).
#' @return A list with null and alternative results.
simulate_eRTe <- function(n_sims = 1000, p_ctrl, p_trt,
                           n_patients = NULL, target_power = 0.80,
                           alpha = 0.05, burn_in = 30, ramp = 50,
                           wager = c("adaptive", "design", "fixed", "oracle"),
                           p_trt_design = NULL, p_ctrl_design = NULL,
                           fixed_lambda = NULL,
                           ramp_fixed = TRUE,
                           kelly_fraction = 1) {
  wager <- match.arg(wager)

  threshold <- 1 / alpha

  # Sample size calculation
  if (is.null(n_patients)) {
    ss <- power.prop.test(p1 = p_ctrl, p2 = p_trt,
                          power = target_power, sig.level = alpha)
    n_freq <- 2 * ceiling(ss$n)
    n_patients <- ceiling(n_freq * 2.5)  # 2.5x inflation for e-RTe
  } else {
    ss <- power.prop.test(p1 = p_ctrl, p2 = p_trt,
                          power = target_power, sig.level = alpha)
    n_freq <- 2 * ceiling(ss$n)
  }

  # Expected events under null and alternative
  n_events_null <- expected_events(n_patients, p_ctrl, p_ctrl)
  n_events_alt <- expected_events(n_patients, p_ctrl, p_trt)

  # Event coin probabilities
  p_coin_null <- 0.5
  p_coin_alt <- event_coin(p_ctrl, p_trt)

  cat(sprintf("e-RTe Simulation\n"))
  cat(sprintf("  Control event rate: %.1f%%\n", p_ctrl * 100))
  cat(sprintf("  Treatment event rate: %.1f%%\n", p_trt * 100))
  cat(sprintf("  ARR: %.1f%%\n", (p_ctrl - p_trt) * 100))
  cat(sprintf("  Event coin (null): %.3f\n", p_coin_null))
  cat(sprintf("  Event coin (alt):  %.3f (tilt = %.1f points)\n",
              p_coin_alt, abs(p_coin_alt - 0.5) * 100))
  cat(sprintf("  Wager policy: %s\n", wager))
  cat(sprintf("  Frequentist N: %d\n", n_freq))
  cat(sprintf("  e-RTe N: %d (%.1fx inflation)\n", n_patients, n_patients / n_freq))
  cat(sprintf("  Expected events (null): %d\n", n_events_null))
  cat(sprintf("  Expected events (alt): %d\n", n_events_alt))
  cat(sprintf("  Threshold: %.0f\n\n", threshold))

  # --- Null simulations ---
  cat("Running null simulations...\n")
  null_rejections <- 0
  null_first_crossing <- numeric(n_sims)
  null_final_evalues <- numeric(n_sims)

  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
  for (sim in 1:n_sims) {
    events <- simulate_events(n_events_null, p_coin_null)
    res <- compute_eRTe(
      events,
      burn_in = burn_in,
      ramp = ramp,
      threshold = threshold,
      wager = wager,
      p_trt_design = p_trt_design,
      p_ctrl_design = p_ctrl_design,
      p_trt_oracle = p_ctrl,
      p_ctrl_oracle = p_ctrl,
      fixed_lambda = fixed_lambda,
      ramp_fixed = ramp_fixed,
      kelly_fraction = kelly_fraction
    )
    null_final_evalues[sim] <- tail(res$wealth, 1)

    if (res$crossed) {
      null_rejections <- null_rejections + 1
      null_first_crossing[sim] <- res$crossed_at
    } else {
      null_first_crossing[sim] <- NA
    }
    setTxtProgressBar(pb, sim)
  }
  close(pb)

  type1 <- null_rejections / n_sims

  # --- Alternative simulations ---
  cat("Running alternative simulations...\n")
  alt_rejections <- 0
  alt_first_crossing <- numeric(n_sims)
  alt_final_evalues <- numeric(n_sims)

  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
  for (sim in 1:n_sims) {
    events <- simulate_events(n_events_alt, p_coin_alt)
    res <- compute_eRTe(
      events,
      burn_in = burn_in,
      ramp = ramp,
      threshold = threshold,
      wager = wager,
      p_trt_design = p_trt_design,
      p_ctrl_design = p_ctrl_design,
      p_trt_oracle = p_trt,
      p_ctrl_oracle = p_ctrl,
      fixed_lambda = fixed_lambda,
      ramp_fixed = ramp_fixed,
      kelly_fraction = kelly_fraction
    )
    alt_final_evalues[sim] <- tail(res$wealth, 1)

    if (res$crossed) {
      alt_rejections <- alt_rejections + 1
      alt_first_crossing[sim] <- res$crossed_at
    } else {
      alt_first_crossing[sim] <- NA
    }
    setTxtProgressBar(pb, sim)
  }
  close(pb)

  power <- alt_rejections / n_sims

  # --- Summary ---
  cat(sprintf("\n==========================================\n"))
  cat(sprintf("  RESULTS\n"))
  cat(sprintf("==========================================\n"))
  cat(sprintf("Type I Error: %.2f%% (SE: %.2f%%)\n",
              type1 * 100,
              sqrt(type1 * (1 - type1) / n_sims) * 100))
  cat(sprintf("Power:        %.1f%% (SE: %.1f%%)\n",
              power * 100,
              sqrt(power * (1 - power) / n_sims) * 100))
  if (sum(!is.na(alt_first_crossing)) > 0) {
    med_cross <- median(alt_first_crossing, na.rm = TRUE)
    cat(sprintf("Median crossing: %.0f events (%.0f%% of total)\n",
                med_cross, med_cross / n_events_alt * 100))
  }

  list(
    # Parameters
    p_ctrl = p_ctrl,
    p_trt = p_trt,
    n_patients = n_patients,
    n_freq = n_freq,
    n_events_null = n_events_null,
    n_events_alt = n_events_alt,
    p_coin_null = p_coin_null,
    p_coin_alt = p_coin_alt,
    wager = wager,
    threshold = threshold,
    # Null results
    type1 = type1,
    type1_se = sqrt(type1 * (1 - type1) / n_sims),
    null_final_evalues = null_final_evalues,
    null_first_crossing = null_first_crossing,
    # Alternative results
    power = power,
    power_se = sqrt(power * (1 - power) / n_sims),
    alt_final_evalues = alt_final_evalues,
    alt_first_crossing = alt_first_crossing,
    median_crossing = median(alt_first_crossing, na.rm = TRUE)
  )
}

# === TRAJECTORY VISUALIZATION ===

#' Plot e-RTe wealth trajectories under null and/or alternative.
#'
#' @param n_trials Number of trajectories to plot (default 30).
#' @param n_events Number of events per trajectory.
#' @param p_trt_given_event P(event from treatment | event occurred).
#' @param burn_in Burn-in period (default 30).
#' @param ramp Ramp period (default 50).
#' @param alpha Significance level (default 0.05).
#' @param title Plot title.
#' @return A ggplot2 object.
plot_trajectories_eRTe <- function(n_trials = 30, n_events,
                                    p_trt_given_event,
                                    burn_in = 30, ramp = 50,
                                    alpha = 0.05,
                                    wager = c("adaptive", "fixed"),
                                    fixed_lambda = NULL,
                                    ramp_fixed = TRUE,
                                    title = NULL) {
  wager <- match.arg(wager)

  threshold <- 1 / alpha

  if (is.null(title)) {
    title <- if (p_trt_given_event == 0.5) {
      sprintf("e-RTe Under Null (p = 0.50, %d events)", n_events)
    } else {
      sprintf("e-RTe Under Alternative (p = %.3f, %d events)", p_trt_given_event, n_events)
    }
  }

  trajectories <- lapply(1:n_trials, function(i) {
    events <- simulate_events(n_events, p_trt_given_event)
    res <- compute_eRTe(
      events,
      burn_in = burn_in,
      ramp = ramp,
      threshold = threshold,
      wager = wager,
      fixed_lambda = fixed_lambda,
      ramp_fixed = ramp_fixed
    )
    data.frame(event = 1:n_events, wealth = res$wealth, trial = i)
  })

  df <- bind_rows(trajectories)

  ggplot(df, aes(x = event, y = wealth, group = trial)) +
    geom_line(alpha = 0.4, color = "steelblue") +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "gray50") +
    scale_y_log10() +
    labs(title = title,
         subtitle = sprintf("P(trt event | event) = %.3f, threshold = %.0f",
                            p_trt_given_event, threshold),
         x = "Event Number", y = "e-value (log scale)") +
    annotate("text", x = n_events * 0.95, y = threshold * 1.5,
             label = sprintf("1/alpha = %.0f", threshold), color = "red", hjust = 1) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
}

# === SIGNAL CONCENTRATION ANALYSIS ===

#' Compute the event-coin tilt for different baseline event rates.
#'
#' Demonstrates the signal concentration property: when baseline event rate
#' is low, the event coin tilt (deviation from 0.5) exceeds the ARR,
#' concentrating the signal in the event stream.
#'
#' @param baseline_rates Vector of baseline (control) event rates.
#' @param arr Absolute risk reduction (default 0.05).
#' @return Data frame with baseline rate, treatment rate, event coin, and tilt.
signal_concentration_table <- function(baseline_rates = c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
                                        arr = 0.05) {
  data.frame(
    baseline = baseline_rates,
    trt_rate = baseline_rates - arr,
    event_coin = sapply(baseline_rates, function(p) event_coin(p, p - arr)),
    tilt_from_half = sapply(baseline_rates, function(p) abs(event_coin(p, p - arr) - 0.5)),
    arr = arr
  )
}

# === EXAMPLE USAGE ===
# Uncomment to run:

# # Signal concentration demonstration
# cat("\n=== Signal Concentration ===\n")
# sc <- signal_concentration_table(arr = 0.05)
# print(sc, digits = 3)
#
# # Quick simulation: 15% baseline event rate, 5pp ARR
# results <- simulate_eRTe(n_sims = 1000, p_ctrl = 0.15, p_trt = 0.10)
#
# # Plot trajectories
# p_alt <- event_coin(0.15, 0.10)
# n_events <- expected_events(results$n_patients, 0.15, 0.10)
#
# p_null <- plot_trajectories_eRTe(n_events = n_events, p_trt_given_event = 0.5)
# p_alt_plot <- plot_trajectories_eRTe(n_events = n_events, p_trt_given_event = p_alt)
#
# print(p_null)
# print(p_alt_plot)
