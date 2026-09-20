# e-RTc: Continuous-Outcome Sequential Monitoring
#
# Adaptive mode preserves the standalone V7 script's effect-size-agnostic
# sign-direction wager. Design/oracle mode uses a prespecified normal-shift
# working model to choose the wager, while validity still follows from
# randomized treatment assignment.

continuous_design_lambda <- function(y, mu_ctrl, mu_trt, sd,
                                     p_random = 0.5) {
  if (!is.finite(sd) || sd <= 0) stop("sd must be positive")
  if (p_random <= 0 || p_random >= 1) stop("p_random must be in (0, 1)")

  log_trt <- log(p_random) + dnorm(y, mean = mu_trt, sd = sd, log = TRUE)
  log_ctrl <- log1p(-p_random) + dnorm(y, mean = mu_ctrl, sd = sd, log = TRUE)

  1 / (1 + exp(log_ctrl - log_trt))
}

cohens_d_estimate <- function(treatment, outcome, upto = length(treatment)) {
  if (is.na(upto) || upto < 1) return(NA_real_)
  upto <- min(length(treatment), as.integer(upto))

  idx <- seq_len(upto)
  trt <- treatment[idx] == 1
  ctrl <- treatment[idx] == 0
  if (sum(trt) < 2 || sum(ctrl) < 2) return(NA_real_)

  y_trt <- outcome[idx][trt]
  y_ctrl <- outcome[idx][ctrl]
  sd_trt <- sd(y_trt)
  sd_ctrl <- sd(y_ctrl)
  if (!is.finite(sd_trt) || !is.finite(sd_ctrl)) return(NA_real_)

  s_pooled <- sqrt((sd_trt^2 + sd_ctrl^2) / 2)
  if (!is.finite(s_pooled) || s_pooled <= 0) return(NA_real_)

  (mean(y_trt) - mean(y_ctrl)) / s_pooled
}

simulate_trial_continuous <- function(n, p_trt = 0.5, mu_ctrl = 0,
                                      mu_trt = 0, sd = 1) {
  treatment <- rbinom(n, 1, p_trt)
  outcome <- numeric(n)
  outcome[treatment == 0] <- rnorm(sum(treatment == 0), mean = mu_ctrl, sd = sd)
  outcome[treatment == 1] <- rnorm(sum(treatment == 1), mean = mu_trt, sd = sd)
  data.frame(treatment = treatment, outcome = outcome)
}

design_n_continuous <- function(design_d, target_power = 0.80,
                                alpha = 0.05, sd = 1) {
  if (!is.finite(design_d) || design_d == 0) {
    stop("design_d must be nonzero")
  }

  ss <- power.t.test(
    delta = abs(design_d) * sd,
    sd = sd,
    power = target_power,
    sig.level = alpha,
    type = "two.sample",
    alternative = "two.sided"
  )
  2 * ceiling(ss$n)
}

compute_eRTc <- function(treatment, outcome, p = 0.5, burn_in = 20,
                         ramp = 50, c_max = 0.6,
                         wager = c("adaptive", "design", "oracle"),
                         adaptive_direction = c("sign", "cohens_d"),
                         mu_ctrl_design = NULL, mu_trt_design = NULL,
                         sd_design = NULL,
                         mu_ctrl_oracle = NULL, mu_trt_oracle = NULL,
                         sd_oracle = NULL,
                         ramp_fixed = TRUE) {
  n <- length(treatment)
  if (length(outcome) != n) stop("treatment and outcome must have the same length")
  if (is.factor(treatment)) treatment <- as.numeric(treatment) - 1
  treatment <- as.numeric(treatment)
  wager <- match.arg(wager)
  adaptive_direction <- match.arg(adaptive_direction)

  if (wager == "design" &&
      (is.null(mu_ctrl_design) || is.null(mu_trt_design) || is.null(sd_design))) {
    stop("Design wager requires mu_ctrl_design, mu_trt_design, and sd_design")
  }
  if (wager == "oracle" &&
      (is.null(mu_ctrl_oracle) || is.null(mu_trt_oracle) || is.null(sd_oracle))) {
    stop("Oracle wager requires mu_ctrl_oracle, mu_trt_oracle, and sd_oracle")
  }

  wealth <- numeric(n)
  w <- 1

  for (i in seq_len(n)) {
    c_i <- min(1, max(0, (i - burn_in) / ramp))

    if (wager == "adaptive") {
      if ((i - 1) < burn_in) {
        lambda <- p
      } else {
        out_prev <- outcome[seq_len(i - 1)]
        trt_prev <- treatment[seq_len(i - 1)]

        med_prev <- median(out_prev, na.rm = TRUE)
        mad_prev <- mad(out_prev, center = med_prev, constant = 1, na.rm = TRUE)
        if (!is.finite(mad_prev) || mad_prev <= 0) mad_prev <- 1

        s_i <- (outcome[[i]] - med_prev) / mad_prev
        g_i <- s_i / (1 + abs(s_i))

        mean_trt <- mean(out_prev[trt_prev == 1])
        mean_ctrl <- mean(out_prev[trt_prev == 0])
        if (is.nan(mean_trt) || is.nan(mean_ctrl)) {
          direction <- 0
        } else if (adaptive_direction == "cohens_d") {
          d_hat <- cohens_d_estimate(trt_prev, out_prev)
          if (!is.finite(d_hat)) d_hat <- 0
          direction <- max(-1, min(1, d_hat))
        } else {
          direction <- sign(mean_trt - mean_ctrl)
        }

        lambda <- p + c_i * c_max * g_i * direction
      }
    } else {
      if (wager == "design") {
        target <- continuous_design_lambda(
          outcome[[i]],
          mu_ctrl = mu_ctrl_design,
          mu_trt = mu_trt_design,
          sd = sd_design,
          p_random = p
        )
      } else {
        target <- continuous_design_lambda(
          outcome[[i]],
          mu_ctrl = mu_ctrl_oracle,
          mu_trt = mu_trt_oracle,
          sd = sd_oracle,
          p_random = p
        )
      }
      lambda <- if (ramp_fixed) p + c_i * (target - p) else target
    }

    lambda <- max(0.001, min(0.999, lambda))
    multiplier <- if (treatment[[i]] == 1) lambda / p else (1 - lambda) / (1 - p)
    w <- w * multiplier
    wealth[[i]] <- w
  }

  wealth
}
