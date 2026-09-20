# FACTT adaptive-only permutation sensitivity figure and table.
#
# This script uses the ignored local validated BioLINCC FACTT CSV and writes
# only aggregate, manuscript-safe outputs. No shuffled patient-level datasets
# are saved.

suppressPackageStartupMessages({
  library(ggplot2)
  library(Rcpp)
})

source("scripts/paths.R")

suppressPackageStartupMessages(source(local_ert_file("ertb.R")))
suppressPackageStartupMessages(source(local_ert_file("erte.R")))
source(local_ert_file("ertc.R"))

OUT_DIR <- file.path("outputs", "factt_permutation")
MANUSCRIPT_DIR <- paper_file("figures")
TABLE_DIR <- paper_file("tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

N_PERM <- as.integer(Sys.getenv("FACTT_N_PERM", "10000"))
SEED <- 20260508L
THRESHOLD <- 20
EPS <- 1e-8

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
using namespace Rcpp;

double weighted_value_at_rank(const std::vector<int>& counts, int rank) {
  int cumulative = 0;
  for (int i = 0; i < static_cast<int>(counts.size()); ++i) {
    cumulative += counts[i];
    if (cumulative >= rank) return static_cast<double>(i);
  }
  return static_cast<double>(counts.size() - 1);
}

double weighted_median_counts(const std::vector<int>& counts, int total) {
  if (total <= 0) return NA_REAL;
  const int lo = (total + 1) / 2;
  const int hi = (total + 2) / 2;
  return 0.5 * (
    weighted_value_at_rank(counts, lo) +
    weighted_value_at_rank(counts, hi)
  );
}

double weighted_deviation_at_rank(
    const std::vector<std::pair<double, int> >& deviations,
    int rank) {
  int cumulative = 0;
  for (size_t i = 0; i < deviations.size(); ++i) {
    cumulative += deviations[i].second;
    if (cumulative >= rank) return deviations[i].first;
  }
  return deviations.back().first;
}

double weighted_mad_counts(
    const std::vector<int>& counts,
    int total,
    double center) {
  if (total <= 0) return NA_REAL;

  std::vector<std::pair<double, int> > deviations;
  deviations.reserve(counts.size());
  for (int i = 0; i < static_cast<int>(counts.size()); ++i) {
    if (counts[i] > 0) {
      deviations.push_back(
        std::make_pair(std::fabs(static_cast<double>(i) - center), counts[i])
      );
    }
  }

  std::sort(
    deviations.begin(),
    deviations.end(),
    [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
      return a.first < b.first;
    }
  );

  const int lo = (total + 1) / 2;
  const int hi = (total + 2) / 2;
  return 0.5 * (
    weighted_deviation_at_rank(deviations, lo) +
    weighted_deviation_at_rank(deviations, hi)
  );
}

double sample_variance(int n, double sum, double sumsq) {
  if (n < 2) return NA_REAL;
  double var = (sumsq - sum * sum / static_cast<double>(n)) /
    static_cast<double>(n - 1);
  if (var < 0 && var > -1e-12) var = 0;
  return var;
}

// [[Rcpp::export]]
NumericVector compute_eRTc_integer_adaptive_fast(
    NumericVector treatment,
    NumericVector outcome,
    int max_value,
    double p = 0.5,
    int burn_in = 50,
    int ramp = 100,
    double c_max = 0.6) {
  const int n = treatment.size();
  if (outcome.size() != n) stop("treatment and outcome must have the same length");
  if (max_value < 0) stop("max_value must be nonnegative");

  std::vector<int> counts(max_value + 1, 0);
  NumericVector wealth(n);

  double w = 1.0;
  int n_trt = 0;
  int n_ctrl = 0;
  double sum_trt = 0.0;
  double sum_ctrl = 0.0;
  double sumsq_trt = 0.0;
  double sumsq_ctrl = 0.0;

  for (int i = 0; i < n; ++i) {
    const double y = outcome[i];
    const int y_int = static_cast<int>(std::round(y));
    if (!R_finite(y) || std::fabs(y - static_cast<double>(y_int)) > 1e-8 ||
        y_int < 0 || y_int > max_value) {
      stop("integer bounded outcome required for fast FACTT e-RTc");
    }

    const double c_i = std::min(
      1.0,
      std::max(0.0, (static_cast<double>(i + 1) - burn_in) / ramp)
    );

    double lambda = p;
    if (i >= burn_in) {
      const double med_prev = weighted_median_counts(counts, i);
      double mad_prev = weighted_mad_counts(counts, i, med_prev);
      if (!R_finite(mad_prev) || mad_prev <= 0) mad_prev = 1.0;

      const double s_i = (y - med_prev) / mad_prev;
      const double g_i = s_i / (1.0 + std::fabs(s_i));

      double direction = 0.0;
      if (n_trt > 0 && n_ctrl > 0) {
        const double var_trt = sample_variance(n_trt, sum_trt, sumsq_trt);
        const double var_ctrl = sample_variance(n_ctrl, sum_ctrl, sumsq_ctrl);

        if (R_finite(var_trt) && R_finite(var_ctrl)) {
          const double pooled = std::sqrt((var_trt + var_ctrl) / 2.0);
          if (R_finite(pooled) && pooled > 0) {
            const double mean_trt = sum_trt / static_cast<double>(n_trt);
            const double mean_ctrl = sum_ctrl / static_cast<double>(n_ctrl);
            double d_hat = (mean_trt - mean_ctrl) / pooled;
            if (R_finite(d_hat)) {
              direction = std::max(-1.0, std::min(1.0, d_hat));
            }
          }
        }
      }

      lambda = p + c_i * c_max * g_i * direction;
    }

    lambda = std::max(0.001, std::min(0.999, lambda));
    if (treatment[i] == 1) {
      w *= lambda / p;
    } else {
      w *= (1.0 - lambda) / (1.0 - p);
    }
    wealth[i] = w;

    counts[y_int] += 1;
    if (treatment[i] == 1) {
      n_trt += 1;
      sum_trt += y;
      sumsq_trt += y * y;
    } else {
      n_ctrl += 1;
      sum_ctrl += y;
      sumsq_ctrl += y * y;
    }
  }

  return wealth;
}
')

fmt_pct <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x)
}

fmt_pp <- function(x, digits = 1) {
  ifelse(is.na(x), "", sprintf(paste0("%+.", digits, "f pp"), 100 * x))
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", format(round(x, digits), nsmall = digits, trim = TRUE))
}

fmt_int <- function(x) {
  ifelse(is.na(x), "", format(round(x), big.mark = ",", trim = TRUE))
}

fmt_effect <- function(effect, mean_diff = NA_real_, scale = c("d", "ARR")) {
  scale <- match.arg(scale)
  if (is.na(effect)) return("")
  if (scale == "d") {
    return(sprintf("d=%.2f (%+.1f VFD)", effect, mean_diff))
  }
  fmt_pp(effect, 1)
}

first_crossing <- function(wealth, threshold = THRESHOLD) {
  idx <- which(wealth >= threshold)
  if (length(idx) == 0L) NA_integer_ else idx[[1L]]
}

arr_estimate <- function(arm, outcome, upto = length(outcome)) {
  idx <- seq_len(min(upto, length(outcome)))
  trt <- outcome[idx][arm[idx] == 1L]
  ctrl <- outcome[idx][arm[idx] == 0L]
  mean(ctrl, na.rm = TRUE) - mean(trt, na.rm = TRUE)
}

mean_diff_estimate <- function(arm, outcome, upto = length(outcome)) {
  idx <- seq_len(min(upto, length(outcome)))
  trt <- outcome[idx][arm[idx] == 1L]
  ctrl <- outcome[idx][arm[idx] == 0L]
  mean(trt, na.rm = TRUE) - mean(ctrl, na.rm = TRUE)
}

d_estimate <- function(arm, outcome, upto = length(outcome)) {
  cohens_d_estimate(arm, outcome, upto = upto)
}

curve_band <- function(mat, panel, update_unit) {
  qs <- apply(mat, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  data.frame(
    panel = panel,
    index = seq_len(ncol(mat)),
    q025 = pmax(as.numeric(qs[1, ]), EPS),
    median = pmax(as.numeric(qs[2, ]), EPS),
    q975 = pmax(as.numeric(qs[3, ]), EPS),
    update_unit = update_unit,
    stringsAsFactors = FALSE
  )
}

observed_curve <- function(wealth, panel, update_unit) {
  data.frame(
    panel = panel,
    index = seq_along(wealth),
    wealth = pmax(as.numeric(wealth), EPS),
    update_unit = update_unit,
    stringsAsFactors = FALSE
  )
}

observed_label <- function(wealth, update_unit) {
  idx <- first_crossing(wealth)
  label <- switch(
    update_unit,
    patient = "masked index",
    event = "death-event index",
    "index"
  )
  if (is.na(idx)) {
    paste0("No crossing; peak E = ", fmt_num(max(wealth), 2))
  } else {
    paste0("Crossed at ", label, " ", fmt_int(idx))
  }
}

summarize_permutation <- function(label, update_unit, effect_scale,
                                  observed_wealth, crossings,
                                  crossing_effects, final_effect,
                                  crossing_mean_diffs = NULL,
                                  final_mean_diff = NA_real_) {
  crossed <- !is.na(crossings)
  effect_cross <- crossing_effects[crossed]
  mean_diff_cross <- if (is.null(crossing_mean_diffs)) {
    rep(NA_real_, length(effect_cross))
  } else {
    crossing_mean_diffs[crossed]
  }
  type_m <- abs(effect_cross) / abs(final_effect)

  data.frame(
    curve = label,
    update_unit = update_unit,
    effect_scale = effect_scale,
    observed = observed_label(observed_wealth, update_unit),
    crossing_rate = mean(crossed),
    median_crossing = if (any(crossed)) median(crossings[crossed]) else NA_real_,
    q025_crossing = if (any(crossed)) unname(stats::quantile(crossings[crossed], 0.025)) else NA_real_,
    q975_crossing = if (any(crossed)) unname(stats::quantile(crossings[crossed], 0.975)) else NA_real_,
    median_effect_at_crossing = if (length(effect_cross)) median(effect_cross) else NA_real_,
    median_mean_diff_at_crossing = if (length(mean_diff_cross)) median(mean_diff_cross) else NA_real_,
    final_effect = final_effect,
    final_mean_difference = final_mean_diff,
    median_type_m = if (length(type_m)) median(type_m) else NA_real_,
    q975_type_m = if (length(type_m)) unname(stats::quantile(type_m, 0.975)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

write_md_table <- function(x, path) {
  x[] <- lapply(x, as.character)
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), path)
}

escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

write_tex_table <- function(x, path) {
  x[] <- lapply(x, escape_tex)
  align <- paste0("l", paste(rep("r", ncol(x) - 1L), collapse = ""))
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{FACTT adaptive e-RT order-sensitivity summary from 10,000 complete random patient-order permutations. Crossing index summaries are among permutations that crossed. Type M is the absolute ratio of the effect at first crossing to the final effect.}",
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste(names(x), collapse = " & "),
    "\\\\ \\hline"
  )
  body <- apply(x, 1, function(row) paste0(paste(row, collapse = " & "), "\\\\"))
  lines <- c(lines, body, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
}

factt <- read.csv(validated_trial_csv("factt"), stringsAsFactors = FALSE)
# The validated enroll_order is the randomized/masked BioLINCC dataset order,
# not verified calendar enrollment order.
factt <- factt[order(factt$enroll_order), ]
factt$arm_fluid <- as.integer(factt$arm_fluid)
factt$arm_catheter <- as.integer(factt$arm_catheter)
factt$death_d60 <- as.integer(factt$death_d60)

max_vfd <- max(factt$vfd28, na.rm = TRUE)
if (any(abs(factt$vfd28 - round(factt$vfd28)) > 1e-8, na.rm = TRUE)) {
  stop("FACTT vfd28 must be integer-valued for the fast adaptive e-RTc path")
}

n <- nrow(factt)
death_events <- sum(factt$death_d60 == 1L)

panel_levels <- c(
  "A. Fluid VFD-28 e-RTc",
  "B. Fluid mortality e-RTb",
  "C. Fluid mortality e-RTe",
  "D. Catheter VFD-28 e-RTc",
  "E. Catheter mortality e-RTb",
  "F. Catheter mortality e-RTe"
)

obs_fluid_vfd <- compute_eRTc_integer_adaptive_fast(
  factt$arm_fluid, factt$vfd28,
  max_value = max_vfd, burn_in = 50, ramp = 100
)
obs_fluid_b <- compute_eRT(
  factt$arm_fluid, factt$death_d60,
  burn_in = 50, ramp = 100, wager = "adaptive"
)
obs_fluid_e <- compute_eRTe(
  factt$arm_fluid[factt$death_d60 == 1L],
  burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
)$wealth

obs_catheter_vfd <- compute_eRTc_integer_adaptive_fast(
  factt$arm_catheter, factt$vfd28,
  max_value = max_vfd, burn_in = 50, ramp = 100
)
obs_catheter_b <- compute_eRT(
  factt$arm_catheter, factt$death_d60,
  burn_in = 50, ramp = 100, wager = "adaptive"
)
obs_catheter_e <- compute_eRTe(
  factt$arm_catheter[factt$death_d60 == 1L],
  burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
)$wealth

check_fluid_vfd <- compute_eRTc(
  factt$arm_fluid,
  factt$vfd28,
  burn_in = 50,
  ramp = 100,
  wager = "adaptive",
  adaptive_direction = "cohens_d"
)
check_catheter_vfd <- compute_eRTc(
  factt$arm_catheter,
  factt$vfd28,
  burn_in = 50,
  ramp = 100,
  wager = "adaptive",
  adaptive_direction = "cohens_d"
)
stopifnot(max(abs(obs_fluid_vfd - check_fluid_vfd)) < 1e-10)
stopifnot(max(abs(obs_catheter_vfd - check_catheter_vfd)) < 1e-10)

mat_fluid_vfd <- matrix(NA_real_, nrow = N_PERM, ncol = n)
mat_fluid_b <- matrix(NA_real_, nrow = N_PERM, ncol = n)
mat_fluid_e <- matrix(NA_real_, nrow = N_PERM, ncol = death_events)
mat_catheter_vfd <- matrix(NA_real_, nrow = N_PERM, ncol = n)
mat_catheter_b <- matrix(NA_real_, nrow = N_PERM, ncol = n)
mat_catheter_e <- matrix(NA_real_, nrow = N_PERM, ncol = death_events)

cross_fluid_vfd <- rep(NA_integer_, N_PERM)
cross_fluid_b <- rep(NA_integer_, N_PERM)
cross_fluid_e <- rep(NA_integer_, N_PERM)
cross_catheter_vfd <- rep(NA_integer_, N_PERM)
cross_catheter_b <- rep(NA_integer_, N_PERM)
cross_catheter_e <- rep(NA_integer_, N_PERM)

effect_fluid_vfd <- rep(NA_real_, N_PERM)
effect_fluid_b <- rep(NA_real_, N_PERM)
effect_fluid_e <- rep(NA_real_, N_PERM)
effect_catheter_vfd <- rep(NA_real_, N_PERM)
effect_catheter_b <- rep(NA_real_, N_PERM)
effect_catheter_e <- rep(NA_real_, N_PERM)

mean_diff_fluid_vfd <- rep(NA_real_, N_PERM)
mean_diff_catheter_vfd <- rep(NA_real_, N_PERM)

set.seed(SEED)
message("Running ", N_PERM, " complete random FACTT patient-order permutations...")
pb <- txtProgressBar(min = 0, max = N_PERM, style = 3)

for (sim in seq_len(N_PERM)) {
  ord <- sample.int(n)
  fluid <- factt$arm_fluid[ord]
  catheter <- factt$arm_catheter[ord]
  death <- factt$death_d60[ord]
  vfd <- factt$vfd28[ord]

  w <- compute_eRTc_integer_adaptive_fast(
    fluid, vfd, max_value = max_vfd, burn_in = 50, ramp = 100
  )
  mat_fluid_vfd[sim, ] <- w
  cross_fluid_vfd[sim] <- first_crossing(w)
  if (!is.na(cross_fluid_vfd[sim])) {
    effect_fluid_vfd[sim] <- d_estimate(fluid, vfd, cross_fluid_vfd[sim])
    mean_diff_fluid_vfd[sim] <- mean_diff_estimate(fluid, vfd, cross_fluid_vfd[sim])
  }

  w <- compute_eRT(fluid, death, burn_in = 50, ramp = 100, wager = "adaptive")
  mat_fluid_b[sim, ] <- w
  cross_fluid_b[sim] <- first_crossing(w)
  if (!is.na(cross_fluid_b[sim])) {
    effect_fluid_b[sim] <- arr_estimate(fluid, death, cross_fluid_b[sim])
  }

  event_pos <- which(death == 1L)
  r <- compute_eRTe(
    fluid[event_pos],
    burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
  )$wealth
  mat_fluid_e[sim, ] <- r
  cross_fluid_e[sim] <- first_crossing(r)
  if (!is.na(cross_fluid_e[sim])) {
    effect_fluid_e[sim] <- arr_estimate(fluid, death, event_pos[cross_fluid_e[sim]])
  }

  w <- compute_eRTc_integer_adaptive_fast(
    catheter, vfd, max_value = max_vfd, burn_in = 50, ramp = 100
  )
  mat_catheter_vfd[sim, ] <- w
  cross_catheter_vfd[sim] <- first_crossing(w)
  if (!is.na(cross_catheter_vfd[sim])) {
    effect_catheter_vfd[sim] <- d_estimate(catheter, vfd, cross_catheter_vfd[sim])
    mean_diff_catheter_vfd[sim] <- mean_diff_estimate(catheter, vfd, cross_catheter_vfd[sim])
  }

  w <- compute_eRT(catheter, death, burn_in = 50, ramp = 100, wager = "adaptive")
  mat_catheter_b[sim, ] <- w
  cross_catheter_b[sim] <- first_crossing(w)
  if (!is.na(cross_catheter_b[sim])) {
    effect_catheter_b[sim] <- arr_estimate(catheter, death, cross_catheter_b[sim])
  }

  r <- compute_eRTe(
    catheter[event_pos],
    burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
  )$wealth
  mat_catheter_e[sim, ] <- r
  cross_catheter_e[sim] <- first_crossing(r)
  if (!is.na(cross_catheter_e[sim])) {
    effect_catheter_e[sim] <- arr_estimate(catheter, death, event_pos[cross_catheter_e[sim]])
  }

  if (sim %% 100L == 0L) setTxtProgressBar(pb, sim)
}
setTxtProgressBar(pb, N_PERM)
close(pb)

bands <- do.call(rbind, list(
  curve_band(mat_fluid_vfd, panel_levels[1], "Masked dataset index"),
  curve_band(mat_fluid_b, panel_levels[2], "Masked dataset index"),
  curve_band(mat_fluid_e, panel_levels[3], "Retrospective death-event index"),
  curve_band(mat_catheter_vfd, panel_levels[4], "Masked dataset index"),
  curve_band(mat_catheter_b, panel_levels[5], "Masked dataset index"),
  curve_band(mat_catheter_e, panel_levels[6], "Retrospective death-event index")
))
bands$panel <- factor(bands$panel, levels = panel_levels)

observed <- do.call(rbind, list(
  observed_curve(obs_fluid_vfd, panel_levels[1], "Masked dataset index"),
  observed_curve(obs_fluid_b, panel_levels[2], "Masked dataset index"),
  observed_curve(obs_fluid_e, panel_levels[3], "Retrospective death-event index"),
  observed_curve(obs_catheter_vfd, panel_levels[4], "Masked dataset index"),
  observed_curve(obs_catheter_b, panel_levels[5], "Masked dataset index"),
  observed_curve(obs_catheter_e, panel_levels[6], "Retrospective death-event index")
))
observed$panel <- factor(observed$panel, levels = panel_levels)

plot_file_png <- file.path(OUT_DIR, "factt_adaptive_permutation_figure.png")
plot_file_pdf <- file.path(OUT_DIR, "factt_adaptive_permutation_figure.pdf")
manuscript_png <- file.path(MANUSCRIPT_DIR, "figS2_factt_permutation.png")
manuscript_pdf <- file.path(MANUSCRIPT_DIR, "figS2_factt_permutation.pdf")

p <- ggplot(bands, aes(x = index)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#9ecae1", alpha = 0.45) +
  geom_line(aes(y = median), color = "#08519c", linewidth = 0.75) +
  geom_line(data = observed, aes(y = wealth), color = "#252525", linewidth = 0.55, alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey50", linewidth = 0.35) +
  geom_hline(yintercept = THRESHOLD, linetype = "dashed", color = "#b23a48", linewidth = 0.5) +
  facet_wrap(~panel, scales = "free_x", ncol = 2) +
  scale_y_log10() +
  labs(
    title = "FACTT adaptive e-RT order sensitivity",
    subtitle = "Blue line/band: median and central 95% band across 10,000 random patient-order permutations.\nBlack line: observed masked BioLINCC order.",
    x = "Retrospective dataset/event index",
    y = "e-value"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 10.2),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.margin = margin(10, 18, 10, 30)
  )

ggsave(plot_file_png, p, width = 10.5, height = 8.6, dpi = 360, bg = "white")
ggsave(plot_file_pdf, p, width = 10.5, height = 8.6, device = "pdf", bg = "white")
ggsave(manuscript_png, p, width = 10.5, height = 8.6, dpi = 360, bg = "white")
ggsave(manuscript_pdf, p, width = 10.5, height = 8.6, device = "pdf", bg = "white")

final_fluid_vfd_d <- d_estimate(factt$arm_fluid, factt$vfd28)
final_fluid_vfd_diff <- mean_diff_estimate(factt$arm_fluid, factt$vfd28)
final_catheter_vfd_d <- d_estimate(factt$arm_catheter, factt$vfd28)
final_catheter_vfd_diff <- mean_diff_estimate(factt$arm_catheter, factt$vfd28)
final_fluid_mort <- arr_estimate(factt$arm_fluid, factt$death_d60)
final_catheter_mort <- arr_estimate(factt$arm_catheter, factt$death_d60)

summary_raw <- do.call(rbind, list(
  summarize_permutation(
    "Fluid VFD-28 e-RTc", "patient", "d", obs_fluid_vfd,
    cross_fluid_vfd, effect_fluid_vfd, final_fluid_vfd_d,
    mean_diff_fluid_vfd, final_fluid_vfd_diff
  ),
  summarize_permutation(
    "Fluid mortality e-RTb", "patient", "ARR", obs_fluid_b,
    cross_fluid_b, effect_fluid_b, final_fluid_mort
  ),
  summarize_permutation(
    "Fluid mortality e-RTe", "event", "ARR", obs_fluid_e,
    cross_fluid_e, effect_fluid_e, final_fluid_mort
  ),
  summarize_permutation(
    "Catheter VFD-28 e-RTc", "patient", "d", obs_catheter_vfd,
    cross_catheter_vfd, effect_catheter_vfd, final_catheter_vfd_d,
    mean_diff_catheter_vfd, final_catheter_vfd_diff
  ),
  summarize_permutation(
    "Catheter mortality e-RTb", "patient", "ARR", obs_catheter_b,
    cross_catheter_b, effect_catheter_b, final_catheter_mort
  ),
  summarize_permutation(
    "Catheter mortality e-RTe", "event", "ARR", obs_catheter_e,
    cross_catheter_e, effect_catheter_e, final_catheter_mort
  )
))

effect_display <- mapply(
  fmt_effect,
  summary_raw$median_effect_at_crossing,
  summary_raw$median_mean_diff_at_crossing,
  summary_raw$effect_scale,
  SIMPLIFY = TRUE
)
final_effect_display <- mapply(
  fmt_effect,
  summary_raw$final_effect,
  summary_raw$final_mean_difference,
  summary_raw$effect_scale,
  SIMPLIFY = TRUE
)

summary_display <- data.frame(
  Curve = summary_raw$curve,
  `Observed masked order` = summary_raw$observed,
  `Permutation crossings` = fmt_pct(summary_raw$crossing_rate, 1),
  `Median crossing index` = fmt_int(summary_raw$median_crossing),
  `Central 95% crossing index` = ifelse(
    is.na(summary_raw$q025_crossing),
    "",
    paste0(fmt_int(summary_raw$q025_crossing), "-", fmt_int(summary_raw$q975_crossing))
  ),
  `Median effect at crossing` = effect_display,
  `Final effect` = final_effect_display,
  `Median Type M` = ifelse(is.na(summary_raw$median_type_m), "", paste0(fmt_num(summary_raw$median_type_m, 2), "x")),
  check.names = FALSE
)

summary_csv <- file.path(OUT_DIR, "factt_adaptive_permutation_summary.csv")
bands_csv <- file.path(OUT_DIR, "factt_adaptive_permutation_bands.csv")
table_md <- file.path(TABLE_DIR, "factt_adaptive_permutation_sensitivity.md")
table_tex <- file.path(TABLE_DIR, "factt_adaptive_permutation_sensitivity.tex")

write.csv(summary_raw, summary_csv, row.names = FALSE)
write.csv(bands, bands_csv, row.names = FALSE)
write_md_table(summary_display, table_md)
write_tex_table(summary_display, table_tex)

message("Saved: ", plot_file_png)
message("Saved: ", plot_file_pdf)
message("Saved: ", manuscript_png)
message("Saved: ", manuscript_pdf)
message("Saved: ", summary_csv)
message("Saved: ", bands_csv)
message("Saved: ", table_md)
message("Saved: ", table_tex)
