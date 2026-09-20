# ARMA adaptive-only permutation sensitivity figure and table.
#
# This script uses the ignored local validated BioLINCC ARMA CSV and writes only
# aggregate, manuscript-safe outputs. No shuffled patient-level datasets are
# saved.

suppressPackageStartupMessages({
  library(ggplot2)
})

source("scripts/paths.R")

suppressPackageStartupMessages(source(local_ert_file("ertb.R")))
suppressPackageStartupMessages(source(local_ert_file("erte.R")))

OUT_DIR <- file.path("outputs", "arma_permutation")
MANUSCRIPT_DIR <- paper_file("figures")
TABLE_DIR <- paper_file("tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 10000L
SEED <- 20260506L
THRESHOLD <- 20
EPS <- 1e-8

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

summarize_permutation <- function(label, update_unit, observed_wealth,
                                  crossings, crossing_effects, final_effect) {
  crossed <- !is.na(crossings)
  effect_cross <- crossing_effects[crossed]
  type_m <- abs(effect_cross) / abs(final_effect)
  observed_crossing_label <- switch(
    update_unit,
    patient = "masked index",
    `death event` = "death-event index",
    `failure event` = "failure-event index",
    "index"
  )

  data.frame(
    curve = label,
    update_unit = update_unit,
    observed = if (any(observed_wealth >= THRESHOLD)) {
      paste0(
        "Crossed at ", observed_crossing_label, " ",
        fmt_int(first_crossing(observed_wealth))
      )
    } else {
      paste0("No crossing; peak E = ", fmt_num(max(observed_wealth), 2))
    },
    crossing_rate = mean(crossed),
    median_crossing = if (any(crossed)) median(crossings[crossed]) else NA_real_,
    q025_crossing = if (any(crossed)) unname(stats::quantile(crossings[crossed], 0.025)) else NA_real_,
    q975_crossing = if (any(crossed)) unname(stats::quantile(crossings[crossed], 0.975)) else NA_real_,
    median_effect_at_crossing = if (length(effect_cross)) median(effect_cross) else NA_real_,
    final_effect = final_effect,
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
    "\\caption{ARMA adaptive e-RT order-sensitivity summary from 10,000 complete random patient-order permutations. Crossing index summaries are among permutations that crossed. Type M is the absolute ratio of the effect at first crossing to the final effect.}",
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste(names(x), collapse = " & "),
    "\\\\ \\hline"
  )
  body <- apply(x, 1, function(row) paste0(paste(row, collapse = " & "), "\\\\"))
  lines <- c(lines, body, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
}

arma <- read.csv(validated_trial_csv("arma"), stringsAsFactors = FALSE)
arma <- arma[order(arma$enroll_order), ]
arma$arm <- as.integer(arma$arm)
arma$death <- as.integer(arma$death)
arma$bwa28 <- as.integer(!is.na(arma$unassist_day) & arma$unassist_day <= 28)
arma$bwa28_fail <- 1L - arma$bwa28

n <- nrow(arma)
death_events <- sum(arma$death == 1L)
bwa_events <- sum(arma$bwa28_fail == 1L)

panel_levels <- c(
  "A. Mortality e-RTb",
  "B. Mortality e-RTe",
  "C. BWA-28 failure e-RTb",
  "D. BWA-28 failure e-RTe"
)

obs_mort_b <- compute_eRT(
  arma$arm, arma$death,
  burn_in = 50, ramp = 100, wager = "adaptive"
)
obs_mort_e <- compute_eRTe(
  arma$arm[arma$death == 1L],
  burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
)$wealth
obs_bwa_b <- compute_eRT(
  arma$arm, arma$bwa28_fail,
  burn_in = 50, ramp = 100, wager = "adaptive"
)
obs_bwa_e <- compute_eRTe(
  arma$arm[arma$bwa28_fail == 1L],
  burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
)$wealth

mat_mort_b <- matrix(NA_real_, nrow = N_PERM, ncol = n)
mat_mort_e <- matrix(NA_real_, nrow = N_PERM, ncol = death_events)
mat_bwa_b <- matrix(NA_real_, nrow = N_PERM, ncol = n)
mat_bwa_e <- matrix(NA_real_, nrow = N_PERM, ncol = bwa_events)

cross_mort_b <- rep(NA_integer_, N_PERM)
cross_mort_e <- rep(NA_integer_, N_PERM)
cross_bwa_b <- rep(NA_integer_, N_PERM)
cross_bwa_e <- rep(NA_integer_, N_PERM)

effect_mort_b <- rep(NA_real_, N_PERM)
effect_mort_e <- rep(NA_real_, N_PERM)
effect_bwa_b <- rep(NA_real_, N_PERM)
effect_bwa_e <- rep(NA_real_, N_PERM)

set.seed(SEED)
message("Running ", N_PERM, " complete random ARMA patient-order permutations...")
pb <- txtProgressBar(min = 0, max = N_PERM, style = 3)

for (sim in seq_len(N_PERM)) {
  ord <- sample.int(n)
  arm <- arma$arm[ord]
  death <- arma$death[ord]
  bwa_fail <- arma$bwa28_fail[ord]

  w <- compute_eRT(arm, death, burn_in = 50, ramp = 100, wager = "adaptive")
  mat_mort_b[sim, ] <- w
  cross_mort_b[sim] <- first_crossing(w)
  if (!is.na(cross_mort_b[sim])) {
    effect_mort_b[sim] <- arr_estimate(arm, death, cross_mort_b[sim])
  }

  event_pos <- which(death == 1L)
  r <- compute_eRTe(
    arm[event_pos],
    burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
  )$wealth
  mat_mort_e[sim, ] <- r
  cross_mort_e[sim] <- first_crossing(r)
  if (!is.na(cross_mort_e[sim])) {
    effect_mort_e[sim] <- arr_estimate(arm, death, event_pos[cross_mort_e[sim]])
  }

  w <- compute_eRT(arm, bwa_fail, burn_in = 50, ramp = 100, wager = "adaptive")
  mat_bwa_b[sim, ] <- w
  cross_bwa_b[sim] <- first_crossing(w)
  if (!is.na(cross_bwa_b[sim])) {
    effect_bwa_b[sim] <- arr_estimate(arm, bwa_fail, cross_bwa_b[sim])
  }

  event_pos <- which(bwa_fail == 1L)
  r <- compute_eRTe(
    arm[event_pos],
    burn_in = 30, ramp = 50, threshold = THRESHOLD, wager = "adaptive"
  )$wealth
  mat_bwa_e[sim, ] <- r
  cross_bwa_e[sim] <- first_crossing(r)
  if (!is.na(cross_bwa_e[sim])) {
    effect_bwa_e[sim] <- arr_estimate(arm, bwa_fail, event_pos[cross_bwa_e[sim]])
  }

  if (sim %% 100L == 0L) setTxtProgressBar(pb, sim)
}
setTxtProgressBar(pb, N_PERM)
close(pb)

bands <- do.call(rbind, list(
  curve_band(mat_mort_b, panel_levels[1], "Masked dataset index"),
  curve_band(mat_mort_e, panel_levels[2], "Retrospective death-event index"),
  curve_band(mat_bwa_b, panel_levels[3], "Masked dataset index"),
  curve_band(mat_bwa_e, panel_levels[4], "Retrospective failure-event index")
))
bands$panel <- factor(bands$panel, levels = panel_levels)

observed <- do.call(rbind, list(
  observed_curve(obs_mort_b, panel_levels[1], "Masked dataset index"),
  observed_curve(obs_mort_e, panel_levels[2], "Retrospective death-event index"),
  observed_curve(obs_bwa_b, panel_levels[3], "Masked dataset index"),
  observed_curve(obs_bwa_e, panel_levels[4], "Retrospective failure-event index")
))
observed$panel <- factor(observed$panel, levels = panel_levels)

plot_file_png <- file.path(OUT_DIR, "arma_adaptive_permutation_figure.png")
plot_file_pdf <- file.path(OUT_DIR, "arma_adaptive_permutation_figure.pdf")
manuscript_png <- file.path(MANUSCRIPT_DIR, "figS1_arma_permutation.png")
manuscript_pdf <- file.path(MANUSCRIPT_DIR, "figS1_arma_permutation.pdf")

p <- ggplot(bands, aes(x = index)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#9ecae1", alpha = 0.45) +
  geom_line(aes(y = median), color = "#08519c", linewidth = 0.75) +
  geom_line(data = observed, aes(y = wealth), color = "#252525", linewidth = 0.55, alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey50", linewidth = 0.35) +
  geom_hline(yintercept = THRESHOLD, linetype = "dashed", color = "#b23a48", linewidth = 0.5) +
  facet_wrap(~panel, scales = "free_x", ncol = 2) +
  scale_y_log10() +
  labs(
    title = "ARMA adaptive e-RT order sensitivity",
    subtitle = "Blue line/band: median and central 95% band across 10,000 random patient-order permutations.\nBlack line: observed masked BioLINCC order.",
    x = "Retrospective dataset/event index",
    y = "e-value"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 10.5),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.margin = margin(10, 18, 10, 30)
  )

ggsave(plot_file_png, p, width = 10.5, height = 6.2, dpi = 360, bg = "white")
ggsave(plot_file_pdf, p, width = 10.5, height = 6.2, device = "pdf", bg = "white")
ggsave(manuscript_png, p, width = 10.5, height = 6.2, dpi = 360, bg = "white")
ggsave(manuscript_pdf, p, width = 10.5, height = 6.2, device = "pdf", bg = "white")

final_mort <- arr_estimate(arma$arm, arma$death)
final_bwa <- arr_estimate(arma$arm, arma$bwa28_fail)

summary_raw <- do.call(rbind, list(
  summarize_permutation(
    "Mortality e-RTb", "patient", obs_mort_b,
    cross_mort_b, effect_mort_b, final_mort
  ),
  summarize_permutation(
    "Mortality e-RTe", "death event", obs_mort_e,
    cross_mort_e, effect_mort_e, final_mort
  ),
  summarize_permutation(
    "BWA-28 failure e-RTb", "patient", obs_bwa_b,
    cross_bwa_b, effect_bwa_b, final_bwa
  ),
  summarize_permutation(
    "BWA-28 failure e-RTe", "failure event", obs_bwa_e,
    cross_bwa_e, effect_bwa_e, final_bwa
  )
))

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
  `Median effect at crossing` = fmt_pp(summary_raw$median_effect_at_crossing, 1),
  `Final effect` = fmt_pp(summary_raw$final_effect, 1),
  `Median Type M` = ifelse(is.na(summary_raw$median_type_m), "", paste0(fmt_num(summary_raw$median_type_m, 2), "x")),
  check.names = FALSE
)

summary_csv <- file.path(OUT_DIR, "arma_adaptive_permutation_summary.csv")
bands_csv <- file.path(OUT_DIR, "arma_adaptive_permutation_bands.csv")
table_md <- file.path(TABLE_DIR, "arma_adaptive_permutation_sensitivity.md")
table_tex <- file.path(TABLE_DIR, "arma_adaptive_permutation_sensitivity.tex")

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
