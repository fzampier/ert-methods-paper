# Binary interruption-curve simulation of the original submission (textbook
# boundaries; ESM Table S19). This script is standalone inside this repository and uses only the
# local R/ertb.R implementation.

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

root <- paper_root()
manuscript_dir <- paper_file("figures", root = root, must_work = FALSE)
table_dir <- paper_file("tables", root = root, must_work = FALSE)

dir.create(manuscript_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260511)

alpha <- 0.05
threshold <- 1 / alpha
n_sims <- 5000
n_patients <- 2690
p_ctrl <- 0.35
p_trt_alt <- 0.30
p_trt_null <- 0.35
look_fracs <- c(1 / 3, 2 / 3, 1)
look_n <- round(n_patients * look_fracs)

obf_bounds <- qnorm(1 - alpha / 2) / sqrt(look_fracs)
hp_interim_alpha <- 0.001
hp_bounds <- c(rep(qnorm(1 - hp_interim_alpha / 2), 2), qnorm(1 - alpha / 2))
# Pocock K=3, two-sided 5%, equally-spaced looks: constant |Z| >= 2.289
# (Pocock 1977; Jennison & Turnbull 2000 Table 2.1).
pocock_bounds <- rep(2.289, 3)

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

  # Positive Z means benefit for treatment: fewer events in treatment than control.
  (p_c - p_t) / se
}

arr_at <- function(treatment, outcome, n_look) {
  idx <- seq_len(n_look)
  trt <- treatment[idx] == 1
  ctrl <- treatment[idx] == 0
  if (!any(trt) || !any(ctrl)) return(NA_real_)
  mean(outcome[idx][ctrl]) - mean(outcome[idx][trt])
}

freq_first_crossing <- function(treatment, outcome, bounds) {
  for (k in seq_along(look_n)) {
    z <- z_at_look(treatment, outcome, look_n[[k]])
    if (is.finite(z) && abs(z) >= bounds[[k]]) {
      return(list(
        first_crossing = look_n[[k]],
        direction = ifelse(z > 0, "benefit", "harm"),
        statistic = z
      ))
    }
  }
  list(first_crossing = NA_integer_, direction = NA_character_, statistic = NA_real_)
}

ert_first_crossing <- function(wealth, treatment, outcome) {
  cross <- which(wealth >= threshold)
  if (!length(cross)) {
    return(list(first_crossing = NA_integer_, direction = NA_character_, statistic = NA_real_))
  }
  fc <- cross[[1]]
  arr <- arr_at(treatment, outcome, fc)
  list(
    first_crossing = fc,
    direction = ifelse(is.finite(arr) && arr >= 0, "benefit", "harm"),
    statistic = arr
  )
}

run_methods <- function(p_trt, scenario, seed) {
  set.seed(seed)
  rows <- vector("list", n_sims * 5L)
  row_i <- 1L

  for (sim in seq_len(n_sims)) {
    trial <- simulate_trial(n_patients, p_trt = p_trt, p_ctrl = p_ctrl)

    obf <- freq_first_crossing(trial$treatment, trial$outcome, obf_bounds)
    hp <- freq_first_crossing(trial$treatment, trial$outcome, hp_bounds)
    pocock <- freq_first_crossing(trial$treatment, trial$outcome, pocock_bounds)

    wealth_adaptive <- compute_eRT(
      trial$treatment,
      trial$outcome,
      burn_in = 50,
      ramp = 100,
      wager = "adaptive",
      kelly_fraction = 0.5
    )
    adaptive <- ert_first_crossing(wealth_adaptive, trial$treatment, trial$outcome)

    wealth_design <- compute_eRT(
      trial$treatment,
      trial$outcome,
      wager = "design",
      p_trt_design = p_trt_alt,
      p_ctrl_design = p_ctrl,
      ramp_fixed = FALSE
    )
    design <- ert_first_crossing(wealth_design, trial$treatment, trial$outcome)

    add_row <- function(method, res) {
      effect_at_crossing <- if (!is.na(res$first_crossing)) {
        arr_at(trial$treatment, trial$outcome, res$first_crossing)
      } else {
        NA_real_
      }
      final_effect <- arr_at(trial$treatment, trial$outcome, n_patients)
      data.frame(
        sim = sim,
        scenario = scenario,
        method = method,
        first_crossing = res$first_crossing,
        direction = res$direction,
        statistic_at_crossing = res$statistic,
        effect_at_crossing = effect_at_crossing,
        final_effect = final_effect,
        type_m = abs(effect_at_crossing) / abs(final_effect)
      )
    }

    rows[[row_i]] <- add_row("OBF two-sided", obf)
    row_i <- row_i + 1L
    rows[[row_i]] <- add_row("Pocock two-sided", pocock)
    row_i <- row_i + 1L
    rows[[row_i]] <- add_row("Haybittle-Peto", hp)
    row_i <- row_i + 1L
    rows[[row_i]] <- add_row("e-RTb adaptive", adaptive)
    row_i <- row_i + 1L
    rows[[row_i]] <- add_row("e-RTb design directional", design)
    row_i <- row_i + 1L
  }

  bind_rows(rows)
}

alternative <- run_methods(p_trt_alt, "Alternative: 35% vs 30%", seed = 20260542)
null <- run_methods(p_trt_null, "Null: 35% vs 35%", seed = 20260543)
crossings <- bind_rows(alternative, null)

method_order <- c(
  "OBF two-sided",
  "Pocock two-sided",
  "Haybittle-Peto",
  "e-RTb adaptive",
  "e-RTb design directional"
)

method_summary <- crossings %>%
  mutate(method = factor(method, levels = method_order)) %>%
  group_by(scenario, method) %>%
  summarise(
    interrupted_any = mean(!is.na(first_crossing)),
    interrupted_benefit = mean(direction == "benefit", na.rm = TRUE),
    interrupted_harm = mean(direction == "harm", na.rm = TRUE),
    se_any = sqrt(interrupted_any * (1 - interrupted_any) / n()),
    median_interruption = median(first_crossing[!is.na(first_crossing)], na.rm = TRUE),
    q25_interruption = quantile(first_crossing[!is.na(first_crossing)], 0.25, na.rm = TRUE),
    q75_interruption = quantile(first_crossing[!is.na(first_crossing)], 0.75, na.rm = TRUE),
    median_type_m = median(type_m[!is.na(first_crossing)], na.rm = TRUE),
    .groups = "drop"
  )

joint_one <- function(wide, base, ert, comparison) {
  base_cross <- !is.na(wide[[base]])
  ert_cross <- !is.na(wide[[ert]])
  base_time <- wide[[base]]
  ert_time <- wide[[ert]]
  both <- base_cross & ert_cross

  data.frame(
    comparison = comparison,
    combined = mean(base_cross | ert_cross),
    base_only = mean(base_cross & !ert_cross),
    ert_only = mean(!base_cross & ert_cross),
    both = mean(both),
    neither = mean(!base_cross & !ert_cross),
    ert_before_base_given_both = ifelse(any(both), mean(ert_time[both] < base_time[both]), NA_real_),
    base_before_ert_given_both = ifelse(any(both), mean(base_time[both] < ert_time[both]), NA_real_),
    same_time_given_both = ifelse(any(both), mean(base_time[both] == ert_time[both]), NA_real_)
  )
}

joint_summary <- function(rows, scenario) {
  wide <- rows %>%
    select(sim, method, first_crossing) %>%
    pivot_wider(names_from = method, values_from = first_crossing)

  bind_rows(
    joint_one(wide, "OBF two-sided", "e-RTb adaptive", "OBF + adaptive e-RTb"),
    joint_one(wide, "OBF two-sided", "e-RTb design directional", "OBF + design e-RTb"),
    joint_one(wide, "Pocock two-sided", "e-RTb adaptive", "Pocock + adaptive e-RTb"),
    joint_one(wide, "Pocock two-sided", "e-RTb design directional", "Pocock + design e-RTb"),
    joint_one(wide, "Haybittle-Peto", "e-RTb adaptive", "HP + adaptive e-RTb"),
    joint_one(wide, "Haybittle-Peto", "e-RTb design directional", "HP + design e-RTb")
  ) %>%
    mutate(scenario = scenario, .before = 1)
}

joint <- bind_rows(
  joint_summary(alternative, "Alternative: 35% vs 30%"),
  joint_summary(null, "Null: 35% vs 35%")
)

curve_grid <- crossings %>%
  mutate(method = factor(method, levels = method_order)) %>%
  group_by(scenario, method) %>%
  reframe(
    patient = seq_len(n_patients),
    cumulative_interruption = sapply(
      patient,
      function(j) mean(!is.na(first_crossing) & first_crossing <= j)
    )
  ) %>%
  ungroup()

method_cols <- c(
  "OBF two-sided" = "#2563eb",
  "Pocock two-sided" = "#dc2626",
  "Haybittle-Peto" = "#7c3aed",
  "e-RTb adaptive" = "#f97316",
  "e-RTb design directional" = "#16a34a"
)

plot <- ggplot(curve_grid, aes(patient, cumulative_interruption, color = method)) +
  geom_step(linewidth = 0.9) +
  geom_vline(
    xintercept = look_n,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey55"
  ) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dotted",
    linewidth = 0.35,
    color = "grey35"
  ) +
  facet_wrap(~scenario, ncol = 1, scales = "free_y") +
  scale_color_manual(values = method_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Interruption curves for binary monitoring",
    subtitle = sprintf(
      "35%% vs 30%% event rates; %s trials per scenario; OBF calibrated at N = %s",
      format(n_sims, big.mark = ","),
      format(n_patients, big.mark = ",")
    ),
    x = "Patients enrolled",
    y = "Cumulative probability of interruption",
    color = NULL,
    caption = sprintf(
      "Looks: %s. OBF |Z|: %s. Pocock |Z|: %s. Haybittle-Peto |Z|: %s. e-RT W=%.0f.",
      paste(format(look_n, big.mark = ","), collapse = ", "),
      paste(sprintf("%.2f", obf_bounds), collapse = ", "),
      paste(sprintf("%.2f", pocock_bounds), collapse = ", "),
      paste(sprintf("%.2f", hp_bounds), collapse = ", "),
      threshold
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", hjust = 0),
    plot.caption = element_text(hjust = 0, size = 9)
  )

ggsave(file.path(manuscript_dir, "legacy_fig2_interruption_curves.png"), plot, width = 9, height = 7.8, dpi = 180)
ggsave(file.path(manuscript_dir, "legacy_fig2_interruption_curves.pdf"), plot, width = 9, height = 7.8)

write.csv(
  crossings,
  file.path(table_dir, "ertb_interruption_first_crossings.csv"),
  row.names = FALSE
)
write.csv(
  method_summary,
  file.path(table_dir, "ertb_interruption_method_summary.csv"),
  row.names = FALSE
)
write.csv(
  joint,
  file.path(table_dir, "ertb_interruption_joint_summary.csv"),
  row.names = FALSE
)

alt_summary <- method_summary %>% filter(scenario == "Alternative: 35% vs 30%")
null_summary <- method_summary %>% filter(scenario == "Null: 35% vs 35%")
joint_alt <- joint %>% filter(scenario == "Alternative: 35% vs 30%")

method_table <- alt_summary %>%
  select(method, interrupted_any, median_interruption, median_type_m) %>%
  left_join(
    null_summary %>% select(method, null_interruption = interrupted_any),
    by = "method"
  ) %>%
  left_join(
    joint_alt %>%
      filter(comparison %in% c("OBF + adaptive e-RTb", "OBF + design e-RTb")) %>%
      transmute(
        method = ifelse(comparison == "OBF + adaptive e-RTb", "e-RTb adaptive", "e-RTb design directional"),
        earlier_than_obf = ert_before_base_given_both
      ),
    by = "method"
  ) %>%
  left_join(
    joint_alt %>%
      filter(comparison %in% c("Pocock + adaptive e-RTb", "Pocock + design e-RTb")) %>%
      transmute(
        method = ifelse(comparison == "Pocock + adaptive e-RTb", "e-RTb adaptive", "e-RTb design directional"),
        earlier_than_pocock = ert_before_base_given_both
      ),
    by = "method"
  ) %>%
  left_join(
    joint_alt %>%
      filter(comparison %in% c("HP + adaptive e-RTb", "HP + design e-RTb")) %>%
      transmute(
        method = ifelse(comparison == "HP + adaptive e-RTb", "e-RTb adaptive", "e-RTb design directional"),
        earlier_than_hp = ert_before_base_given_both
      ),
    by = "method"
  ) %>%
  mutate(
    Method = as.character(method),
    `Alternative interruption` = fmt_pct(interrupted_any),
    `Median interruption` = fmt_num(median_interruption),
    `Null interruption` = fmt_pct(null_interruption),
    `Median Type M` = ifelse(is.na(median_type_m), "--", sprintf("%.2fx", median_type_m)),
    `Earlier than OBF when both interrupt` = fmt_pct(earlier_than_obf),
    `Earlier than Pocock when both interrupt` = fmt_pct(earlier_than_pocock),
    `Earlier than HP when both interrupt` = fmt_pct(earlier_than_hp)
  ) %>%
  select(
    Method,
    `Alternative interruption`,
    `Median interruption`,
    `Null interruption`,
    `Median Type M`,
    `Earlier than OBF when both interrupt`,
    `Earlier than Pocock when both interrupt`,
    `Earlier than HP when both interrupt`
  )

write_md_table(
  method_table,
  file.path(table_dir, "ertb_interruption_operating_characteristics.md")
)

joint_table <- joint_alt %>%
  filter(comparison %in% c("OBF + design e-RTb", "HP + design e-RTb")) %>%
  mutate(
    Combination = comparison,
    `Either trigger` = fmt_pct(combined),
    `Primary only` = fmt_pct(base_only),
    `e-RT only` = fmt_pct(ert_only),
    `Both` = fmt_pct(both),
    `e-RT earlier when both interrupt` = fmt_pct(ert_before_base_given_both)
  ) %>%
  select(Combination, `Either trigger`, `Primary only`, `e-RT only`, Both, `e-RT earlier when both interrupt`)

write_md_table(
  joint_table,
  file.path(table_dir, "ertb_interruption_joint_design_summary.md")
)

cat("Wrote interruption figure and tables.\n")
cat("Figure: ", file.path(manuscript_dir, "legacy_fig2_interruption_curves.png"), "\n", sep = "")
cat("Table:  ", file.path(table_dir, "ertb_interruption_operating_characteristics.md"), "\n", sep = "")
