# Design-assumption versus actual-effect power table for e-RTb.

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260514)

p_ctrls <- c(0.40, 0.25)
alpha <- 0.05
threshold <- 1 / alpha
n_sims <- 5000
design_arrs <- c(0.05, 0.10)
actual_arrs <- c(0.05, 0.10)
target_powers <- c(0.80, 0.90)

fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
fmt_rate_pair <- function(p_ctrl, arr) sprintf("%.0f%% vs %.0f%%", 100 * p_ctrl, 100 * (p_ctrl - arr))

write_md_table <- function(x, path) {
  x[] <- lapply(x, as.character)
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |")
  rows <- apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), path)
  invisible(path)
}

sample_size_total <- function(p_ctrl, arr, target_power) {
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
  # Delegates to compute_eRT(wager="design", ramp_fixed=FALSE) so the design
  # wealth path in this table tracks any future change to the canonical kernel.
  # ramp_fixed = FALSE applies the design-implied lambda from patient 2 onward
  # without the c_i ramp-up; this matches the previous local implementation.
  compute_eRT(
    treatment = treatment,
    outcome = outcome,
    p = 0.5,
    wager = "design",
    p_trt_design = p_trt_design,
    p_ctrl_design = p_ctrl_design,
    ramp_fixed = FALSE
  )
}

crossed_threshold <- function(wealth) any(wealth >= threshold)

rows <- list()
row_i <- 1L

for (p_ctrl in p_ctrls) {
  for (design_arr in design_arrs) {
    for (target_power in target_powers) {
      n <- sample_size_total(p_ctrl, design_arr, target_power)
      actual_freq_power <- function(actual_arr) {
        power.prop.test(
          n = n / 2,
          p1 = p_ctrl,
          p2 = p_ctrl - actual_arr,
          sig.level = alpha
        )$power
      }

      for (actual_arr in actual_arrs) {
        adaptive_crossed <- logical(n_sims)
        design_crossed <- logical(n_sims)

        for (sim in seq_len(n_sims)) {
          trial <- simulate_trial_fast(n, p_ctrl - actual_arr, p_ctrl)

          adaptive_wealth <- compute_eRT(
            trial$treatment,
            trial$outcome,
            burn_in = 50,
            ramp = 100,
            wager = "adaptive",
            kelly_fraction = 0.5
          )

          design_wealth <- compute_design_wealth(
            trial$treatment,
            trial$outcome,
            p_trt_design = p_ctrl - design_arr,
            p_ctrl_design = p_ctrl
          )

          adaptive_crossed[[sim]] <- crossed_threshold(adaptive_wealth)
          design_crossed[[sim]] <- crossed_threshold(design_wealth)
        }

        rows[[row_i]] <- data.frame(
          p_ctrl = p_ctrl,
          design_arr = design_arr,
          actual_arr = actual_arr,
          target_power = target_power,
          n_patients = n,
          actual_frequentist_power = actual_freq_power(actual_arr),
          ertb_adaptive_power = mean(adaptive_crossed),
          ertb_design_power = mean(design_crossed)
        )
        row_i <- row_i + 1L
      }
    }
  }
}

# ---------------------------------------------------------------------------
# T4 reversal block (BMC revision, R2.4/R4.m8): actual effect in the OPPOSITE
# direction (treatment worse than control), while the design wager stays
# calibrated to the anticipated benefit. Encoded as negative actual_arr, so
# actual event rates are p_ctrl (control) vs p_ctrl + |arr| (treatment).
# APPENDED AFTER the original grid under a FRESH seed -- the original 16 cells
# above run first on the untouched RNG stream and must reproduce exactly.
# ---------------------------------------------------------------------------

set.seed(20260546)
reversed_arrs <- c(-0.05, -0.10)

for (p_ctrl in p_ctrls) {
  for (design_arr in design_arrs) {
    for (target_power in target_powers) {
      n <- sample_size_total(p_ctrl, design_arr, target_power)
      actual_freq_power <- function(actual_arr) {
        power.prop.test(
          n = n / 2,
          p1 = p_ctrl,
          p2 = p_ctrl - actual_arr,
          sig.level = alpha
        )$power
      }

      for (actual_arr in reversed_arrs) {
        adaptive_crossed <- logical(n_sims)
        design_crossed <- logical(n_sims)

        for (sim in seq_len(n_sims)) {
          trial <- simulate_trial_fast(n, p_ctrl - actual_arr, p_ctrl)

          adaptive_wealth <- compute_eRT(
            trial$treatment,
            trial$outcome,
            burn_in = 50,
            ramp = 100,
            wager = "adaptive",
            kelly_fraction = 0.5
          )

          design_wealth <- compute_design_wealth(
            trial$treatment,
            trial$outcome,
            p_trt_design = p_ctrl - design_arr,
            p_ctrl_design = p_ctrl
          )

          adaptive_crossed[[sim]] <- crossed_threshold(adaptive_wealth)
          design_crossed[[sim]] <- crossed_threshold(design_wealth)
        }

        rows[[row_i]] <- data.frame(
          p_ctrl = p_ctrl,
          design_arr = design_arr,
          actual_arr = actual_arr,
          target_power = target_power,
          n_patients = n,
          actual_frequentist_power = actual_freq_power(actual_arr),
          ertb_adaptive_power = mean(adaptive_crossed),
          ertb_design_power = mean(design_crossed)
        )
        row_i <- row_i + 1L
      }
    }
  }
}

out <- do.call(rbind, rows)
out$direction <- ifelse(out$actual_arr >= 0, "as-designed", "reversed")
write.csv(
  out,
  file.path(table_dir, "ertb_design_actual_power.csv"),
  row.names = FALSE
)

md <- data.frame(
  Design = fmt_rate_pair(out$p_ctrl, out$design_arr),
  `Actual effect` = fmt_rate_pair(out$p_ctrl, out$actual_arr),
  `Design power` = fmt_pct(out$target_power, digits = 0),
  `e-RTb adaptive power` = fmt_pct(out$ertb_adaptive_power),
  `e-RTb design power` = fmt_pct(out$ertb_design_power),
  `Actual frequentist power` = fmt_pct(out$actual_frequentist_power),
  check.names = FALSE
)

write_md_table(md, file.path(table_dir, "ertb_design_actual_power.md"))

cat("Wrote design-versus-actual e-RTb power table.\n")
cat("Table: ", file.path(table_dir, "ertb_design_actual_power.md"), "\n", sep = "")
