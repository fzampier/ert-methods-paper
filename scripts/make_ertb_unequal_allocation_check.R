# Freeze amendment 4 (R2.m1 support): e-RTb validity under 2:1 allocation.
#
# R2.m1 asks the e-RTb section to be generalized to unequal allocation
# (lambda_i^* = Pr_design(T_i = 1 | Y_i) with allocation p != 0.5). This
# script substantiates the generalized formula with numbers: under the null
# with 2:1 allocation (p = 2/3), both the generalized design wager and the
# shipped adaptive wager keep the crossing rate at or below 5% (any
# predictable lambda yields a martingale when the update divides by the true
# p and 1-p); under the alternative the design wager retains material power.
#
# Note: the shipped adaptive wager centers its bet at 0.5 regardless of p.
# That is VALID for any p (predictability is all the martingale needs) but
# mis-centered relative to the allocation, so its 2:1 power is not
# representative; the design row is the informative one. Stated in the table
# footnote and, in the manuscript, in the generalized-e-RTb passage.
#
# Output: tables/ertb_unequal_allocation_check.{csv,md}

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

set.seed(20260549)
n_sims <- as.integer(Sys.getenv("UNEQ_NSIMS", "5000"))
n_patients <- 2690
p_alloc <- 2 / 3
p_ctrl <- 0.35
p_trt_alt <- 0.30
threshold <- 20

simulate_trial <- function(n, p_trt, p_ctrl, p_random) {
  treatment <- rbinom(n, 1, p_random)
  outcome <- numeric(n)
  outcome[treatment == 1] <- rbinom(sum(treatment == 1), 1, p_trt)
  outcome[treatment == 0] <- rbinom(sum(treatment == 0), 1, p_ctrl)
  list(treatment = treatment, outcome = outcome)
}

run_cell <- function(p_trt, scenario) {
  crossed_design <- logical(n_sims)
  crossed_adaptive <- logical(n_sims)
  for (sim in seq_len(n_sims)) {
    trial <- simulate_trial(n_patients, p_trt, p_ctrl, p_alloc)
    w_design <- compute_eRT(trial$treatment, trial$outcome, p = p_alloc,
                            wager = "design",
                            p_trt_design = p_trt_alt, p_ctrl_design = p_ctrl,
                            ramp_fixed = FALSE)
    w_adaptive <- compute_eRT(trial$treatment, trial$outcome, p = p_alloc,
                              burn_in = 50, ramp = 100,
                              wager = "adaptive", kelly_fraction = 0.5)
    crossed_design[sim] <- any(w_design >= threshold)
    crossed_adaptive[sim] <- any(w_adaptive >= threshold)
  }
  rbind(
    data.frame(scenario = scenario, allocation = "2:1 (p=2/3)",
               method = "e-RTb design (generalized lambda*)",
               crossing_rate = mean(crossed_design),
               se = sqrt(mean(crossed_design) * (1 - mean(crossed_design)) / n_sims)),
    data.frame(scenario = scenario, allocation = "2:1 (p=2/3)",
               method = "e-RTb adaptive (as shipped; 0.5-centered bet)",
               crossing_rate = mean(crossed_adaptive),
               se = sqrt(mean(crossed_adaptive) * (1 - mean(crossed_adaptive)) / n_sims))
  )
}

t0 <- Sys.time()
out <- rbind(
  run_cell(p_ctrl, "Null: 35% vs 35%"),
  run_cell(p_trt_alt, "Alternative: 35% vs 30%")
)
write.csv(out, file.path(table_dir, "ertb_unequal_allocation_check.csv"), row.names = FALSE)

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
md <- data.frame(
  Scenario = out$scenario, Allocation = out$allocation, Method = out$method,
  `Crossing rate (W>=20)` = sprintf("%s (%.2fpp)", fmt_pct(out$crossing_rate), 100 * out$se),
  check.names = FALSE
)
writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  sprintf("N=%s, %s simulated trials per cell, seed 20260549. Wealth updates divide by the true allocation probabilities (p=2/3, 1-p=1/3). Any predictable wager is valid at any allocation; the shipped adaptive wager's bet is centered at 0.5 rather than at p, so its power under 2:1 is not representative -- the generalized design row is the informative one (R2.m1).",
          format(n_patients, big.mark = ","), format(n_sims, big.mark = ","))
), file.path(table_dir, "ertb_unequal_allocation_check.md"))

cat(sprintf("Done in %.1f min.\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(md, row.names = FALSE)
