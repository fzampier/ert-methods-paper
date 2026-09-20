# T7: Unconditional cumulative crossing probabilities at the scheduled looks
# (all trials as denominator) for every method in the fair-comparator run.
# Answers R4.5's request for time-anchored unconditional metrics; post-processing
# of the per-trial first crossings -- no rerun, no RNG.
#
# Input : tables/fair_comparators_first_crossings.csv
# Output: tables/fair_comparators_cumulative_at_looks.{csv,md}

source("scripts/paths.R")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

cross <- read.csv(file.path(table_dir, "fair_comparators_first_crossings.csv"))
looks <- c(897, 1793, 2690)   # K=3 schedule: 1/3, 2/3, full information

cum <- cross %>%
  group_by(scenario, method) %>%
  summarise(
    n = n(),
    !!!setNames(
      lapply(looks, function(L) quo(mean(!is.na(first_crossing) & first_crossing <= !!L))),
      paste0("cum_", looks)
    ),
    .groups = "drop"
  )

write.csv(cum, file.path(table_dir, "fair_comparators_cumulative_at_looks.csv"), row.names = FALSE)

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
method_order <- c(
  "OBF two-sided", "Pocock two-sided", "Haybittle-Peto",
  "LD-OBF K3 two-sided a=0.05", "LD-OBF K5 two-sided a=0.05",
  "LD-OBF K3 one-sided a=0.05", "LD-OBF K5 one-sided a=0.05",
  "LD-OBF K3 one-sided a=0.025", "LD-OBF K5 one-sided a=0.025",
  "KD rho=1 K5 one-sided a=0.05", "KD rho=3 K5 one-sided a=0.05",
  "e-RTb adaptive", "e-RTb design directional", "e-RTb mixture two-sided"
)

md <- cum %>%
  mutate(method = factor(method, levels = method_order)) %>%
  arrange(scenario, method) %>%
  transmute(
    Scenario = scenario,
    Method = as.character(method),
    `P(crossed by 897)` = fmt_pct(cum_897),
    `P(crossed by 1793)` = fmt_pct(cum_1793),
    `P(crossed by 2690)` = fmt_pct(cum_2690)
  )

writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |")),
  "",
  "Cumulative probability of first crossing by the stated patient count, all simulated trials as denominator (n=5,000 per scenario). Looks at patients 897/1793/2690 correspond to the K=3 schedule at 1/3, 2/3, and full information."
), file.path(table_dir, "fair_comparators_cumulative_at_looks.md"))

cat("Alternative scenario, cumulative by first scheduled look (897):\n")
print(md %>% filter(Scenario == "Alternative: 35% vs 30%") %>% select(-Scenario), row.names = FALSE)
cat("\nWrote tables/fair_comparators_cumulative_at_looks.{csv,md}\n")
