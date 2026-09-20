# Worked-alert replay for the practical-deployment section.
# DETERMINISTIC: selects one representative trial from the frozen
# per-trial output of the union-doctrine simulation
# (tables/union_doctrine_first_crossings.csv) and replays
# its stream under the frozen seed to recover the full wealth and z paths.
# No new randomness; the replayed trial must reproduce its CSV row exactly
# (asserted below) or the script stops.
#
# Selection (Alternative: 35% vs 30%): alert fired before look 2 (fc20 < 1,793),
# companion z at alert below the Repair-A bar (z < 3.0) and near the frozen
# median (2.35-2.65), GS stop at look 2 (gs_time == 1,793); among those, the
# trial whose lead (1,793 - fc20) is closest to the subset median; smallest
# sim index breaks ties.
#
# Outputs: figures/fig3_worked_alert.{png,pdf}
#          tables/worked_alert_replay.md   (ALERT_N / ALERT_Z / LEAD_N / LEAD_PCT)

source("scripts/paths.R")
suppressPackageStartupMessages(source(local_ert_file("ertb.R")))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)
fig_dir <- paper_file("figures", root = root, must_work = FALSE)

# frozen design constants (mirror make_ertb_fair_comparators.R)
seed_alt <- 20260542
n_patients <- 2690
p_ctrl <- 0.35
p_trt_alt <- 0.30
look3 <- round(n_patients * (1:3) / 3)      # 897 1793 2690
bounds_2s <- c(3.710, 2.511, 1.993)         # LD-OBF K3 two-sided a=0.05 (pinned)
z_bar <- 3.0                                # Repair A confirmation bar
threshold <- 20

# trial machinery (identical to make_ertb_fair_comparators.R)
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
  (p_c - p_t) / se
}

# ---------------------------------------------------------------------------
# 1. Select the representative trial from the frozen per-trial CSV
# ---------------------------------------------------------------------------

fc <- read.csv(file.path(table_dir, "union_doctrine_first_crossings.csv"))

cand <- fc %>%
  filter(scenario == "Alternative: 35% vs 30%",
         gs_time == look3[2],
         !is.na(fc20), fc20 < look3[2],
         !is.na(z_at_alert), z_at_alert < z_bar,
         z_at_alert >= 2.35, z_at_alert <= 2.65) %>%
  mutate(lead = look3[2] - fc20)

stopifnot(nrow(cand) > 0)
med_lead <- median(cand$lead)
pick <- cand %>%
  arrange(abs(lead - med_lead), sim) %>%
  slice(1)
cat(sprintf("Picked sim %d: alert at %d, z at alert %.3f, lead %d (subset median lead %.0f, n candidates %d)\n",
            pick$sim, pick$fc20, pick$z_at_alert, pick$lead, med_lead, nrow(cand)))

# ---------------------------------------------------------------------------
# 2. Replay the stream to the picked trial (kernels consume no RNG)
# ---------------------------------------------------------------------------

set.seed(seed_alt)
for (s in seq_len(pick$sim)) {
  trial <- simulate_trial(n_patients, p_trt = p_trt_alt, p_ctrl = p_ctrl)
}

w_up <- compute_eRT(trial$treatment, trial$outcome,
                    wager = "design",
                    p_trt_design = p_trt_alt, p_ctrl_design = p_ctrl,
                    ramp_fixed = FALSE)
z_seq <- vapply(seq_len(look3[2]),
                function(n) z_at_look(trial$treatment, trial$outcome, n),
                numeric(1))

# reproduction checks: the replayed trial must match its frozen CSV row
alert_n <- which(w_up >= threshold)[1]
stopifnot(identical(as.integer(alert_n), as.integer(pick$fc20)))
stopifnot(abs(z_seq[alert_n] - pick$z_at_alert) < 1e-9)
stopifnot(abs(z_seq[look3[1]]) < bounds_2s[1])   # no GS stop at look 1
stopifnot(abs(z_seq[look3[2]]) >= bounds_2s[2])  # GS stop at look 2
cat("Reproduction checks passed: replayed trial matches the frozen CSV row.\n")

alert_z <- z_seq[alert_n]
lead_n <- look3[2] - alert_n
lead_pct <- 100 * lead_n / n_patients

# ---------------------------------------------------------------------------
# 3. Values for the manuscript box
# ---------------------------------------------------------------------------

md <- c(
  "# Worked-alert replay (practical-deployment section)",
  "",
  sprintf("Frozen stream: seed %d, sim %d (Alternative: 35%% vs 30%%, N = %d).",
          seed_alt, pick$sim, n_patients),
  "",
  "| Placeholder | Value |",
  "| --- | --- |",
  sprintf("| ALERT_N | %d |", alert_n),
  sprintf("| ALERT_Z | %.2f |", alert_z),
  sprintf("| LEAD_N | %d |", lead_n),
  sprintf("| LEAD_PCT | %.0f%% |", lead_pct),
  "",
  sprintf("Candidates matching the selection window: %d; subset median lead %.0f.",
          nrow(cand), med_lead),
  "Reproduction checks (CSV row, look-1 no-cross, look-2 cross) all passed."
)
writeLines(md, file.path(table_dir, "worked_alert_replay.md"))

# ---------------------------------------------------------------------------
# 4. Figure: wealth path (A) and companion z path (B), aligned on patients
# ---------------------------------------------------------------------------

x_max <- look3[2] + 60
df_w <- data.frame(n = seq_len(look3[2]), w = w_up[seq_len(look3[2])])
df_z <- data.frame(n = seq_len(look3[2]), z = z_seq) %>% filter(is.finite(z), n >= 30)

look_lines <- data.frame(n = look3[1:2], lab = c("look 1", "look 2 - trial stops"))

pA <- ggplot(df_w, aes(n, w)) +
  geom_vline(xintercept = look3[1:2], linetype = "dotted", colour = "#9a9a9a", linewidth = 0.4) +
  geom_hline(yintercept = threshold, linetype = "dashed", colour = "#b23a48", linewidth = 0.45) +
  geom_line(colour = "#0072B2", linewidth = 0.55) +
  geom_point(data = data.frame(n = alert_n, w = w_up[alert_n]),
             colour = "#0072B2", size = 2.2) +
  annotate("text", x = alert_n, y = w_up[alert_n], label = sprintf("alert (patient %d)", alert_n),
           hjust = 1.05, vjust = -0.8, size = 2.9, colour = "#1a1a1a") +
  annotate("text", x = 20, y = threshold, label = "W = 20", hjust = 0, vjust = -0.6,
           size = 2.7, colour = "#b23a48") +
  scale_y_log10() +
  scale_x_continuous(breaks = c(0, 500, look3[1], 1500, look3[2])) +
  coord_cartesian(xlim = c(0, x_max)) +
  labs(x = NULL, y = "Wealth W (log scale)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#ececec", linewidth = 0.3),
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

df_bound <- data.frame(n = look3[1:2], z = bounds_2s[1:2])

pB <- ggplot(df_z, aes(n, z)) +
  geom_vline(xintercept = look3[1:2], linetype = "dotted", colour = "#9a9a9a", linewidth = 0.4) +
  geom_hline(yintercept = z_bar, linetype = "dashed", colour = "#b23a48", linewidth = 0.45) +
  geom_line(colour = "#4d4d4d", linewidth = 0.5) +
  geom_point(data = df_bound, aes(n, z), shape = 21, fill = "white",
             colour = "#E69F00", size = 2.6, stroke = 0.9) +
  geom_point(data = data.frame(n = alert_n, z = alert_z),
             colour = "#0072B2", size = 2.2) +
  annotate("text", x = alert_n, y = alert_z,
           label = sprintf("z = %.2f at alert", alert_z),
           hjust = 1.08, vjust = 2.6, size = 2.9, colour = "#1a1a1a") +
  annotate("text", x = 20, y = z_bar, label = "Repair A bar, z = 3.0", hjust = 0,
           vjust = -0.6, size = 2.7, colour = "#b23a48") +
  annotate("text", x = look3[1], y = bounds_2s[1], label = "scheduled boundary",
           hjust = -0.08, vjust = 0.4, size = 2.7, colour = "#8a5b00") +
  scale_x_continuous(breaks = c(0, 500, look3[1], 1500, look3[2])) +
  coord_cartesian(xlim = c(0, x_max)) +
  labs(x = "Patients enrolled", y = "One-sided z") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#ececec", linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA))

p <- pA / pB + plot_annotation(tag_levels = "A")

ggsave(file.path(fig_dir, "fig3_worked_alert.png"), p, width = 7.2, height = 5.6,
       dpi = 360, bg = "white")
ggsave(file.path(fig_dir, "fig3_worked_alert.pdf"), p, width = 7.2, height = 5.6,
       device = cairo_pdf, bg = "white")
cat("Wrote fig3_worked_alert.png/.pdf and tables/worked_alert_replay.md\n")
