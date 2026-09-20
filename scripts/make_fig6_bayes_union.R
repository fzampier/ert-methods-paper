# Figure 6: Bayesian/e-RT union calibration (carried exhibit, rebuilt).
# A scheduled Bayesian rule (independent Beta(1,1) priors, posterior threshold
# checked at the three K=3 looks) alone, and combined with design e-RTb as an
# automatic union, across posterior thresholds 0.975-0.995. Panel A:
# alternative signaling. Panel B: null signaling, y-axis 0-8% with a dashed 5%
# reference (same convention as Figure 4B).
# Deterministic post-processing of frozen summaries -- no rerun, no RNG.
#
# Inputs : tables/bayesian_monitor_method_summary.csv,
#          tables/bayesian_ert_union_summary.csv
# Output : figures/fig6_bayes_union.{png,pdf}

source("scripts/paths.R")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)
fig_dir <- paper_file("figures", root = root, must_work = FALSE)

thresholds <- c(0.975, 0.98, 0.985, 0.99, 0.995)

alone <- read.csv(file.path(table_dir, "bayesian_monitor_method_summary.csv")) %>%
  filter(threshold %in% thresholds) %>%
  transmute(scenario, threshold, series = "Bayesian rule alone",
            rate = interrupted_any)

union <- read.csv(file.path(table_dir, "bayesian_ert_union_summary.csv")) %>%
  filter(threshold %in% thresholds,
         ert_method == "e-RTb design directional") %>%
  transmute(scenario, threshold,
            series = "Bayesian rule + design e-RTb (automatic union)",
            rate = combined_signal)

panel_lab <- c("Alternative: 35% vs 30%" = "A. Alternative: 35% vs 30%",
               "Null: 35% vs 35%" = "B. Null: 35% vs 35%")

df <- bind_rows(alone, union) %>%
  filter(scenario %in% names(panel_lab)) %>%
  mutate(panel = panel_lab[scenario],
         series = factor(series, levels = c(
           "Bayesian rule + design e-RTb (automatic union)",
           "Bayesian rule alone")))

pal <- c("Bayesian rule + design e-RTb (automatic union)" = "#D55E00",
         "Bayesian rule alone" = "#000000")
lty <- c("Bayesian rule + design e-RTb (automatic union)" = "solid",
         "Bayesian rule alone" = "dashed")

ref <- data.frame(panel = panel_lab[2], y = 0.05)
caps <- data.frame(panel = panel_lab, cap = c(1, 0.08))

p <- ggplot(df, aes(threshold, rate, colour = series, linetype = series)) +
  geom_hline(data = ref, aes(yintercept = y), linetype = "22",
             linewidth = 0.35, colour = "grey45") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8, shape = 16, show.legend = FALSE) +
  geom_blank(data = caps, aes(y = cap), inherit.aes = FALSE) +
  geom_blank(data = caps, aes(y = 0), inherit.aes = FALSE) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(values = lty, name = NULL) +
  scale_x_continuous(breaks = thresholds) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  labs(x = "Posterior probability threshold", y = "Signaling probability") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    strip.text = element_text(face = "bold", hjust = 0, size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.margin = margin(4, 8, 4, 4)
  )

ggsave(file.path(fig_dir, "fig6_bayes_union.png"), p,
       width = 180, height = 85, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "fig6_bayes_union.pdf"), p,
       width = 180, height = 85, units = "mm", device = cairo_pdf)

cat("fig6 written\n")
