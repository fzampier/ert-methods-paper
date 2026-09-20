# Figure 5 (round-2 rebuild, 2026-09-17): panel B relabelled 'expected time to
# first boundary crossing' per Reviewer 4 (round 2). Same frozen inputs, same
# deterministic post-processing.
# Panel A: power of exactly calibrated LD-OBF error-spending designs as K grows
# (analytic values from power_vs_K.csv; zero MC error), one- and two-sided at
# alpha=0.05, with the e-RTb crossing probabilities as horizontal references.
# Panel B: expected stopping time E[min(tau,N)] under the alternative, same
# layout. e-RT levels are recomputed from the frozen per-trial CSV (not typed
# in), so the figure cannot drift from the tables.
# Deterministic post-processing -- no rerun, no RNG.
#
# Inputs : tables/power_vs_K.csv, tables/fair_comparators_first_crossings.csv
# Output : figures/fig5_price_anytime.{png,pdf}
#
# Colors follow the same entities as Figure 4 (Okabe-Ito).

source("scripts/paths.R")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)
fig_dir <- paper_file("figures", root = root, must_work = FALSE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

N <- 2690

pk <- read.csv(file.path(table_dir, "power_vs_K.csv")) %>%
  mutate(series = ifelse(sided == 1, "LD-OBF one-sided", "LD-OBF two-sided"))

ert_methods <- c(
  "e-RTb design directional" = "e-RTb design",
  "e-RTb mixture two-sided"  = "e-RTb mixture",
  "e-RTb adaptive"           = "e-RTb adaptive"
)
ert <- read.csv(file.path(table_dir, "fair_comparators_first_crossings.csv")) %>%
  filter(scenario == "Alternative: 35% vs 30%", method %in% names(ert_methods)) %>%
  group_by(series = ert_methods[method]) %>%
  summarise(power = mean(!is.na(first_crossing)),
            expected_n = mean(pmin(first_crossing, N, na.rm = FALSE) %>%
                                ifelse(is.na(.), N, .)),
            .groups = "drop")

levels_all <- c("LD-OBF one-sided", "LD-OBF two-sided",
                "e-RTb design", "e-RTb mixture", "e-RTb adaptive")
pal <- c("LD-OBF one-sided" = "#0072B2", "LD-OBF two-sided" = "#000000",
         "e-RTb design" = "#D55E00", "e-RTb mixture" = "#CC79A7",
         "e-RTb adaptive" = "#E69F00")
lty <- c("LD-OBF one-sided" = "solid", "LD-OBF two-sided" = "solid",
         "e-RTb design" = "solid", "e-RTb mixture" = "longdash",
         "e-RTb adaptive" = "dotdash")

long_pk <- pk %>%
  select(K, series, Power = power, `Expected stopping time` = expected_n) %>%
  pivot_longer(c(Power, `Expected stopping time`),
               names_to = "panel", values_to = "value")
long_ert <- ert %>%
  select(series, Power = power, `Expected stopping time` = expected_n) %>%
  pivot_longer(c(Power, `Expected stopping time`),
               names_to = "panel", values_to = "value")

panel_lab <- c("Power" = "A. Power under the alternative",
               "Expected stopping time" = "B. Expected time to first crossing (alternative)")
long_pk$panel <- panel_lab[long_pk$panel]
long_ert$panel <- panel_lab[long_ert$panel]
long_pk$series <- factor(long_pk$series, levels = levels_all)
long_ert$series <- factor(long_ert$series, levels = levels_all)

fmt_y <- function(panel) if (grepl("Power", panel)) scales::percent else scales::comma

p <- ggplot(long_pk, aes(K, value, colour = series, linetype = series)) +
  geom_hline(data = long_ert, aes(yintercept = value, colour = series,
                                  linetype = series), linewidth = 0.7) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8, shape = 16, show.legend = FALSE) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = lty, name = NULL, drop = FALSE) +
  scale_x_continuous(trans = "log2", breaks = c(1, 2, 3, 5, 10, 20)) +
  labs(x = "Number of scheduled looks K (log scale)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    strip.text = element_text(face = "bold", hjust = 0, size = 9),
    legend.position = "bottom",
    legend.key.width = unit(1.6, "lines"),
    legend.text = element_text(size = 8),
    plot.margin = margin(4, 8, 4, 4)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# Percent labels for panel A, comma labels for panel B: format after faceting
# via per-panel scales is not native in ggplot2; use labeller-free approach with
# facetted_pos_scales if available, else keep raw proportions for A.
if (requireNamespace("ggh4x", quietly = TRUE)) {
  p <- p + ggh4x::facetted_pos_scales(y = list(
    panel == panel_lab["Power"] ~
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)),
    panel == panel_lab["Expected stopping time"] ~
      scale_y_continuous(labels = scales::comma)
  ))
} else {
  p <- p + scale_y_continuous(labels = function(x)
    ifelse(x <= 1, sprintf("%.0f%%", 100 * x), scales::comma(x)))
}

ggsave(file.path(fig_dir, "fig5_price_anytime.png"), p,
       width = 180, height = 85, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "fig5_price_anytime.pdf"), p,
       width = 180, height = 85, units = "mm", device = cairo_pdf)

cat("fig5 written; e-RT levels:\n"); print(as.data.frame(ert))
