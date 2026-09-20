# Figure 4: cumulative first-crossing (alert) curves, fair-comparator run.
# Primary method set only (Table 1's first seven rows): calibrated one-sided
# LD-OBF K3/K5, KD rho=1, calibrated two-sided LD-OBF K3, and the three e-RTb
# wager policies. Panel A: alternative (35% vs 30%). Panel B: null (35% vs 35%),
# y-axis 0-8% with a dashed 5% reference line (a 0-5% axis would clip nothing in
# this method set, but the 8% ceiling keeps the panel comparable with the
# worst-case union rate discussed in the text).
# Post-processing of the frozen per-trial CSV -- no rerun, no RNG.
#
# Input : tables/fair_comparators_first_crossings.csv
# Output: figures/fig4_crossing_curves.{png,pdf}
#
# Okabe-Ito palette; adjacent-pair CVD separation checked (OKLab distance under
# Machado protan/deutan/tritan simulation; all pairs >= 8, normal >= 15).
# Scheduled designs draw as step functions (they can only cross at their looks);
# e-RTb draws as a dense curve (it can cross at any patient).

source("scripts/paths.R")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)
fig_dir <- paper_file("figures", root = root, must_work = FALSE)

cross <- read.csv(file.path(table_dir, "fair_comparators_first_crossings.csv"))
N <- 2690
n_sims <- 5000

methods <- c(
  "e-RTb design directional"     = "e-RTb design",
  "e-RTb mixture two-sided"      = "e-RTb mixture",
  "e-RTb adaptive"               = "e-RTb adaptive",
  "LD-OBF K3 one-sided a=0.05"   = "LD-OBF K=3 one-sided",
  "LD-OBF K5 one-sided a=0.05"   = "LD-OBF K=5 one-sided",
  "KD rho=1 K5 one-sided a=0.05" = "KD ρ=1 K=5 one-sided",
  "LD-OBF K3 two-sided a=0.05"   = "LD-OBF K=3 two-sided"
)
scenarios <- c(
  "Alternative: 35% vs 30%" = "A. Alternative: 35% vs 30%",
  "Null: 35% vs 35%"        = "B. Null: 35% vs 35%"
)

pal <- c(
  "e-RTb design"           = "#D55E00",
  "e-RTb mixture"          = "#CC79A7",
  "e-RTb adaptive"         = "#E69F00",
  "LD-OBF K=3 one-sided"   = "#0072B2",
  "LD-OBF K=5 one-sided"   = "#56B4E9",
  "KD ρ=1 K=5 one-sided" = "#009E73",
  "LD-OBF K=3 two-sided"   = "#000000"
)
lty <- c(
  "e-RTb design"           = "solid",
  "e-RTb mixture"          = "longdash",
  "e-RTb adaptive"         = "dotdash",
  "LD-OBF K=3 one-sided"   = "solid",
  "LD-OBF K=5 one-sided"   = "dashed",
  "KD ρ=1 K=5 one-sided" = "solid",
  "LD-OBF K=3 two-sided"   = "dashed"
)
is_ert <- function(m) grepl("^e-RTb", m)

# Cumulative curve on a common patient grid (right-continuous step values).
curve_for <- function(d) {
  x <- sort(d$first_crossing[!is.na(d$first_crossing)])
  grid <- 0:N
  data.frame(patient = grid,
             cum = findInterval(grid, x) / n_sims)
}

df <- cross %>%
  filter(scenario %in% names(scenarios), method %in% names(methods)) %>%
  mutate(panel = scenarios[scenario], label = methods[method]) %>%
  group_by(panel, label) %>%
  group_modify(~ curve_for(.x)) %>%
  ungroup() %>%
  mutate(label = factor(label, levels = unname(methods)),
         family = ifelse(is_ert(as.character(label)), "e-RTb", "scheduled"))

# Free y per panel via facet; cap the null panel at 8% with a 5% line.
ref <- data.frame(panel = scenarios[2], y = 0.05)

p <- ggplot(df, aes(patient, cum, colour = label, linetype = label)) +
  geom_hline(data = ref, aes(yintercept = y), linetype = "22",
             linewidth = 0.35, colour = "grey45") +
  geom_step(data = ~ filter(.x, family == "scheduled"), linewidth = 0.55) +
  geom_line(data = ~ filter(.x, family == "e-RTb"), linewidth = 0.8) +
  facet_wrap(~panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(values = lty, name = NULL) +
  scale_x_continuous(breaks = c(0, 897, 1793, 2690), expand = c(0.01, 0)) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  labs(x = "Patient", y = "Cumulative crossing (alert) probability") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    strip.text = element_text(face = "bold", hjust = 0, size = 9),
    axis.title = element_text(size = 9),
    legend.position = "bottom",
    legend.key.width = unit(1.6, "lines"),
    legend.text = element_text(size = 8),
    plot.margin = margin(4, 8, 4, 4)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# Null panel ceiling: 8%.
panel_limits <- df %>%
  mutate(cap = ifelse(panel == scenarios[2], 0.08, 1))
p <- p + geom_blank(data = panel_limits, aes(y = cap))

ggsave(file.path(fig_dir, "fig4_crossing_curves.png"), p,
       width = 180, height = 95, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "fig4_crossing_curves.pdf"), p,
       width = 180, height = 95, units = "mm", device = cairo_pdf)

cat("fig4 written; methods:", paste(unname(methods), collapse = " | "), "\n")
