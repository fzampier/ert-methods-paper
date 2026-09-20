# NEW Figure 2 (author redesign 2026-08-10): six wealth-trajectory panels for
# the e-RT family, null and effect scenario per flavour, uniform anchors:
#   A/B  e-RTb  binary, control 40%: null vs 5 pp benefit; N = 2,942
#               (chi-square N for 80% power, as in ESM Tables S1-S2)
#   C/D  e-RTe  the events of the same trials (updates at events only)
#   E/F  e-RTc  continuous: null vs +0.2 SD shift; N = 788 (t-test N, 80%)
# 30 paths per panel, adaptive wagers at shipped defaults.
# DISPLAY seed only (figures are deterministic replay under the freeze).
# Output: figures/fig2_family.{png,pdf}

source("scripts/paths.R")
suppressPackageStartupMessages({
  source(local_ert_file("ertb.R"))
  source(local_ert_file("erte.R"))
  source(local_ert_file("ertc.R"))
  library(ggplot2)
})

root <- paper_root()
fig_dir <- paper_file("figures", root = root, must_work = FALSE)

set.seed(20260552)  # display seed
N_PATHS <- 30L
THRESHOLD <- 20
N_B <- 2942L; P_CTRL <- 0.40; P_TRT_EFF <- 0.35
N_C <- 788L; DELTA <- 0.2

sim_ertb <- function(p_trt) {
  treatment <- rbinom(N_B, 1, 0.5)
  rate <- ifelse(treatment == 1, p_trt, P_CTRL)
  outcome <- rbinom(N_B, 1, rate)
  compute_eRT(treatment, outcome)
}
sim_erte <- function(p_trt) {
  pooled <- (P_CTRL + p_trt) / 2
  q <- p_trt / (p_trt + P_CTRL)
  d <- rbinom(1, N_B, pooled)
  arms <- rbinom(d, 1, q)
  compute_eRTe(arms, wager = "adaptive")$wealth
}
sim_ertc <- function(delta) {
  treatment <- rbinom(N_C, 1, 0.5)
  outcome <- rnorm(N_C, mean = delta * treatment, sd = 1)
  # conservative setting (also used for the FACTT illustration): paths persist
  # over the whole stream; the default sign policy drains far below the
  # display range after its early peaks (volatility drag) and panels go empty
  compute_eRTc(treatment, outcome, burn_in = 50L, ramp = 100L, c_max = 0.6,
               wager = "adaptive", adaptive_direction = "cohens_d")
}

panels <- list(
  list(key = "A", label = "A   e-RTb: no effect",             sim = function() sim_ertb(P_CTRL)),
  list(key = "B", label = "B   e-RTb: 5 pp benefit",          sim = function() sim_ertb(P_TRT_EFF)),
  list(key = "C", label = "C   e-RTe: no effect",             sim = function() sim_erte(P_CTRL)),
  list(key = "D", label = "D   e-RTe: 5 pp benefit",          sim = function() sim_erte(P_TRT_EFF)),
  list(key = "E", label = "E   e-RTc: no effect",             sim = function() sim_ertc(0)),
  list(key = "F", label = "F   e-RTc: +0.2 SD shift",         sim = function() sim_ertc(DELTA))
)

rows <- list(); ri <- 1L
for (pn in panels) {
  n_crossed <- 0L
  for (path in seq_len(N_PATHS)) {
    w <- pn$sim()
    crossed <- any(w >= THRESHOLD)
    n_crossed <- n_crossed + crossed
    rows[[ri]] <- data.frame(
      panel = pn$label, path = paste0(pn$key, path),
      update = seq_along(w), wealth = pmax(w, 1e-6),
      status = if (crossed) "Crossed" else "Did not cross"
    )
    ri <- ri + 1L
  }
  cat(sprintf("%s: %d/%d crossed\n", pn$label, n_crossed, N_PATHS))
}
df <- do.call(rbind, rows)
df$panel <- factor(df$panel, levels = vapply(panels, `[[`, "", "label"))
df$status <- factor(df$status, levels = c("Did not cross", "Crossed"))
df <- df[order(df$status), ]  # crossed paths drawn on top

lab_df <- data.frame(panel = factor(vapply(panels, `[[`, "", "label"),
                                    levels = levels(df$panel)),
                     x = 1, y = THRESHOLD)

p <- ggplot(df, aes(update, wealth, group = path, colour = status)) +
  geom_hline(yintercept = THRESHOLD, linetype = "dashed", colour = "#b23a48", linewidth = 0.45) +
  geom_line(aes(linewidth = status, alpha = status)) +
  geom_text(data = lab_df, aes(x = x, y = y), label = "W = 20", inherit.aes = FALSE,
            hjust = 0, vjust = -0.6, size = 2.7, colour = "#b23a48") +
  scale_colour_manual(values = c("Did not cross" = "#b8b8b8", "Crossed" = "#0072B2"), name = NULL) +
  scale_linewidth_manual(values = c("Did not cross" = 0.3, "Crossed" = 0.5), guide = "none") +
  scale_alpha_manual(values = c("Did not cross" = 0.65, "Crossed" = 0.9), guide = "none") +
  scale_y_log10(limits = c(1e-2, NA)) +
  facet_wrap(~panel, ncol = 2, scales = "free_x") +
  labs(x = "Update (patient for e-RTb/e-RTc; event for e-RTe)",
       y = "Wealth W (log scale)") +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.9, alpha = 1))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#ececec", linewidth = 0.3),
        strip.text = element_text(face = "bold", hjust = 0, size = 9),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(fig_dir, "fig2_family.png"), p, width = 8.6, height = 8.8, dpi = 360, bg = "white")
ggsave(file.path(fig_dir, "fig2_family.pdf"), p, width = 8.6, height = 8.8, device = cairo_pdf, bg = "white")
cat("Wrote fig2_family.png/.pdf\n")
