# Path attribution for the union-doctrine rules (post-processing; no rerun).
#
# Answers the author's question (2026-08-08): Repair A's power is carried by
# the GS component -- so decompose every stop by WHICH path caused it, report
# the earliness delivered by the confirmed-alert path, count alerts whose
# confirmation failed, and quantify the pure surveillance value: among trials
# the scheduled monitor stopped, how often had e-RT already alerted earlier,
# and by how many patients.  Note the z>=3.0 check at the alert is not an
# independent "agreeing method": it is evaluated at that moment BECAUSE the
# alert fired; attribution, not agreement, is the meaningful bookkeeping.
#
# Input : tables/union_doctrine_first_crossings.csv  (from make_ertb_fair_comparators.R)
# Output: tables/union_path_decomposition.md / .csv
#         tables/alert_lead_time.md   (role-1 surveillance value)

source("scripts/paths.R")
root <- paper_root()
table_dir <- paper_file("tables", root = root, must_work = FALSE)

u <- read.csv(file.path(table_dir, "union_doctrine_first_crossings.csv"))

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)
fmt_n <- function(x) ifelse(is.na(x), "--", format(round(x), big.mark = ","))

decompose_repair_a <- function(d, zthr) {
  n <- nrow(d)
  alert <- !is.na(d$fc20)
  confirmed <- alert & is.finite(d$z_at_alert) & d$z_at_alert >= zthr
  conf_time <- ifelse(confirmed, d$fc20, NA)
  gs <- d$gs_time

  stop_any <- confirmed | !is.na(gs)
  # attribution: confirmed-alert path strictly earlier than the scheduled stop
  # (or scheduled path never stops)
  via_alert <- confirmed & (is.na(gs) | conf_time < gs)
  via_gs <- !via_alert & !is.na(gs)

  # earliness delivered by the alert path in the trials it won
  lead <- (gs - conf_time)[via_alert & !is.na(gs)]

  data.frame(
    rule = sprintf("Repair A z>=%.1f", zthr),
    stop_rate = mean(stop_any),
    via_confirmed_alert = mean(via_alert),
    via_scheduled_gs = mean(via_gs),
    alert_fired = mean(alert),
    alert_unconfirmed = mean(alert & !confirmed),
    alert_won_median_lead = ifelse(any(via_alert & !is.na(gs)), median(lead), NA),
    alert_won_vs_no_gs_stop = mean(via_alert & is.na(gs))
  )
}

decompose_repair_b <- function(d, wcol, gcol, label) {
  trig <- d[[wcol]]
  gs <- d[[gcol]]
  stop_any <- !is.na(trig) | !is.na(gs)
  via_trig <- !is.na(trig) & (is.na(gs) | trig < gs)
  via_gs <- !via_trig & !is.na(gs)
  lead <- (gs - trig)[via_trig & !is.na(gs)]
  data.frame(
    rule = label,
    stop_rate = mean(stop_any),
    via_confirmed_alert = mean(via_trig),
    via_scheduled_gs = mean(via_gs),
    alert_fired = mean(!is.na(trig)),
    alert_unconfirmed = NA,
    alert_won_median_lead = ifelse(any(via_trig & !is.na(gs)), median(lead), NA),
    alert_won_vs_no_gs_stop = mean(via_trig & is.na(gs))
  )
}

# Repair B needs the recalibrated GS times, which were not stored as columns;
# reconstruct conservatively from stored stop times: repair_bXXX = min(gsRecal,
# fcXXX), so the trigger won exactly when repair_bXXX == fcXXX and (repair time
# < stored default-gs proxy is NOT valid). Instead attribute directly:
# trigger won iff !is.na(fcXXX) and repair_bXXX == fcXXX and
# (fcXXX != scheduled look time OR repair equals fc at a non-look index).
# Cleanest available: via_trigger := !is.na(fc) & repair == fc & (is.na(gs_recal_proxy) | fc <= gs_recal_proxy)
# We recover gs_recal times: repair_b == fc when trigger first, else repair_b == gs_recal.
recover_b <- function(d, wcol, rcol) {
  fc <- d[[wcol]]; rp <- d[[rcol]]
  via_trig <- !is.na(rp) & !is.na(fc) & rp == fc
  gs_recal <- ifelse(!is.na(rp) & (is.na(fc) | rp != fc), rp, NA)  # stop time when GS-recal path won
  list(via_trig = via_trig, gs_recal_won_time = gs_recal)
}

out <- list()
for (sc in unique(u$scenario)) {
  d <- u[u$scenario == sc, ]
  a30 <- decompose_repair_a(d, 3.0)
  a25 <- decompose_repair_a(d, 2.5)

  b100r <- recover_b(d, "fc100", "repair_b100")
  b200r <- recover_b(d, "fc200", "repair_b200")
  b100 <- data.frame(
    rule = "Repair B W>=100 (+GS a=0.04)",
    stop_rate = mean(!is.na(d$repair_b100)),
    via_confirmed_alert = mean(b100r$via_trig),
    via_scheduled_gs = mean(!is.na(d$repair_b100)) - mean(b100r$via_trig),
    alert_fired = mean(!is.na(d$fc100)),
    alert_unconfirmed = NA,
    alert_won_median_lead = NA,
    alert_won_vs_no_gs_stop = NA
  )
  b200 <- data.frame(
    rule = "Repair B W>=200 (+GS a=0.045)",
    stop_rate = mean(!is.na(d$repair_b200)),
    via_confirmed_alert = mean(b200r$via_trig),
    via_scheduled_gs = mean(!is.na(d$repair_b200)) - mean(b200r$via_trig),
    alert_fired = mean(!is.na(d$fc200)),
    alert_unconfirmed = NA,
    alert_won_median_lead = NA,
    alert_won_vs_no_gs_stop = NA
  )
  block <- rbind(a30, a25, b100, b200)
  block$scenario <- sc
  out[[sc]] <- block
}
dec <- do.call(rbind, out)
write.csv(dec, file.path(table_dir, "union_path_decomposition.csv"), row.names = FALSE)

md <- data.frame(
  Scenario = dec$scenario,
  Rule = dec$rule,
  `Stop rate` = fmt_pct(dec$stop_rate),
  `Stops via alert path` = fmt_pct(dec$via_confirmed_alert),
  `Stops via scheduled GS` = fmt_pct(dec$via_scheduled_gs),
  `Alert fired` = fmt_pct(dec$alert_fired),
  `Alert fired, unconfirmed` = ifelse(is.na(dec$alert_unconfirmed), "--", fmt_pct(dec$alert_unconfirmed)),
  `Median lead when alert won (patients)` = fmt_n(dec$alert_won_median_lead),
  check.names = FALSE
)
writeLines(c(
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |"),
  apply(md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |"))
), file.path(table_dir, "union_path_decomposition.md"))

# ---- Role-1 surveillance value: alert precedes the scheduled stop ----
lead_rows <- lapply(unique(u$scenario), function(sc) {
  d <- u[u$scenario == sc, ]
  gs_stopped <- !is.na(d$gs_time)
  pre <- gs_stopped & !is.na(d$fc20) & d$fc20 < d$gs_time
  lead <- (d$gs_time - d$fc20)[pre]
  data.frame(
    Scenario = sc,
    `GS-stopped trials` = fmt_pct(mean(gs_stopped)),
    `...preceded by e-RT alert` = fmt_pct(ifelse(any(gs_stopped), mean(pre[gs_stopped]), NA)),
    `Median lead (patients)` = fmt_n(ifelse(any(pre), median(lead), NA)),
    `Q1--Q3 lead` = if (any(pre)) sprintf("%s--%s", fmt_n(quantile(lead, .25, names = FALSE)),
                                          fmt_n(quantile(lead, .75, names = FALSE))) else "--",
    `Median z at alert` = sprintf("%.2f", median(d$z_at_alert[pre], na.rm = TRUE)),
    check.names = FALSE
  )
})
lead_md <- do.call(rbind, lead_rows)
writeLines(c(
  paste0("| ", paste(names(lead_md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(lead_md)), collapse = " | "), " |"),
  apply(lead_md, 1, function(r) paste0("| ", paste(as.character(r), collapse = " | "), " |"))
), file.path(table_dir, "alert_lead_time.md"))

cat("=== Path decomposition ===\n"); print(md, row.names = FALSE)
cat("\n=== Role-1 surveillance value (design e-RT W>=20 vs scheduled LD-OBF K3 2s) ===\n")
print(lead_md, row.names = FALSE)
