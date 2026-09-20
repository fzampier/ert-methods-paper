# Run every target from the repository root.
#
# `make check` reads the frozen outputs only; it reruns nothing. The other
# targets rewrite files under tables/ and figures/ in place, so `git status`
# afterwards shows any drift from the frozen record.

.PHONY: all check kernels-check manifest-check numbers-check figures postprocess simulations trials session-info

# Default = the checks only. NEVER a regeneration target.
all: check

check: kernels-check manifest-check numbers-check

kernels-check:
	Rscript R/check_kernels.R

manifest-check:
	@if command -v sha256sum >/dev/null 2>&1; then sha256sum --quiet -c MANIFEST.sha256; else shasum -a 256 -q -c MANIFEST.sha256; fi
	@echo "MANIFEST.sha256: frozen tables, figures, and supplementary code are intact."

numbers-check:
	Rscript R/check_paper_numbers.R

# Figures 1-6. Deterministic: fixed display seeds or frozen CSV inputs.
figures:
	python3 scripts/make_fig1_infoflow.py
	Rscript scripts/make_fig2_family_trajectories.R
	Rscript scripts/make_worked_alert_replay.R
	Rscript scripts/make_fig4_crossing_curves.R
	Rscript scripts/make_fig5_price_anytime.R
	Rscript scripts/make_fig6_bayes_union.R

# Tables that are post-processing of frozen per-trial outputs (seconds, no RNG;
# make_power_vs_K.R carries ~1e-6 quadrature jitter, see README).
postprocess:
	Rscript scripts/make_cumulative_at_looks.R
	Rscript scripts/make_union_path_decomposition.R
	Rscript scripts/make_power_vs_K.R

# Every simulation, seeds fixed in the scripts. About 1.5 hours in sequence;
# the scripts are independent of each other and can run in parallel
# (fair comparators first if you also want to rebuild its post-processing).
simulations:
	Rscript scripts/make_ertb_section_3_tables.R
	Rscript scripts/make_blocked_randomization_sensitivity.R
	Rscript scripts/make_ertb_unequal_allocation_check.R
	Rscript scripts/make_erte_censoring_null.R
	Rscript scripts/make_erte_power_check.R
	Rscript scripts/make_ertc_type1_simulation.R
	Rscript scripts/make_ertc_power_check.R
	Rscript scripts/make_ertb_interruption_curve.R
	Rscript scripts/make_pocock_matched_power.R
	Rscript scripts/make_ertb_fair_comparators.R
	Rscript scripts/make_batched_updating_sensitivity.R
	Rscript scripts/make_ertb_design_actual_power_table.R
	Rscript scripts/make_bayesian_monitor_sensitivity.R
	Rscript scripts/make_premium_decomposition.R
	Rscript scripts/make_pi_out_of_order_check.R

# ESM Figures S1-S2 and Tables S21-S22. Needs the BioLINCC-derived files
# (data/README.md).
trials:
	Rscript scripts/make_arma_permutation_figure.R
	Rscript scripts/make_factt_permutation_figure.R

session-info:
	Rscript R/capture_session_info.R
