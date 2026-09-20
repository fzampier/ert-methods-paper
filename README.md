# e-RT: code and frozen results for the article

Code, simulation outputs, and figures for

> Zampieri FG, Albuquerque AM, Ramdas A. Sequential randomization tests using
> e-values for clinical trial monitoring across endpoint types.
> (Under review. The citation and DOI will be added on publication.)

e-RT is a randomization-based e-value method for continuous trial monitoring.
The article describes three variants: e-RTb for binary endpoints, e-RTe for
event-only streams, and e-RTc for continuous endpoints. This repository holds
the three kernels, every simulation and figure script behind the article and
its Electronic Supplementary Material (ESM), and the frozen outputs those
scripts produced.

## Layout

| Path | Content |
| --- | --- |
| `R/ertb.R`, `R/erte.R`, `R/ertc.R` | The three kernels. Every script sources them from here. |
| `scripts/` | One script per simulation, table, or figure. `scripts/paths.R` holds the path helpers. |
| `tables/` | Frozen outputs: CSV files and the markdown tables printed in the article and ESM. |
| `figures/` | Figures 1 to 6 and ESM Figures S1 and S2 (PNG and PDF). |
| `supplementary_code/` | The supplementary code archive of the article, as submitted to the journal. |
| `R/check_kernels.R`, `R/check_paper_numbers.R` | The checks behind `make check`. |
| `data/README.md` | How to obtain the ARMA and FACTT data, and where the scripts expect them. |
| `MANIFEST.sha256` | SHA-256 of every frozen table, figure, and supplementary file. |
| `sessionInfo.txt` | R and package versions of the frozen runs. |
| `CITATION.cff` | Citation metadata for this repository and its archived releases. |

## Requirements

R 4.6.1 with `tidyverse`, `rpact`, `mvtnorm`, `patchwork`, `ggh4x`, `scales`,
and `Rcpp` with a C++ compiler (FACTT script only); `jsonlite` for the trial
pins of the number check. Exact versions are in `sessionInfo.txt`. Figure 1 needs Python 3 with
`matplotlib`.

## Checks

```sh
make check
```

`make check` reruns nothing. It reads the frozen outputs and verifies three
things:

- the kernels in `R/` match the supplementary code of the article
  (`R/check_kernels.R`);
- every frozen table, figure, and supplementary file matches
  `MANIFEST.sha256`;
- the numbers quoted in the article and ESM agree with the frozen tables
  (`R/check_paper_numbers.R`). The pins that need the ARMA and FACTT data are
  reported as skipped when the data are absent.

## Rerunning

Run everything from the repository root. The scripts rewrite `tables/` and
`figures/` in place, so `git status` afterwards shows any drift from the frozen
record.

```sh
make figures       # Figures 1 to 6 (deterministic)
make postprocess   # tables that post-process frozen per-trial outputs
make simulations   # every simulation, seeds fixed in the scripts
make trials        # ESM Figures S1-S2, Tables S21-S22 (needs the BioLINCC data)
```

Every seed is fixed in its script. The kernels consume no random numbers, so
scripts that share a seed replay the same simulated trials: the batched-updating
tables, the worked alert, and the premium decomposition all replay the streams
of the main comparison (seeds 20260542 and 20260543) and check that they
reproduce its frozen rows.

`scripts/make_power_vs_K.R` integrates multivariate normal probabilities with
`mvtnorm`, whose default algorithm uses randomized quadrature. Reruns move
`tables/power_vs_K.csv` in the sixth decimal of power (up to 4e-5 in the
achieved-alpha column and 0.01 patients in the expected stopping time at 20
looks). The rounded table `tables/power_vs_K.md` and every pinned number stay
the same.

Runtimes on one core of an Apple-silicon desktop: 21 to 24 minutes for each of
the two e-RTc simulations, 7 to 9 minutes for the blocked-randomization and the
design-versus-actual simulations, and under 3 minutes for every other script;
about 80 minutes for everything in sequence.

## Verification of this copy

On 2026-09-19 all 26 scripts were rerun from a copy of this repository (R 4.6.1,
macOS, Apple silicon), and `make check`, `make figures`, and `make trials` were
run from a fresh clone. All 59 files under `tables/` came back byte for byte,
except `tables/power_vs_K.csv` (the quadrature jitter described above). All PNG
figures came back byte for byte from the frozen inputs. PDF figures differ in
bytes on every run while their PNG versions are identical.

## Where each result comes from

Main text:

| Item | Script | Seeds | Frozen output |
| --- | --- | --- | --- |
| Figure 1 | `make_fig1_infoflow.py` | none (drawing) | `figures/fig1_infoflow` |
| Figure 2 | `make_fig2_family_trajectories.R` | 20260552 (display) | `figures/fig2_family` |
| Figure 3 | `make_worked_alert_replay.R` | replays one trial of stream 20260542 | `figures/fig3_worked_alert`, `tables/worked_alert_replay.md` |
| Tables 1 and 2 | `make_ertb_fair_comparators.R` | 20260542 (alternative), 20260543 (null), 20260544 (harm) | `tables/fair_comparators_*`, `tables/union_doctrine*`, `tables/detectable_directions.md` |
| Figure 4 | `make_fig4_crossing_curves.R` | post-processing | `figures/fig4_crossing_curves` |
| Figure 5 | `make_power_vs_K.R`, then `make_fig5_price_anytime.R` | analytic | `tables/power_vs_K.*`, `figures/fig5_price_anytime` |
| Figure 6 | `make_bayesian_monitor_sensitivity.R`, then `make_fig6_bayes_union.R` | 20260514; scenarios 20260552, 20260553 | `tables/bayesian_*`, `figures/fig6_bayes_union` |

ESM:

| Item | Script | Seeds | Frozen output |
| --- | --- | --- | --- |
| Tables S1 to S3 | `make_ertb_section_3_tables.R` | 20260502 | `tables/ertb_section3_*` |
| Table S4 | `make_blocked_randomization_sensitivity.R` | 20260545 | `tables/ertb_blocked_randomization_sensitivity.*` |
| Table S5 | `make_ertb_unequal_allocation_check.R` | 20260549 | `tables/ertb_unequal_allocation_check.*` |
| Table S6 | `make_erte_censoring_null.R` | 20260548 | `tables/erte_censoring_null.*` |
| Table S7 | `make_erte_power_check.R` | 20260550 | `tables/erte_power_check.*` |
| Tables S8 and S10 | `make_ertc_type1_simulation.R` | 20260518; 20260547 (cap rows) | `tables/ertc_type1_simulation.*` |
| Table S9 | `make_ertc_power_check.R` | 20260551 | `tables/ertc_power_check.*` |
| Tables S11 and S12 | `make_batched_updating_sensitivity.R` | streams 20260542, 20260543 | `tables/batched_updating_sensitivity.*`, `tables/ascertainment_delay.*` |
| Tables S13 and S14 | `make_ertb_fair_comparators.R` | as above | `tables/fair_comparators_boundaries.md`, `tables/fair_comparators_table2.md` |
| Table S15 | `make_cumulative_at_looks.R` | post-processing | `tables/fair_comparators_cumulative_at_looks.*` |
| Tables S16 and S17 | `make_union_path_decomposition.R` | post-processing | `tables/union_path_decomposition.*`, `tables/alert_lead_time.md` |
| Table S18 | `make_ertb_design_actual_power_table.R` | 20260514; 20260546 (reversals) | `tables/ertb_design_actual_power.*` |
| Table S19 | `make_ertb_interruption_curve.R`, `make_pocock_matched_power.R` | 20260511; streams 20260542, 20260543; Pocock scenarios 20260544, 20260545 | `tables/ertb_interruption_*`, `tables/pocock_matched_power_summary.csv` |
| Table S20 | `make_premium_decomposition.R` | streams 20260542, 20260543; calibration 20260554 | `tables/premium_decomposition*.csv` |
| Figure S1, Table S21 | `make_arma_permutation_figure.R` | 20260506 | `figures/figS1_arma_permutation`, `tables/arma_adaptive_permutation_sensitivity.*` |
| Figure S2, Table S22 | `make_factt_permutation_figure.R` | 20260508 | `figures/figS2_factt_permutation`, `tables/factt_adaptive_permutation_sensitivity.*` |

Two outputs are not printed in the article. `make_pi_out_of_order_check.R`
(seed 20260555) backs the statement on which assignment probability prices a
wager when outcomes are confirmed out of enrollment order.
`make_ertb_interruption_curve.R` also redraws the comparison figure of the
original submission (`figures/legacy_fig2_interruption_curves`), which is not
part of the article.

## Kernels and the supplementary code

`supplementary_code/` is the code archive of the article as submitted: the
three kernels, its README, and the session record. The ESM prints the same
kernels in Sections 1.5, 2.4, and 3.5. `R/erte.R` and `R/ertc.R` are
byte-identical to their supplementary files. `R/ertb.R` differs in one place:
the supplementary file ends with a demonstration block that runs a full
simulation and writes four PDF files whenever the file is sourced in an
interactive R session. Here that block runs only after
`options(ert.run_demo = TRUE)`. `make check` verifies both statements.

## Trial data

The ARMA and FACTT analyses use deidentified patient-level data from BioLINCC.
The data are not redistributed here; see `data/README.md`. The shared datasets
do not include verified calendar enrollment order, so crossing indices are
positions in the dataset order.

## Tags in the code comments

Script headers keep the tags of the peer-review record. `R2.1` means Reviewer
2, comment 1 of the first round; an `m` marks a minor comment (`R4.m9`).
`T1` to `T10` are the task numbers of the revision's simulation plan, and
other short codes (`D1`, `A9`) are item numbers of the same plan. A "freeze
amendment" is a simulation added after the main simulation freeze. Each one
uses a fresh seed or replays the frozen streams, so earlier results did not
move.

## License

MIT for the source code and scripts; see `LICENSE`.
