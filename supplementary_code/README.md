# Supplementary code: e-RT kernels (R)

Code for "Sequential randomization tests using e-values for clinical trial
monitoring across endpoint types" (Zampieri, Albuquerque, Ramdas).

## Contents

- `SuppCode_eRTb.R` (binary endpoints, e-RTb): the design-wager calculation
  (`binary_design_lambdas`), the per-patient wealth update with the adaptive
  wager (`compute_eRT`), and the trial-simulation harness (`simulate_eRT`).
- `SuppCode_eRTe.R` (event-only streams, e-RTe): the event coin
  (`event_coin`), the per-event wealth update (`compute_eRTe`), and the
  event-stream simulation harness (`simulate_eRTe`).
- `SuppCode_eRTc.R` (continuous endpoints, e-RTc): the design wager
  (`continuous_design_lambda`), the running effect-size estimate
  (`cohens_d_estimate`), and the per-patient update (`compute_eRTc`).
- `sessionInfo.txt`: R and package versions.

Each kernel is also printed in full in the Electronic Supplementary Material
(ESM Sections 1.5, 2.4, and 3.5), with implementation notes. Simulation seeds
are fixed and documented per result in the ESM (the "Code:" line at the end of
each subsection). The kernels implement simple randomization at a fixed
allocation probability; the conditional-probability rule for restricted
randomization is stated in the main text and demonstrated in ESM Section 1.3.

## Requirements

R with the tidyverse. Exact versions: `sessionInfo.txt`.
