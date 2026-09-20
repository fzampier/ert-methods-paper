# Trial data (not in this repository)

The completed-trial illustrations (ESM Section 6: Figures S1 and S2, Tables S21
and S22, and the ARMA and FACTT numbers in the main text) use deidentified
patient-level data from the ARMA and FACTT trials. The data come from the NHLBI
Biologic Specimen and Data Repository Information Coordinating Center
(BioLINCC, https://biolincc.nhlbi.nih.gov). They are available from BioLINCC
subject to BioLINCC approval and a data-use agreement, and they cannot be
redistributed here. Nothing under `data/` except this file is tracked by git.

## Where the scripts look

```text
data/derived/biolincc_validated/csv_for_analysis/arma.csv
data/derived/biolincc_validated/csv_for_analysis/factt.csv
data/derived/biolincc_validated/results/arma_ert_summary_R.json    (number check only)
data/derived/biolincc_validated/results/factt_ert_summary_R.json   (number check only)
```

To keep the files outside the repository, set `ERT_BIOLINCC_DIR` to the folder
that holds `csv_for_analysis/` and `results/`.

## Columns the scripts read

`arma.csv`, one row per randomized patient (861 rows):

| Column | Meaning |
| --- | --- |
| `enroll_order` | Order of the patient in the BioLINCC dataset (1 to 861) |
| `arm` | 0 = 12 ml/kg tidal volume (control), 1 = 6 ml/kg (intervention) |
| `death` | 0 = alive, 1 = dead (BioLINCC `STATUS = 2`) |
| `unassist_day` | Day of first unassisted breathing (missing when not recorded) |

`factt.csv`, one row per randomized patient (1,000 rows):

| Column | Meaning |
| --- | --- |
| `enroll_order` | Order of the patient in the BioLINCC dataset (1 to 1,000) |
| `arm_fluid` | 0 = liberal (control), 1 = conservative (intervention) |
| `arm_catheter` | 0 = central venous catheter (control), 1 = pulmonary artery catheter |
| `death_d60` | Death by day 60 (0/1) |
| `vfd28` | Ventilator-free days to day 28 (integer, 0 to 28) |

## Order of patients

The BioLINCC datasets are deidentified and do not include verified calendar
enrollment order. `enroll_order` is the order of the shared dataset. Crossing
indices in the article are therefore positions in that order, not prospective
monitoring times, and the permutation analyses show how much each result
depends on the order.

## What the trial scripts write

`scripts/make_arma_permutation_figure.R` and
`scripts/make_factt_permutation_figure.R` write aggregate results only: the two
figures, the two permutation tables under `tables/`, and working files under
`outputs/` (ignored by git). No shuffled patient-level dataset is saved.
