# `published_reference/` — published baselines and ground-truth labels

Inputs, not outputs, of this repository: 33 files — 30 previously published baseline metric tables
and 3 ground-truth `label.csv` tables (31 committed; the two oversized RMS tables noted below are
gitignored and available on request) — carved out of the original working directories so that the
figures which compare new results against previously published numbers can be rebuilt from this
repository alone. Scripts reach this directory through the `PUBLISHED_REF` configuration variable
(`config/config.sh`), which defaults to this location.

Layout is one directory per dataset key:

| Directory | Contents |
|---|---|
| `pbmc3k/` | published baseline tables for the Fig 2B splice (`metrics/`, `benchmark_matrix/`, `knn_test/`, `peak_similarity/`) |
| `pbmc10k/` | `knn_test/` k-sensitivity baseline and `peakdist_adj_random.csv` for Fig S2 |
| `bmmc_d1/` | `label.csv` ground truth, `benchmark_matrix/`, published metric tables for Fig 6A / Fig S11 |
| `brca/` | aggregated, de-identified tables only (`HT243B1-S1H4*/`) for Fig 4 / Fig S6 — no raw or per-cell BRCA data |
| `rms/Mast607A/` | `label.csv` ground truth, `benchmark_matrix/`, `knn_test/` for Fig 5 / Figs S8–S10 |

Only two plotting scripts actually splice these values into a figure:
`datasets/pbmc3k/04_figures/fig2b/plot_metrics_matrix.R` (Fig 2B) and
`datasets/brca/04_figures/plot_metrics_matrix.R` (Fig S6). Every other `SPLICE_PUBLISHED` branch in
the tree is dead code retained for provenance. See README section 7 and DATA.md for details.

Two oversized RMS tables are **not** in the public repository (gitignored;
`predicted_ataclabel_foxo1value.csv`, 88.7 MB, and `predicted_ataclabel_myod1value.csv`, 18.1 MB).
They are required by Fig 5B and Figs S8/S9 and are available from the authors on request.
