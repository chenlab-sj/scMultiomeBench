# `common/` — shared helpers

Helper code shared across datasets, reached through the `BENCHMARK_FUN_DIR` configuration variable
(default `$PROJECT_ROOT/common`).

| File | Used by |
|---|---|
| `benchmark_fun.py` | imported by the `00_compute_metrics.py` scripts — the exact metric definitions |
| `0_create10x_datafmt.R` | sourced by the preprocessing scripts to write 10x-format h5 |
| `export_groupbwg.R` | sourced by every `*peak_similarity.R` — pseudobulk coverage-track export |
| `plot_regionpeak_fun.R` | sourced by `datasets/rms/03_metrics/04_compute_new_pileups.R` |
| `1_pbmc10k_annot0818.csv` | cell-type list for the pbmc10k scJoint input prep (covers cells absent from `results/pbmc10k/label.csv`) |

Two metric scripts create output directories under `common/` at run time (`BMMC_d1/`,
`pbmc3k/knn_test/`); their absence in the shipped tree is expected.
