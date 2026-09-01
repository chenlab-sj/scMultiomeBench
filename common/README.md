# `common/` — vendored helpers shared across datasets

On the original cluster these files lived in one shared analysis directory that many scripts
referenced by absolute path; the path scrub maps that directory onto `common/`, reached through the
`BENCHMARK_FUN_DIR` configuration variable (`config/config.sh`, default `$PROJECT_ROOT/common`).
That is why cross-dataset files sit here rather than under any one `datasets/<dataset>/`.

## Code helpers

| File | Used by |
|---|---|
| `benchmark_fun.py` | imported (`import benchmark_fun`) by the `00_compute_metrics.py` scripts — the exact metric definitions (`lat2pdist`, `find_louvain_res`, KNN helpers) |
| `0_create10x_datafmt.R` | sourced by 11 preprocessing scripts (pbmc3k, pbmc10k-derived, pbmc_parse, bmmc_d1, brca) to write 10x-format h5 — the single canonical copy; a byte-identical `pbmc10k/` clone was removed 2026-09-01 and its consumers repointed |
| `peak_similarity/export_groupbwg.R` | sourced by every `*peak_similarity.R` — pseudobulk coverage-track export |
| `SJRHB013758_X2/plot_predpeak/plot_regionpeak_fun.R` | sourced by `datasets/rms/03_metrics/04_compute_new_pileups.R` |

## Data inputs

| File | Read by |
|---|---|
| `1_pbmc10k_annot0818.csv` | `datasets/pbmc10k/02_methods/scJoint/00_make_scjoint_h5.R` (cell-type list; unique — a different schema from `results/pbmc10k/label.csv`) |

The pbmc3k ground-truth labels formerly duplicated here now live only at
`results/pbmc3k/fig2b/label.csv`; the two consumers
(`datasets/pbmc_parse/03_metrics/00_build_label.py` and
`datasets/pbmc3k/03_metrics/ksensitivity/00_compute_ktest_legacy18.py`) point there.

`datasets/bmmc_d1/03_metrics/benchmark_metrics_lib.py` additionally reads `common/BMMC_d1/label.csv`;
that file is byte-identical to the shipped `results/bmmc_d1/published_reference/label.csv` — point the
placeholder there (or copy it in) when rerunning the BMMC metrics.

## Output locations (created at run time, intentionally not committed)

Two metric scripts write their results under `common/` because that is where the original shared
directory collected them: `common/BMMC_d1/benchmark_matrix/` (BMMC metric matrices, from
`benchmark_metrics_lib.py`) and `common/pbmc3k/knn_test/` (pbmc3k k-sensitivity outputs, from
`ksensitivity/00_compute_ktest_legacy18.py`). Both scripts `os.makedirs` their output directory, so
the empty state is fine; the published copies of these outputs ship under `results/` (including
each dataset's `published_reference/` subdirectory).
