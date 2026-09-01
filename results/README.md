# `results/` — collected metric tables

Copies of the processed metric tables the figures are drawn from, one directory
per dataset key (matching `datasets/`), named for the manuscript figure each
table set feeds. The plotting scripts read their own `03_metrics/` working
directories; these copies exist so every figure input ships with the repo.
FIGURES.md maps each figure to its generating script and inputs.

| Directory | Feeds |
|---|---|
| `pbmc3k/fig2b/` | Fig 2B (23-method matrix; also read by Fig 3A, 6C, 7, S13) |
| `pbmc3k/fig4ab/` | Fig 4A/B, pbmc3k side |
| `pbmc3k/ksensitivity/` | Fig S2A |
| `pbmc3k/geneactivity_archr/` | Fig S3 |
| `pbmc3k/sum_metrics_clean.csv` | pbmc3k composite input read by Fig 7 / Fig S13 |
| `pbmc10k/` | Fig S2C |
| `pbmc_parse/` | Fig 6C, Fig S12 |
| `bmmc_d1/` (top-level tables, incl. `tableS4_batch_correction.csv`) | Fig 6A metric tables; Table S4 |
| `bmmc_d1/metrics/` | Fig 6A, Fig S11 (3-site matrix) |
| `bmmc_d1/crosssite/`, `bmmc_d1/s1d1_paired/` | Fig 6B, Fig S11 |
| `bmmc_d1/fig6/` | Fig 6 assembly |
| `brca/fig4ab/` | Fig 4A/B, BRCA side |
| `brca/figS4_S5/` | Fig S4/S5 confusion grids |
| `brca/figS6/` | Fig S6 (BRCA metrics matrix) |
| `rms/fig5a/` | Fig 5A (RMS metrics matrix) |
| `rms/figS10/` | Fig S10 (RMS reproducibility) |
| `summary/` | Fig 7, Fig S13 |
