# `results/` — figure inputs

All tables the figures are drawn from, one directory per dataset key (matching
`datasets/`), named for the manuscript figure each table set feeds. The merged
final tables (`fig2b_matrix.csv`, `metrics_matrix.csv`, `fig6_scores.csv`) are
single-file figure inputs: each metrics-matrix script loads its merged table
directly by default and regenerates the published figure from that one file.
FIGURES.md maps each figure to its generating script and inputs.

| Directory | Feeds |
|---|---|
| `pbmc3k/fig2b/` | Fig 2B (23-method matrix; also read by Fig 3A, 6C, 7, S13) and the `knn_pred_label__<method>.csv` files Fig 2A reads |
| `pbmc3k/latents/` | per-method latent embeddings for Fig 2A and Fig 3B–D |
| `pbmc3k/fig4ab/` | Fig 4A/B, pbmc3k side |
| `pbmc3k/ksensitivity/` | Fig S2A |
| `pbmc3k/geneactivity_archr/` | Fig S3 |
| `pbmc3k/sum_metrics_clean.csv` | pbmc3k composite input read by Fig 7 / Fig S13 |
| `pbmc10k/` | Fig S2C |
| `pbmc_parse/` | Fig 6C, Fig S12 (incl. `latents/` + `knn_pred_label/` for the per-method UMAPs) |
| `bmmc_d1/` | Fig 6A metric tables; Table S4 |
| `bmmc_d1/metrics/` | Fig 6A, Fig S11 (3-site matrix) |
| `bmmc_d1/crosssite/`, `bmmc_d1/s1d1_paired/` | Fig 6B, Fig S11 |
| `bmmc_d1/fig6/` | Fig 6 assembly |
| `brca/fig4ab/` | Fig 4A/B, BRCA side |
| `brca/figS4_S5/` | Fig S4/S5 confusion grids |
| `brca/figS6/` | Fig S6 (BRCA metrics matrix) |
| `rms/fig5a/` | Fig 5A (RMS metrics matrix) |
| `rms/figS10/` | Fig S10 (RMS reproducibility) |
| `summary/` | Fig 7, Fig S13 |
| `brca/label.csv`, `brca/fig4d/`, `brca/reproducbility.csv`, `brca/subset_stability.csv`, `brca/sum_metrics_clean.csv` | Fig 4C truth labels; Fig 4D rep-1 accuracy series; Fig 7 BRCA inputs |
| `bmmc_d1/label.csv`, `bmmc_d1/sum_metrics_clean.csv` | BMMC ground truth (48k cells); Fig 6A / Fig 7 input |
| `rms/label.csv`, `rms/knn_test/`, `rms/sum_metrics_clean.csv` | RMS ground truth + KNN prediction tables (Figs 5B/5C/S8/S9; the two oversized `predicted_ataclabel_*value.csv` under `rms/knn_test/` are gitignored, available on request); Fig 7 input |
