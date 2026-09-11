# `results/` — figure inputs

**This directory is not published in this repository** (kept out of version control by policy) —
it is **available directly from the authors on request**. This file documents what it contains for
anyone who has obtained a copy.

All tables the manuscript figures and tables are drawn from, one directory per dataset key. The
merged final tables (`fig2b_matrix.csv`, `metrics_matrix.csv`, `fig6_scores.csv`) are single-file
figure inputs: each metrics-matrix script loads its merged table directly and regenerates the
published figure. FIGURES.md maps each figure to its generating script and inputs.

| Path | Feeds |
|---|---|
| `pbmc3k/fig2b/` | Fig 2B matrix (`fig2b_matrix.csv`); Fig 2A (`label.csv`, `knn_pred_label__*.csv`); Fig 3A (`celltype_metrics.csv`, `adj_atac_predaccu.csv`); Fig S3 (`sum_metrics_clean.csv`) |
| `pbmc3k/latents/` | Fig 2A + Fig 3B–D per-method embeddings |
| `pbmc3k/fig4ab/` | Fig 4A/B pbmc3k side (`nmi_df.csv`, `sd_df.csv`) |
| `pbmc3k/ksensitivity/` | Fig S2A/B |
| `pbmc3k/geneactivity_archr/compare_geneactivity.csv` | Fig S3 comparison record |
| `pbmc3k/sum_metrics_clean.csv` | Fig 7 / Fig S13 pbmc3k input |
| `pbmc10k/` | Fig S2C (`fig2b_matrix.csv`); Fig S2A/B (`ksensitivity/`); ground truth (`label.csv`) |
| `pbmc_parse/` | Fig S12 + Fig 6C (`metrics_matrix.csv`, `label.csv`, `latents/`, `knn_pred_label/`) |
| `bmmc_d1/` | Table S4 (`sum_metrics.csv`, `celltype_metrics.csv`, `adj_atac_predaccu.csv`, `peakdist.csv` → `tableS4_batch_correction.csv`); Fig 7 (`sum_metrics_clean.csv`); ground truth (`label.csv`) |
| `bmmc_d1/metrics/`, `bmmc_d1/crosssite/`, `bmmc_d1/s1d1_paired/` | Fig S11 matrices; Fig 6B (crosssite + s1d1 `metrics_matrix.csv`); s1d1 `peakdist.csv` feeds the Fig S11 9-metric table |
| `bmmc_d1/fig6/fig6_scores.csv` | Fig 6A |
| `brca/` | Fig 4C (`label.csv`); Fig 4D (`fig4d/`); Fig 4A/B BRCA side (`fig4ab/`); Fig S6 (`figS6/fig2b_matrix.csv`); Fig 7 (`sum_metrics_clean.csv`, `subset_stability.csv`, `reproducbility.csv`) |
| `rms/` | Fig 5A (`fig5a/metrics_matrix.csv` + `knn_pred_accu.csv`, `label.csv`); Figs 5B/5C/S8/S9 (`knn_test/`, `label.csv` — the two oversized `predicted_ataclabel_*value.csv` are gitignored, available on request); Fig S10 (`figS10/`); Fig 7 (`sum_metrics_clean.csv`) |
| `summary/` | Fig 7 (`runtime_memory_sel_plus.csv`, `new_methods_grouped.csv`); Fig S13 (`weighting_sensitivity_ranks.csv`) |
