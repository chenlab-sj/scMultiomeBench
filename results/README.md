# `results/` — figure inputs

The processed tables the manuscript's figures and tables are drawn from, one directory per dataset
key. **Every metrics table for all six datasets is in this repository**, together with the 23
per-method integrated embeddings behind the pbmc3k Fig 2B matrix.

The merged final tables (`fig2b_matrix.csv`, `metrics_matrix.csv`, `fig6_scores.csv`) are
single-file figure inputs: each metrics-matrix script loads its merged table directly and
regenerates the published figure. `FIGURES.md` maps each figure to its generating script and
inputs.

## Integrated embeddings — `pbmc3k/latents/`

The 23 per-method latent representations that the pbmc3k Fig 2B metrics were computed from, one
CSV per method (rows = cells, columns = latent dimensions; the first column is the cell barcode).
Cell counts differ between methods because some emit only the unpaired test cells and others also
carry the paired training cells — this is the input each method actually produced, unmodified.

`SOURCES.tsv` records, for every method: the staged filename, cell and dimension counts, the md5
of the file as shipped, and the original path it was produced at, written with the same
`${TOOLS_ROOT}` / `${PROJECT_ROOT}` placeholders used by
`datasets/pbmc3k/03_metrics/00_compute_metrics_sub.sh`, the script that consumed them. Every file
was md5-verified against that original.

## What feeds what

| Path | Feeds |
|---|---|
| `pbmc3k/fig2b/` | Fig 2B matrix (`fig2b_matrix.csv`); Fig 3A (`celltype_metrics.csv`, `adj_atac_predaccu.csv`); Fig S3 (`sum_metrics_clean.csv`) |
| `pbmc3k/latents/` | Fig 2A + Fig 3B–D per-method embeddings (23 methods) |
| `pbmc3k/fig4ab/` | Fig 4A/B pbmc3k side (`nmi_df.csv`, `sd_df.csv`) |
| `pbmc3k/ksensitivity/` | Fig S2A/B |
| `pbmc3k/geneactivity_archr/compare_geneactivity.csv` | Fig S3 comparison record |
| `pbmc3k/sum_metrics_clean.csv` | Fig 7 / Fig S13 pbmc3k input |
| `pbmc10k/` | Fig S2C (`fig2b_matrix.csv`); Fig S2A/B (`ksensitivity/`) |
| `pbmc_parse/` | Fig S12 + Fig 6C (`metrics_matrix.csv`) |
| `bmmc_d1/` | Table S4 (`sum_metrics.csv`, `celltype_metrics.csv`, `adj_atac_predaccu.csv`, `peakdist.csv` → `tableS4_batch_correction.csv`); Fig 7 (`sum_metrics_clean.csv`) |
| `bmmc_d1/metrics/`, `bmmc_d1/crosssite/`, `bmmc_d1/s1d1_paired/` | Fig S11 matrices; Fig 6B (crosssite + s1d1 `metrics_matrix.csv`); s1d1 `peakdist.csv` feeds the Fig S11 9-metric table |
| `bmmc_d1/fig6/fig6_scores.csv` | Fig 6A |
| `brca/` | Fig 4D (`fig4d/`); Fig 4A/B BRCA side (`fig4ab/`); Fig S6 (`figS6/fig2b_matrix.csv`); Fig 7 (`sum_metrics_clean.csv`, `subset_stability.csv`, `reproducbility.csv`) |
| `rms/` | Fig 5A (`fig5a/metrics_matrix.csv` + `knn_pred_accu.csv`); Fig S10 (`figS10/`); Fig 7 (`sum_metrics_clean.csv`) |
| `summary/` | Fig 7 (`runtime_memory_sel_plus.csv`, `new_methods_grouped.csv`); Fig S13 (`weighting_sensitivity_ranks.csv`) |

The tables listed above are what the published figures are drawn from: every metrics-matrix figure
renders directly from the merged table in its dataset directory, and the pbmc3k embeddings carry
the integration results behind Fig 2A and Fig 3B–D.
