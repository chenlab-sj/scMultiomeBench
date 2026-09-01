# `results/` — collected metric tables

Copies of the processed metric tables the figures are drawn from, one directory
per dataset key (matching `datasets/`), named for the manuscript figure each
table set feeds. FIGURES.md maps each figure to its generating script and inputs.

The merged final tables (`fig2b_matrix.csv`, `metrics_matrix.csv`,
`fig6_scores.csv`) are **single-file figure inputs**: each metrics-matrix
plotting script loads its merged table directly by default and regenerates the
published figure from that one file (verified to round-trip identically). Set
`REBUILD_MATRIX=1` to rebuild the merged table from the per-metric tables plus
the `published_reference/` splice instead.

| Directory | Feeds |
|---|---|
| `pbmc3k/fig2b/` | Fig 2B (23-method matrix; also read by Fig 3A, 6C, 7, S13) — includes the `knn_pred_label__<method>.csv` files Fig 2A reads (run `make_umap_grid.py` with `KNN_DIR` pointed here) |
| `pbmc3k/latents/` | per-method latent embeddings for Fig 2A (scglue_multiome, scVI, Seurat_CCA, scJoint, Conos) and Fig 3B–D (scglue_multiome, scJoint, scBridge, Cobolt) |
| `pbmc3k/fig4ab/` | Fig 4A/B, pbmc3k side |
| `pbmc3k/ksensitivity/` | Fig S2A |
| `pbmc3k/geneactivity_archr/` | Fig S3 |
| `pbmc3k/sum_metrics_clean.csv` | pbmc3k composite input read by Fig 7 / Fig S13 |
| `pbmc10k/` | Fig S2C |
| `pbmc_parse/` | Fig 6C, Fig S12 |
| `pbmc_parse/latents/`, `pbmc_parse/knn_pred_label/` | per-method latents + KNN-predicted labels for the Fig S12 per-method UMAP panels (`plot_umap.py`) |
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
| `<dataset>/published_reference/` | previously published baseline values and ground-truth labels for that dataset — **inputs, not outputs** of this repository, spliced into Fig 2B and Fig S6 and read by several other figures (README section 7). 33 files total; the two oversized RMS tables under `rms/published_reference/Mast607A/knn_test/` are gitignored, available on request |
