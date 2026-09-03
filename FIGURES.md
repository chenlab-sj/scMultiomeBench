# FIGURES.md — figure-to-script index

Every manuscript and supplementary figure, the script in this repository that draws it, that
script's inputs, and the dataset it comes from.

Paths are relative to the repository root. `<metrics>` denotes the metric working directory for
that dataset (`datasets/<dataset>/03_metrics/`, or the sub-experiment's own `03_metrics/`), which
is where the `00_*`–`05_*` metric scripts write `sum_metrics.csv`, `celltype_metrics.csv`,
`adj_atac_predaccu.csv` and `peakdist.csv`. Those intermediate CSVs are produced by the metric
step; they are not tracked in `results/`.

Every metrics-matrix figure regenerates from its shipped merged table in `results/` by default
(single-file mode). Rebuilding a merged table from scratch (`REBUILD_MATRIX=1`) additionally needs
the reference tables of the original run, which are not shipped (available from the authors on
request); `SPLICE_PUBLISHED` branches other than pbmc3k Fig 2B and BRCA Fig S6 are dead code kept
for provenance.

Datasets: `pbmc3k`, `pbmc10k`, `pbmc_parse` (Parse Evercode), `bmmc_d1`, `brca` (HT243B1-S1H4),
`rms` (Mast607A).

---

## Main figures

### Fig 1 — study design and benchmarking workflow

Schematic. **No generating script.** The file census contains no entry for Fig 1; it was drawn by
hand, not plotted.

### Fig 2 — pbmc3k benchmark overview

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| 2A | `datasets/pbmc3k/04_figures/fig2a/make_umap_grid.py` | per-method latent CSVs for scglue(multiome), scVI, Seurat(CCA), scJoint, Conos — **shipped in `results/pbmc3k/latents/`**; `label.csv` and `knn_pred_label__<method>.csv` — shipped in `results/pbmc3k/fig2b/` (set `KNN_DIR` there) | pbmc3k |
| 2B | `datasets/pbmc3k/04_figures/fig2b/plot_metrics_matrix.R` | `<metrics>/sum_metrics.csv`, `celltype_metrics.csv`, `adj_atac_predaccu.csv`, `peakdist.csv`; rebuild mode (`REBUILD_MATRIX=1`) additionally needs the reference tables — not shipped, available from the authors on request | pbmc3k |

Fig 2B is the 23-method × 9-metric matrix and the composite ranking. It writes
`fig2b_matrix.csv` and `sum_metrics_clean.csv`, which several later figures (Fig 3A, Fig 6C,
Fig 7, Fig S13) read.

Upstream metric steps for both panels: `datasets/pbmc3k/03_metrics/00_compute_metrics.py` →
`01_adjust_accuracy.py` → `02_peak_similarity.R`.

Supporting, not a published panel: `datasets/pbmc3k/04_figures/fig2b/plot_umap.py` renders a
single-method UMAP panel (`--latent`, `--label`, `--knn-pred`) for the methods added in revision.

### Fig 3 — per-cell-type behaviour, pbmc3k

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| 3A | `datasets/pbmc3k/04_figures/fig3a/plot_fig3a.R` | `<fig2b dir>/celltype_metrics.csv`, `adj_atac_predaccu.csv`, `label.csv`, `fig2b_matrix.csv` (supplies the method order, restricted to the 14 Fig 3A methods) | pbmc3k |
| 3B–3D | `datasets/pbmc3k/04_figures/fig3bcd/make_celltype_dist.py` | latent CSVs for scglue(multiome), scJoint, scBridge, Cobolt — **shipped in `results/pbmc3k/latents/`**; `label.csv` | pbmc3k |

3B–3D are the CD16 monocyte, naive CD8 T and memory T inter-omics distance panels for the four
curated methods.

### Fig 4 — reproducibility (BRCA + pbmc3k)

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| 4A (pairwise NMI, pbmc3k \| BRCA) | `datasets/brca/04_figures/plot_fig4_combined_brca_pbmc.py` | `nmi_df_7type.csv`, `sd_df_7type.csv` (BRCA, from `03_metrics/reproducibility/`); `nmi_df.csv`, `sd_df.csv` (pbmc3k, from `datasets/pbmc3k/03_metrics/reproducibility/00_reproduce_curated.py`) | brca + pbmc3k |
| 4B (per-cell-type SD, pbmc3k \| BRCA) | same script — it draws all four subpanels in one row | same | brca + pbmc3k |
| 4C (confusion matrices for Portal / Seurat(CCA) / scVI + cell-count bar) | `datasets/brca/04_figures/plot_fig4c_confusion.py` (ported from the original analysis notebook `macro_sub-Fig4C1.ipynb`; the published `Fig4C.pdf` is md5-identical to that notebook's output) | `rep_knn_k10_pred_label_7type.csv` (from `03_metrics/reproducibility/03_fig4_reproduce_7type.ipynb`) + `results/brca/label.csv` | brca |
| 4D (macrophage subsampling) | `datasets/brca/05_macrophage_subsample/04_figures/plot_fig4d_macrophage.py` | `results/brca/fig4d/knn_k10sub{0..5}_pred_accu.csv` (rep 1); `knn_k10sub{1..5}_pred_accu.csv` for reps 2–5 from `05_macrophage_subsample`; `new_methods_macro_accu.csv` from `05_macrophage_subsample/03_metrics/compute_new_methods_accu.py` | brca |

Upstream for 4A/4B on the BRCA side: `datasets/brca/03_metrics/reproducibility/00_cobolt_nmi.py`,
`01_midas_pick.py`, `02_reproduce_metrics.py`, `03_fig4_reproduce_7type.ipynb` (the notebook writes
`nmi_df_7type.csv`, `sd_df_7type.csv` and `rep_knn_k10_pred_label_7type.csv`), plus
`results/brca/label.csv` and `reproducbility.csv` (the legacy-method Louvain clusters the
notebook reuses, `louvain_cluster_reproduce.csv`, are not shipped — available on request).
Upstream for 4D: `datasets/brca/05_macrophage_subsample/01_preprocess/macro_subset_rep{1..5}.R`
then the MaxFuse / MIDAS / scButterfly runs in that sub-experiment.

Superseded, not included in this repository: `plot_fig4.py` (BRCA and pbmc3k preview plotters),
`plot_fig4_manuscript_7type.py` (BRCA-only version of 4A/4B), `make_repro_fig.py` (slope charts, a
different chart type that appears nowhere), and `macrophage_sub_fig4d.ipynb` (the scratch notebook
behind 4D). They exist only in the pre-revision working directories.

### Fig 5 — RMS (Mast607A)

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| 5A (metrics matrix + ranking) | `datasets/rms/04_figures/fig5a/plot_metrics_matrix.R` | `<metrics>/sum_metrics.csv`, `celltype_metrics.csv`, `adj_atac_predaccu.csv`, `peakdist.csv`; rebuild mode (`REBUILD_MATRIX=1`) additionally needs the reference tables — not shipped, available from the authors on request | rms |
| 5B (Euclidean distance of predicted vs true coverage tracks) | `datasets/rms/04_figures/fig5b/plot_fig5b.py` | `results/rms/knn_test/predicted_ataclabel_{myod1,foxo1}value.csv` (published methods) + `new_methods_{myod1,foxo1}value.csv` from `datasets/rms/03_metrics/04_compute_new_pileups.R` (new methods) | rms |
| 5C (MYOD1 / FOXO1 pileup tracks, 7 rows: Annotation + 5 methods + Random) | `datasets/rms/04_figures/figS7_S9/plot_pileup_grid_2page.R` with the Fig 5C method selection — see note below | same inputs as Figs S8/S9 | rms |

Upstream for 5A: `datasets/rms/03_metrics/00_prep_latents.py` → `01_compute_metrics.py` →
`02_adjust_accuracy.py` → `03_peak_similarity.R`.

**Fig 5C method selection.** Fig 5C is the same pileup pipeline and the same inputs as Figs S8/S9;
only the method selection differs. To reproduce it, replace the `PAGES` list (lines 30–32 of
`plot_pileup_grid_2page.R`) with the single entry
`fig5c = c("scglue \n (multiome)", "Seurat \n (CCA)", "scVI", "scBridge", "Portal")` — the
Annotation (top) and Random (bottom) rows are added automatically by `build_region()`, giving the
published 7-row page. The original scripts carry this exact 7-row list as a superseded assignment
(`plot_predpeak_final.R` line 185, `plot_predpeak_final_foxo1.R` line 180, in the pre-revision
working directory), which is the fossil of the state that drew the published panel.

### Fig 6 — generalization to independently sequenced data

`datasets/bmmc_d1/04_figures/fig6_combined.py` draws all three panels as one figure.
`datasets/bmmc_d1/04_figures/fig6_panels.py` re-exports A/B/C as separate vector PDFs and reuses
the combined script.

| Panel | Content | Key inputs | Dataset |
|---|---|---|---|
| 6A | BMMC benchmark, four metric facets | `<metrics>/sum_metrics.csv`, `celltype_metrics.csv`, `adj_atac_predaccu.csv`, `peakdist.csv`; `results/bmmc_d1/sum_metrics_clean.csv`; `results/bmmc_d1/label.csv` | bmmc_d1 |
| 6B | BMMC same-donor cross-site (s1d1 RNA + s4d1 ATAC vs s1d1 multiome) | matrices from `datasets/bmmc_d1/05_crosssite/03_metrics/` and `datasets/bmmc_d1/06_s1d1_paired/03_metrics/`, assembled by the `plot_metrics_matrix.R` in each of those sub-experiments' `04_figures/` | bmmc_d1 |
| 6C | PBMC cross-platform (Parse RNA + 10x ATAC vs 10x multiome) | `datasets/pbmc_parse/03_metrics/` matrix + pbmc3k `fig2b_matrix.csv` | pbmc_parse + pbmc3k |

Panel C also exists as a standalone plotter,
`datasets/pbmc_parse/04_figures/plot_generalization_scatter.py`, which reads
`pbmc_parse` `metrics_matrix.csv` and pbmc3k `fig2b_matrix.csv` directly.

Both scatters recompute the 6-metric unpaired composite so the multiome baseline is not credited
for the same-cell metrics (ks.statistic / ARI / AMI) it alone can have.

Upstream for 6A: `datasets/bmmc_d1/03_metrics/00_fix_bmmc_barcodes.py` → `01_prep_latents.py` →
`02_run_bmmc_metrics.py` → `03_knn_bmmc.py` → `04_adjust_accuracy.py` → `05_peak_similarity.R`,
with `benchmark_metrics_lib.py` as the shared library. The six `(batch)` method variants
(BindSC_batch, MIDAS_batch, scJoint_batch, scVI_batch, scglue_batch, scglue_multiome_batch) enter
here.

**Dependency note (resolved):** `fig6_combined.py` imports `bmmc_crosssite_data`, `pbmc_data` and
`plot_pbmc` from `make_r1_brca_pbmc_composite.py`, which ships alongside it in
`datasets/bmmc_d1/04_figures/`, so the import resolves from a clean checkout. That module also
contains reviewer-response plotting code for the excluded HT137 / cross-donor matrices; only the
three imported functions are used here.

### Fig 7 — cross-scenario summary

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| Fig 7 | `summary/04_figures/plot_fig7_summary.R` | `sum_metrics_clean.csv` for pbmc3k, Mast607A (rms), HT243B1-S1H4 (brca) and BMMC_d1; `HT243B1-S1H4/subset_stability.csv`; `HT243B1-S1H4/reproducbility.csv`; `runtime_memory_sel_plus.csv`; `new_methods_grouped.csv` | all four benchmark datasets |

Driver: `summary/04_figures/plot_fig7_summary_sub.sh`. The 9 originally published methods keep
their original values; the 3 methods added in revision are appended from
`new_methods_grouped.csv`, and their runtime/memory come from the LSF logs via
`runtime_memory_sel_plus.csv`.

---

## Supplementary figures

### Fig S1 — validation of cell-type labels

`datasets/pbmc3k/04_figures/figS1/make_figS1_label_validation.R` assembles the figure by sourcing
the three panel scripts, so the panels can never drift from their standalone outputs.

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S1A | `datasets/pbmc3k/04_figures/figS1/marker_plot_label.R` | pbmc3k h5 + `label.csv` (10x ground-truth labels); `LABELS=singleR` switches it to `label_singleR.csv` | pbmc3k |
| S1B | `datasets/rms/04_figures/label_validation/plot_marker_label_validation.R` | RMS h5 + Mast607A sub-stage `label.csv` | rms |
| S1C | `datasets/pbmc3k/03_metrics/label_sensitivity/02_remap_labels.R` | 10x labels vs SingleR labels (`00_singleR_annotate.R`, `01_singleR_confidence.R`, `fine_label_map.R`); exposes the confusion plot object | pbmc3k |

Supporting, not a panel:
`datasets/pbmc3k/03_metrics/label_sensitivity/rank_stability/rank_stability.py` computes the
Spearman correlation of the composite ranking under 10x vs SingleR labels. That number lives in
the response text; it is not drawn in Fig S1.

### Fig S2 — k-sensitivity and the pbmc10k benchmark

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S2A + S2B (two-row k-sensitivity, one shared legend and x-axis) | `summary/04_figures/plot_figS2a_ksens_combined.py` | pbmc3k `kNN_celltype_accu_sum.csv` + `new_methods_ktest_long.csv` (from `datasets/pbmc3k/03_metrics/ksensitivity/00_compute_ktest_legacy18.py` and `01_compute_ktest_new_methods.py`); pbmc10k `results/pbmc10k/ksensitivity/kNN_celltype_accu_sum.csv` + `new_methods_ktest_long.csv` (same dir) | pbmc3k + pbmc10k |
| S2C (pbmc10k metrics matrix) | `datasets/pbmc10k/04_figures/plot_metrics_matrix.R` | `<metrics>/sum_metrics.csv`, `celltype_metrics.csv`, `adj_atac_predaccu.csv`, `peakdist.csv`; rebuild mode (`REBUILD_MATRIX=1`) additionally needs the reference tables — not shipped, available from the authors on request | pbmc10k |

k = 5, 10, 20, 40, 80 for 23 methods per dataset. The published benchmark uses k = 10, cosine
distance, distance-weighted.

Upstream for S2C: `datasets/pbmc10k/01_preprocess/00_data_prep.R`, then
`03_metrics/00_compute_metrics.py` → `01_adjust_accuracy.py` → `02_peak_similarity.R` →
`03_compute_ktest_new_methods.py` → `04_compute_ktest_reverse.py`.

Superseded, not included in this repository: `plot_ksens.py` (pbmc3k only), `plot_figS1a.py`
(pbmc10k only) and `plot_figS1b.py` (the reverse RNA-label direction, which is not in Fig S2 at
all). "figS1" was the local working name for this figure in the pre-revision directories.

### Fig S3 — Signac vs ArchR gene activity

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S3 | `datasets/pbmc3k/05_geneactivity_archr/04_figures/compare_geneactivity.py` | `datasets/pbmc3k/04_figures/fig2b/sum_metrics_clean.csv` (Signac) and the ArchR-side `metrics/sum_metrics_clean.csv` | pbmc3k |

The ArchR gene-activity matrices come from
`05_geneactivity_archr/01_archr/archr_genescores_pbmc3k.R` and `archr_GeneActivity.R`; the methods
are re-run against them under `05_geneactivity_archr/02_methods/`, latents normalised by
`03_metrics/00_fix_latent_format.py`, then scored with the same metric scripts. The ArchR-side
metrics matrix is built by re-running `datasets/pbmc3k/04_figures/fig2b/plot_metrics_matrix.R`
against the ArchR metric directory — that is what
`05_geneactivity_archr/04_figures/plot_metrics_matrix_sub.sh` does; there is no second copy of the
R script.

Note for readers: scJoint is excluded from the Fig S3 headline comparison
(`in_headline = False`, 0.676 → 0.101).

### Fig S4 and Fig S5 — BRCA triplicate reproducibility confusion grids

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S4 (methods ranked 1–7) | `datasets/brca/04_figures/plot_figS4_S5_confusion.py` | `rep_knn_k10_pred_label_7type.csv` (columns `<method>-1/-2/-3`) + `results/brca/label.csv` | brca |
| S5 (methods ranked 8–14) | same script, second page | same | brca |

One script writes both pages; the split is the `FIGS2` / `FIGS3A` lists inside it (again, local
working names). Method order follows the BRCA composite-metrics ranking so the reproducibility
grid reads in the same order as Fig S6.

### Fig S6 — BRCA metrics matrix

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S6 | `datasets/brca/04_figures/plot_metrics_matrix.R` | `<metrics>/major/{sum_metrics,celltype_metrics,adj_atac_predaccu}.csv` and `<metrics>/peakdist.csv`; rebuild mode (`REBUILD_MATRIX=1`) additionally needs the reference tables — not shipped, available from the authors on request | brca |

In rebuild mode its `_sub.sh` driver sets `SPLICE_PUBLISHED=1` (the reference tables above are
only read then). Upstream: `datasets/brca/01_preprocess/00_data_prep.R` and
`01_prep_HT263_train_commonpeaks.R`, then `03_metrics/00_compute_metrics.py` →
`02_adjust_accuracy.py` → `03_peak_similarity.R`.

### Fig S7 — MYOD1 / FOXO1 / MEOX2 multiome vs annotation tracks

Recovered from the original pre-revision analysis directory (2026-08-31). One script per locus;
each draws the two-row (Multiome vs Annotation) pileup strip for its locus, faceted over the three
myogenic cell types, and the three strips were assembled into the published figure manually. The
scATAC side is the RMS sample **SJRHB013758_X2**; the Multiome row comes from the paired Mast607A
data.

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S7, MYOD1 strip | `datasets/rms/04_figures/figS7_multiome_vs_annotation/plot_predpeak_final_v2.R` | `SJRHB013758_X2/knn_test/{knn_k10_pred_label,knn_k10_pred_accu,predicted_ataclabel_myod1value}.csv`; `SJRHB013758_X2_clusters.csv`; Mast607A `lca_label.csv` + `predicted_ataclabel_myod1value.csv` | rms |
| S7, FOXO1 strip | `.../plot_predpeak_final_foxo1_v2.R` | same set with the `foxo1` value tables | rms |
| S7, MEOX2 strip | `.../plot_predpeak_final_meox2_v2.R` | same set with the `meox2` value tables | rms |

The per-cell `predicted_ataclabel_*value.csv` inputs are not committed (size) and are **available
from the authors on request**, like the Fig S8/S9 inputs. The two pileup scripts in
`datasets/rms/04_figures/figS7_S9/` are named for the S7–S9 block but generate S8 and S9 only.

### Fig S8 and Fig S9 — RMS pileup grids (14 methods, two pages)

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S8 (page 1) | `datasets/rms/04_figures/figS7_S9/plot_pileup_grid_2page.R` | see below | rms |
| S9 (page 2) | same script, second page | see below | rms |

Inputs:
`results/rms/knn_test/predicted_ataclabel_myod1value.csv`,
`predicted_ataclabel_foxo1value.csv`, `knn_k10_pred_label.csv`, `knn_k10_pred_accu.csv`,
`results/rms/label.csv`, plus the new methods'
`new_methods_{myod1,foxo1}value.csv`, `knn_pred_accu.csv` and `knn_pred_label__<method>.csv` from
`datasets/rms/03_metrics/04_compute_new_pileups.R`. The script caches the assembled tracks in
`region_long.rds` and reuses it on later runs.

> **Two inputs are gitignored and are not in the public repository:**
> `results/rms/knn_test/predicted_ataclabel_foxo1value.csv` (88.7 MB) and
> `predicted_ataclabel_myod1value.csv` (18.1 MB). They are required for Fig 5B, Fig S8 and Fig S9.
> **Available from the authors on request.**

Superseded, not included in this repository: `plot_figS6.R`, the single-page 14-row version. It
was too tall to publish and was split into the two-page script three days later; the two-page
script reuses its label and track recipe.

### Fig S10 — RMS reproducibility

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S10 (NMI boxplot + per-cell-type SD, combined) | `datasets/rms/04_figures/figS10/plot_figS10_combined.py` | `nmi_df.csv`, `sd_df.csv` from `datasets/rms/03_metrics/reproducibility/00_pick3_reps.py` → `01_select_pick3_triplets.py` → `02_reproduce_metrics.py` | rms |

Superseded, not included in this repository: `plot_fig4_manuscript.py` (emits the two panels as
separate PDFs; merged six minutes later into the published single figure) and `plot_fig4.py`
(draft preview plotter). Both exist only in the pre-revision working directories.

### Fig S11 — BMMC supplementary matrices

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| BMMC 3-site metrics matrix | `datasets/bmmc_d1/04_figures/plot_metrics_matrix.R` | `<metrics>/{sum_metrics,celltype_metrics,adj_atac_predaccu,peakdist}.csv`; rebuild mode (`REBUILD_MATRIX=1`) additionally needs the reference tables — not shipped, available from the authors on request | bmmc_d1 |
| cross-site matrix | `datasets/bmmc_d1/05_crosssite/04_figures/plot_metrics_matrix.R` | `05_crosssite/03_metrics/` outputs (`00_build_label.py` → `01_prep_latents.py` → `02_compute_metrics.py` → `04_adjust_accuracy.py` → `05_peak_similarity.R`) | bmmc_d1 |
| s1d1 paired-multiome baseline matrix | `datasets/bmmc_d1/06_s1d1_paired/04_figures/plot_metrics_matrix.R` | `06_s1d1_paired/03_metrics/` outputs | bmmc_d1 |
| s1d1 paired 9-metric table | `datasets/bmmc_d1/06_s1d1_paired/04_figures/figS11_table_9metric.R` | `06_s1d1_paired/03_metrics/major/{sum_metrics,celltype_metrics,adj_atac_predaccu}.csv` + `peakdist.csv` | bmmc_d1 |

The s1d1 baseline is paired multiome, so it carries the full 9 metrics; the cross-site and
3-site experiments are independently sequenced and carry the 6 unpaired metrics.

### Fig S12 — Parse / 10x cross-platform matrix

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S12 (metrics matrix) | `datasets/pbmc_parse/04_figures/plot_metrics_matrix.R` | `<metrics>/{sum_metrics,celltype_metrics,adj_atac_predaccu,peakdist}.csv` from `03_metrics/00_build_label.py` → `01_prep_latents.py` → `02_compute_metrics.py` → `03_adjust_accuracy.py` → `04_peak_similarity.R` | pbmc_parse |
| per-method UMAPs | `datasets/pbmc_parse/04_figures/plot_umap.py` | per-method latent CSV (**shipped in `results/pbmc_parse/latents/`**) + `label.csv` + `knn_pred_label__<method>.csv` (**shipped in `results/pbmc_parse/knn_pred_label/`**) | pbmc_parse |

Preprocessing: `01_preprocess/00_data_prep.R` → `01_azimuth_annotate_parse.R` →
`02_compare_donor_composition.R` → `03_qc_confidence_donors.R` → `04_build_integration_dataset.R`.

### Fig S13 — scalability and ranking-weight sensitivity

| Panel | Script | Key inputs | Dataset |
|---|---|---|---|
| S13A (runtime and memory vs test-cell count) | `summary/04_figures/plot_figS13_scalability.R` | `results/summary/runtime_memory_sel_plus.csv` | all four |
| S13B (weighted rank over the 50–90% weight range) | `summary/04_figures/plot_figS13_weighting_rank.R` | `results/summary/weighting_sensitivity_ranks.csv` (shipped copy; regenerated by `summary/03_metrics/00_weighting_sensitivity.R`) | all four |

`00_weighting_sensitivity.R` reads the four datasets' `sum_metrics_clean.csv`,
`runtime_memory_sel_plus.csv` and `new_methods_grouped.csv`.

Superseded, not included in this repository: `plot_weighting_sensitivity_integer.R`. Its y-axis
reads "Final rank (1 = best)" with integer stepped lines; the published S13B reads "Weighted rank
(lower = better)" with continuous linear ranks, which is the shipped script's label verbatim.

---

## Table S4

| Item | Script | Key inputs | Dataset |
|---|---|---|---|
| Table S4 (paired Wilcoxon, batch vs no-batch, 11 metrics) | `datasets/bmmc_d1/04_figures/tableS4_batch_wilcoxon.R` | `<metrics>/{sum_metrics,celltype_metrics,adj_atac_predaccu,peakdist}.csv` covering the six `(batch)` variants and their no-batch counterparts | bmmc_d1 |

Writes `tableS4_batch_correction.csv`.

---

## Known gaps

The three panels that previously had no generating script here were resolved on 2026-08-31 by
recovering the originals from the pre-revision analysis directory:

| Figure | Resolution |
|---|---|
| Fig S7 | the three per-locus strip scripts now ship in `datasets/rms/04_figures/figS7_multiome_vs_annotation/`; the published figure is a manual assembly of their outputs |
| Fig 4C | ported to `datasets/brca/04_figures/plot_fig4c_confusion.py` from the original notebook, whose `Fig4C.pdf` output is md5-identical to the published asset |
| Fig 5C | same script and inputs as Figs S8/S9 (`figS7_S9/plot_pileup_grid_2page.R`) with the 7-row method selection documented under Fig 5 above |

Remaining caveats: the Fig S7 and Fig S8/S9 per-cell pileup inputs are not committed (size;
available from the authors on request), and the Fig S7 figure assembly itself (three strips into
one page) was manual, so no single script emits the composed figure.

---

## Composite score and metric definitions

The composite score is the mean of 6 metrics for independently sequenced (unpaired) data:
`ks_celltype_mean`, `asw`, `omics_asw`, `ks_inter_celltype_mean`, `average_accu`, `peakdist_adj`.
Paired data adds `ks.statistic`, `ari` and `ami`, for 9 metrics. ARI, AMI and the omics KS
statistic are same-cell metrics and are therefore undefined for independently sequenced datasets —
that is why the matrices differ between 9 and 6 columns.

KNN label transfer throughout: k = 10, cosine distance, distance-weighted.

Random adjustment is a normalisation, not a subtraction:
`(accu - random_accu) / (1 - random_accu)` for accuracy, and `(random - method) / random` for the
peak distance.

---

## How to regenerate a figure

1. **Set your paths.** From the repository root:

   ```bash
   cp config/config.local.sh.example config/config.local.sh
   # edit config/config.local.sh: DATA_ROOT, PROJECT_ROOT
   ```

   `config/config.sh` defines `DATA_ROOT`, `PROJECT_ROOT`, `BENCHMARK_FUN_DIR` (defaults to
   `common/`) and `PUBLISHED_REF` (legacy, defaults to `results/`), all with `/path/to/data`
   placeholders, and sources `config/config.local.sh` if present. `config/config.yaml` and
   `config/config.R` carry the same values for Python and R. R scripts read them through
   `Sys.getenv()`.

   `config/config.local.sh` is gitignored. Do not commit it.

2. **Make sure the inputs exist — usually they already do.** Every metrics-matrix figure script
   (pbmc3k Fig 2B, pbmc10k Fig S2C, pbmc_parse Fig S12, BRCA Fig S6, RMS Fig 5A, the three BMMC
   matrices, and the Fig 6A panel of `fig6_combined.py`) runs in **single-file mode by default**:
   it loads the shipped final merged table (`fig2b_matrix.csv` / `metrics_matrix.csv` /
   `fig6_scores.csv` in the dataset's `results/` folder) and renders the figure directly — no
   per-metric tables needed. Set `REBUILD_MATRIX=1` to instead rebuild the merged
   table from the per-metric CSVs plus the not-shipped reference tables (available on request). For other figures, if the
   inputs are not already on disk, run that dataset's `01_preprocess`, `02_methods` and
   `03_metrics` steps first; the numbered filenames give the order, and each dataset's
   `03_metrics/run_all.sh` chains them where one exists. Every figure input table ships under
   `results/`, except the two large RMS files noted under Fig S8/S9 (gitignored, available on request).

3. **Run the script named in the table.**

   ```bash
   . config/config.sh
   Rscript datasets/pbmc3k/04_figures/fig2b/plot_metrics_matrix.R <metrics_dir> <out_dir>
   python  datasets/bmmc_d1/04_figures/fig6_combined.py
   ```

   Scripts with a sibling `*_sub.sh` are normally submitted through it; that wrapper is an LSF
   `bsub` driver and sets the environment variables (for example `SPLICE_PUBLISHED`) that the
   script expects. Read the wrapper before running the script by hand.

**Status note.** Absolute cluster paths have been scrubbed from every script. Shell scripts source
`config/config.sh`; `.py` and `.R` scripts carry literal `/path/to/...` placeholders at the top of
the file rather than live config lookups, so edit that placeholder block (or set the matching
environment variables where the script reads them via `Sys.getenv()`) before running one by hand.
Scripts are syntax-checked but were not re-executed after scrubbing. Treat this file as the map of
what generates what, not as a claim that every script runs unmodified from a fresh checkout.
