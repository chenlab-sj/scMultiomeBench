#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
Fig4D macrophage-subset: compute KNN ATAC-label accuracy for the 3 NEW methods
(MaxFuse, MIDAS, scButterfly) across the 5 subsampling replicate series (macrosub,
macrosub2..5) x subset fractions (sub0=full ... sub5=1/32).

The 11 OLD methods already have their per-(rep,sub) accuracy CSVs
(knn_k10sub{0..5}_pred_accu.csv) and are NOT recomputed -- the plotting notebook reads
those directly. This script only fills the gap for the 3 new methods, which have subset
latents under BRCA/HT243-S1H4_subsample/HT243_S1H4_macrosub{rep}/{method}/sub{s}/latent.csv
(+ full-data latent under BRCA/HT243B1-S1H4/{method}/latent.csv for the sub0 point).

Metric = exactly the old macophage_sub_fig4.ipynb knn_ataclabel:
  KNeighborsClassifier(n_neighbors=10, metric="cosine", weights="distance"),
  fit on the test RNA cells present in the latent, predict the test ATAC cells,
  class-wise accuracy for a cell type = confusion-matrix diagonal / row-sum
  = recall among that type's true ATAC cells. Subsampling is encoded purely in
  which barcodes are present in each latent (global label file, unchanged).

Runs anywhere with pandas + scikit-learn (no scanpy needed). Writes:
  new_methods_accu_full.csv  -- overall + every cell type, long format (archival)
  new_methods_macro_accu.csv -- Macrophages only, columns matched to the old macro_accu
"""
import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, confusion_matrix

G = "/path/to/multiomeBench"
LABEL      = f"{G}/BRCA/benchmark/fig4/label.csv"          # == notebook's HT243B1-S1H4_adj/label.csv
SUBSAMPLE  = f"{G}/BRCA/HT243-S1H4_subsample"
FULLDIR    = f"{G}/BRCA/HT243B1-S1H4"                       # full-data latents (sub0)
OUT        = f"{G}/BRCA/benchmark/fig4d_macrophage_subset"
K = 10

# on-disk dir name -> display name (matches 02_reproduce_metrics.py naming)
NEW = {"maxfuse": "MaxFuse", "midas": "MIDAS", "scbutterfly": "scButterfly"}

# ---- global label (test rows are <bc>-1_rna / <bc>-1_atac; train rows are plain <bc>-1) ----
labels_annot = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test = labels_annot.index.values[labels_annot["modality"] != "train multiomics"]
bc_test_rna  = [b for b in bc_test if str(b).endswith("_rna")]
bc_test_atac = [b for b in bc_test if str(b).endswith("_atac")]
ALL_TYPES = sorted(labels_annot["cell_type"].unique())     # explicit label order -> robust per-type recall


def rep_dir(rep):
    return f"{SUBSAMPLE}/HT243_S1H4_macrosub" if rep == 1 else f"{SUBSAMPLE}/HT243_S1H4_macrosub{rep}"


def latent_path(method_dir, rep, sub):
    """sub0 = full-data latent (rep1 only); sub1..5 = the subsampled latent for that replicate."""
    if sub == 0:
        return f"{FULLDIR}/{method_dir}/latent.csv"
    return f"{rep_dir(rep)}/{method_dir}/sub{sub}/latent.csv"


def knn_accu(lat_df):
    """Return (overall_acc, {cell_type: recall}) for one latent, matching knn_ataclabel."""
    rna  = [b for b in bc_test_rna  if b in lat_df.index]
    atac = [b for b in bc_test_atac if b in lat_df.index]
    clf = KNeighborsClassifier(n_neighbors=K, metric="cosine", weights="distance")
    clf.fit(lat_df.loc[rna], labels_annot.loc[rna, "cell_type"])
    pred = clf.predict(lat_df.loc[atac])
    true = labels_annot.loc[atac, "cell_type"].values
    overall = accuracy_score(true, pred)
    cm = confusion_matrix(true, pred, labels=ALL_TYPES).astype(float)
    rowsum = cm.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        recall = np.divide(np.diagonal(cm), rowsum, out=np.full(len(ALL_TYPES), np.nan), where=rowsum > 0)
    return overall, dict(zip(ALL_TYPES, recall))


full_rows, macro_rows = [], []
for method_dir, disp in NEW.items():
    for rep in range(1, 6):
        subs = range(0, 6) if rep == 1 else range(1, 6)   # sub0 (full) only contributes once, like the old rep1
        for sub in subs:
            path = latent_path(method_dir, rep, sub)
            if not os.path.exists(path):
                print(f"  MISSING {disp} rep{rep} sub{sub}: {path}")
                continue
            lat = pd.read_csv(path, index_col=0)
            overall, recall = knn_accu(lat)
            n_atac = sum(b in lat.index for b in bc_test_atac)
            n_macro = int(((labels_annot.loc[[b for b in bc_test_atac if b in lat.index], "cell_type"]
                            == "Macrophages").sum()))
            pct = 1 / 2 ** sub
            full_rows.append([disp, overall, "overall", sub, rep, pct])
            for ct, val in recall.items():
                full_rows.append([disp, val, ct, sub, rep, pct])
            macro_rows.append([disp, recall["Macrophages"], "Macrophages", sub, rep, pct])
            print(f"  {disp:11s} rep{rep} sub{sub}: overall={overall:.3f} "
                  f"macro_recall={recall['Macrophages']:.3f} (n_atac={n_atac}, n_macro_atac={n_macro})")

cols = ["pipeline", "accuracy", "type", "sub", "rep", "percentage"]
pd.DataFrame(full_rows, columns=cols).to_csv(f"{OUT}/new_methods_accu_full.csv", index=False)
pd.DataFrame(macro_rows, columns=cols).to_csv(f"{OUT}/new_methods_macro_accu.csv", index=False)
print(f"\nDONE: wrote new_methods_accu_full.csv ({len(full_rows)} rows) and "
      f"new_methods_macro_accu.csv ({len(macro_rows)} rows) to {OUT}")
