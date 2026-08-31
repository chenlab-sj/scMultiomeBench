#!/usr/bin/env python3
"""
pbmc3k k-sensitivity new-method k-test: per-cell-type KNN ATAC-label-prediction accuracy at k=5/10/20/40/80
for the 5 new methods (MaxFuse/MIDAS/scButterfly/MIRA/Multigrate), EXACTLY as 00_compute_ktest_legacy18.py::knn_ataclabel
(fit KNN on test RNA cells, predict test ATAC cells; class-wise accuracy = confusion-matrix diagonal /
row-sums; KNeighborsClassifier(n_neighbors=k, metric='cosine', weights='distance')). Appends to
new_methods_ktest_long.csv so plot_ksens.py splices them onto the published 18-method kNN_celltype_accu_sum.csv.

Reads each method's latent from pbmc/pbmc3k/scripts/<dir>/latent.csv (skips ones not present). Needs sklearn/
pandas (benchmark_env on cluster, or Mac).  python 01_compute_ktest_new_methods.py
"""
import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, confusion_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
P3K = os.path.abspath(os.path.join(HERE, "..", ".."))          # pbmc/pbmc3k (portable Mac/cluster)
LABEL = os.path.join(HERE, "label.csv")                         # Barcode,cell_type,modality,... (pbmc3k)
KS = [5, 10, 20, 40, 80]
METHODS = {"MaxFuse": "maxfuse", "MIDAS": "midas", "scButterfly": "scbutterfly",
           "MIRA": "MIRA", "Multigrate": "multigrate"}
OUT = os.path.join(HERE, "new_methods_ktest_long.csv")

labels_annot = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test = labels_annot[labels_annot["modality"] != "train multiomics"].index.values
bc_test_rna = [b for b in bc_test if b.endswith("_rna")]
bc_test_atac = [b for b in bc_test if b.endswith("_atac")]

def knn_ataclabel(k, lat_df, pipeline):
    """VERBATIM 00_compute_ktest_legacy18.py recipe so new methods match the published 18."""
    rna_valid = [i for i in bc_test_rna if i in lat_df.index]
    atac_valid = [i for i in bc_test_atac if i in lat_df.index]
    lat_rna, lat_atac = lat_df.loc[rna_valid], lat_df.loc[atac_valid]
    lab_rna = labels_annot.loc[labels_annot.index.isin(lat_rna.index)][["cell_type"]]
    lab_atac = labels_annot.loc[labels_annot.index.isin(lat_atac.index)][["cell_type"]]
    clf = KNeighborsClassifier(n_neighbors=k, metric="cosine", weights="distance")
    clf.fit(lat_rna, lab_rna.values.ravel())
    pred = clf.predict(lat_atac)
    rows = [[pipeline, accuracy_score(lab_atac, pred), "overall"]]
    conf = confusion_matrix(lab_atac, pred)
    cwa = np.diagonal(conf) / np.sum(conf, axis=1)
    types = np.unique(pred)
    for j in range(len(types)):
        rows.append([pipeline, cwa[j], types[j]])
    return rows

records = []
for method, sub in METHODS.items():
    lat_path = os.path.join(P3K, "scripts", sub, "latent.csv")
    if not os.path.exists(lat_path):
        print(f"SKIP {method}: no latent ({lat_path})"); continue
    lat = pd.read_csv(lat_path, index_col=0)
    for k in KS:
        for pl, acc, typ in knn_ataclabel(k, lat, method):
            records.append({"pipeline": pl, "accuracy": acc, "type": typ, "k value": k})
    print(f"OK   {method}: {lat.shape[0]} cells")

if records:
    out = pd.DataFrame(records)[["pipeline", "accuracy", "type", "k value"]]
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}  ({out.pipeline.nunique()} methods x {len(KS)} k)")
else:
    print("no latents found yet.")
