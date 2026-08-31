#!/usr/bin/env python3
"""BMMC Fig2b step 3a: KNN ATAC label transfer (k=10, cosine, distance-weighted) per method, over the
same staged latents + label as 02_run_bmmc_metrics.py (proven 36198/36198 overlap). Writes one
knn_pred_label__<method>.csv per method (index = ATAC barcode, single column = method -> predicted
cell type) for 01_adjust_accuracy.py (average_accu) AND the peak step, plus knn_pred_accu.csv
(method, accuracy, type) for reference.  python 03_knn_bmmc.py"""
import os
import re
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, confusion_matrix

HERE  = os.path.dirname(os.path.abspath(__file__))
STAGE = os.path.join(HERE, "staged")
LABEL = os.environ.get("LABEL", os.path.join(HERE, "old", "BMMC_d1", "label.csv"))
OUT   = os.environ.get("KNN_OUT", HERE)
K     = 10

PIPELINES = ["BindSC", "Seurat(CCA)", "MaxFuse", "MIDAS", "scBridge", "scVI", "scButterfly", "scDART",
             "scglue", "scglue(multiome)", "scJoint", "simba", "Portal",
             "BindSC(batch)", "scVI(batch)", "scglue(batch)", "scglue(multiome,batch)",
             "scJoint(batch)", "MIDAS(batch)"]

labels = pd.read_csv(LABEL, index_col=0)
bc_test = labels[labels["set"] != "train multiomics"].index.values
bc_rna  = [b for b in bc_test if re.search(r"_rna.*$",  str(b))]
bc_atac = [b for b in bc_test if re.search(r"_atac.*$", str(b))]
print(f"label: {len(labels)} rows; test rna={len(bc_rna)} atac={len(bc_atac)}", flush=True)

accu_rows = []
for m in PIPELINES:
    f = os.path.join(STAGE, m, "latent.csv")
    if not os.path.exists(f):
        print("SKIP (missing staged):", m); continue
    lat = pd.read_csv(f, index_col=0)
    rna  = [b for b in bc_rna  if b in lat.index]
    atac = [b for b in bc_atac if b in lat.index]
    clf = KNeighborsClassifier(n_neighbors=K, metric="cosine", weights="distance")
    clf.fit(lat.loc[rna], labels.loc[rna, "cell_type"])
    pred = clf.predict(lat.loc[atac])
    pd.DataFrame({m: pred}, index=atac).to_csv(os.path.join(OUT, f"knn_pred_label__{m}.csv"))

    true = labels.loc[atac, "cell_type"].values
    accu_rows.append([m, float(accuracy_score(true, pred)), "overall"])
    cats = sorted(pd.unique(np.concatenate([true, pred])))
    cm = confusion_matrix(true, pred, labels=cats).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        recall = np.diagonal(cm) / cm.sum(axis=1)
    present = set(np.unique(true))
    for ct, r in zip(cats, recall):
        if ct in present:
            accu_rows.append([m, float(r), ct])
    print(f"  {m}: acc={accu_rows[-len(present)-1][1]:.3f}  (rna {len(rna)}, atac {len(atac)})", flush=True)

pd.DataFrame(accu_rows, columns=["pipeline", "accuracy", "type"]).to_csv(
    os.path.join(OUT, "knn_pred_accu.csv"), index=False)
print(f"DONE: {sum(os.path.exists(os.path.join(OUT, f'knn_pred_label__{m}.csv')) for m in PIPELINES)} "
      f"knn_pred_label__*.csv + knn_pred_accu.csv in {OUT}", flush=True)
