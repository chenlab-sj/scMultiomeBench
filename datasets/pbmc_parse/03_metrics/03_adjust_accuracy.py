#!/usr/bin/env python
"""
01_adjust_accuracy.py -- random-adjusted ATAC cell-type prediction accuracy (average_accu).

Faithful port of old/pbmc3k/benchmarkmatrix.ipynb (cells 42-47), parameterized.
Combines the per-method KNN predicted-label files written by 00_compute_metrics.py
(knn_pred_label__<method>.csv, index = ATAC barcodes) with label.csv's ground-truth
'cell_type' and the 'random_atac' random-baseline column, then per cell type:
    accu(method, ct) = accuracy_score(true==ct subset:  cell_type vs predicted)   [per-class recall]
    adjusted        = (accu - random_accu) / (1 - random_accu)     # random_atac as baseline
    average_accu    = mean over cell types of the adjusted values
Output: adj_atac_predaccu.csv  (index method; per-celltype adjusted cols + average_accu).

    python 01_adjust_accuracy.py --label label.csv --out .
"""
import argparse
import glob
import os
from functools import reduce

import pandas as pd
from sklearn.metrics import accuracy_score


def cls_accu(knn_atac, method, celltypes):
    row = [method]
    # score each method only on cells it actually predicted -> drop NaN (the outer-join fills
    # NaN for a method, e.g. Conos, on barcodes it lacks; mixing NaN(float) + str labels otherwise
    # crashes accuracy_score). nan accuracy where the method predicted no cell of a true type.
    sub = knn_atac[["cell_type", method]].dropna(subset=[method])
    for cls in celltypes:
        mask = sub["cell_type"] == cls
        row.append(accuracy_score(sub["cell_type"][mask], sub[method][mask])
                   if mask.any() else float("nan"))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--pred-glob", default="knn_pred_label__*.csv",
                    help="glob for per-method predicted-label CSVs from 00_compute_metrics.py")
    ap.add_argument("--out", default=".", help="dir to write adj_atac_predaccu.csv")
    args = ap.parse_args()

    lab = pd.read_csv(args.label, index_col=0)
    # test-ATAC ground truth + the random baseline column (matches the old notebook)
    atac = lab[lab["modality"] == "test scATAC"][["cell_type", "random_atac"]]

    pred_files = sorted(glob.glob(os.path.join(args.out, args.pred_glob)))
    if not pred_files:
        raise SystemExit(f"no predicted-label files matching {args.pred_glob} in {args.out}")
    preds = [pd.read_csv(f, index_col=0) for f in pred_files]
    pred_wide = reduce(lambda l, r: l.join(r, how="outer"), preds)   # one column per method

    knn_atac = atac.join(pred_wide, how="inner")
    methods = ["random_atac"] + list(pred_wide.columns)
    celltypes = knn_atac["cell_type"].unique().tolist()

    accu = pd.DataFrame([cls_accu(knn_atac, m, celltypes) for m in methods],
                        columns=["method"] + celltypes).set_index("method")
    random_accu = accu.iloc[0].values                       # the random_atac row
    adj = accu.apply(lambda r: (r - random_accu) / (1 - random_accu), axis=1)
    adj["average_accu"] = adj[celltypes].mean(axis=1)       # mean over celltype cols only
    out = os.path.join(args.out, "adj_atac_predaccu.csv")
    adj.to_csv(out)
    print(f"wrote {out}  ({len(methods)} methods incl. random_atac baseline)")


if __name__ == "__main__":
    main()
