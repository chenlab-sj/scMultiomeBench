#!/usr/bin/env python
"""
umap.py -- per-method UMAP panels for the pbmc3k benchmark (standalone, one method/call).

Faithful to the umap cell of benchmarkmatrix_sel.ipynb, but a single 1x5 row per method:
panels = [multiomics, scRNA, scATAC, Cell type, Predicted ATAC cell type].
The "Predicted ATAC cell type" panel uses this method's KNN predicted ATAC labels
(knn_pred_label__<method>.csv from compute_metrics.py).

    python umap.py --latent scripts/MIRA/latent.csv --method MIRA \
                   --label metrics/label.csv --knn-pred metrics/knn_pred_label__MIRA.csv --out metrics
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

OMICS = ["multiomics", "scRNA", "scATAC"]
PANEL = ["multiomics", "scRNA", "scATAC", "Cell type", "Predicted ATAC cell type"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--knn-pred", required=True,
                    help="this method's predicted-ATAC-label CSV (index=atac bc, 1 col)")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    labels = pd.read_csv(args.label, index_col=0).dropna(subset=["cell_type"])
    bc_sel = labels.index.values
    labels_test = labels[labels["modality"] != "train multiomics"]
    celltype = labels_test["cell_type"].unique().tolist()

    lat_df = pd.read_csv(args.latent, index_col=0)
    lat_df = lat_df.loc[lat_df.index.isin(bc_sel)]

    # this method's predicted ATAC labels -> a Series named like the panel
    knn = pd.read_csv(args.knn_pred, index_col=0)
    knn_atac = knn.iloc[:, 0]
    knn_atac.name = "Predicted ATAC cell type"

    adata = sc.AnnData(lat_df)
    sc.pp.neighbors(adata)
    adata.obs = pd.merge(adata.obs, knn_atac, how="left", left_index=True, right_index=True)
    sc.tl.umap(adata)

    ann = labels.reindex(adata.obs.index)
    adata.obs["Cell type"] = pd.Categorical(ann["cell_type"], categories=celltype, ordered=True)
    adata.obs["Predicted ATAC cell type"] = pd.Categorical(
        adata.obs["Predicted ATAC cell type"], categories=celltype, ordered=True)
    mod = ann["modality"]
    adata.obs["multiomics"] = pd.Categorical(
        mod.apply(lambda x: "multiomics" if x == "train multiomics" else ""), categories=OMICS, ordered=True)
    adata.obs["scRNA"] = pd.Categorical(
        mod.apply(lambda x: "scRNA" if x == "test scRNA" else ""), categories=OMICS, ordered=True)
    adata.obs["scATAC"] = pd.Categorical(
        mod.apply(lambda x: "scATAC" if x == "test scATAC" else ""), categories=OMICS, ordered=True)

    fig, axs = plt.subplots(1, len(PANEL), figsize=(16, 3.2))
    plt.rc("font", size=14)
    handles1 = handles2 = labels1 = labels2 = None
    for j, name in enumerate(PANEL):
        ax = axs[j]
        if j <= 2:
            sc.pl.umap(adata, color=name, palette=sns.color_palette("Set2"),
                       show=False, ax=ax, title=None, frameon=False)
        else:
            sc.pl.umap(adata, color=name, show=False, ax=ax, title=None, frameon=False)
        if j == 0:
            ax.axis("on")
            ax.set_ylabel("\n" + args.method, rotation=90, fontsize=16)
            ax.set_xlabel("")
            ax.set(frame_on=False)
            handles1, labels1 = ax.get_legend_handles_labels()
        if j == 3:
            handles2, labels2 = ax.get_legend_handles_labels()
        ax.set_title(name, fontsize=14)
        if ax.legend_ is not None:
            ax.legend_.remove()

    fig.text(0.5, 0.01, "UMAP1", ha="center", fontsize=14)
    fig.text(0.01, 0.5, "UMAP2", va="center", rotation="vertical", fontsize=14)
    if handles1:
        fig.legend(handles1, labels1, loc="upper left", bbox_to_anchor=(0.02, 0.99),
                   ncol=2, fontsize=12, title=PANEL[0])
    if handles2:
        fig.legend(handles2, labels2, loc="upper left", bbox_to_anchor=(0.40, 0.99),
                   ncol=3, fontsize=12, title=PANEL[3])
    fig.tight_layout(rect=(0.02, 0.05, 0.98, 0.80))
    out = os.path.join(args.out, f"{args.method}_umap.png")
    plt.savefig(out, dpi=300)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
