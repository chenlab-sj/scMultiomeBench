#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
Fig4C = BRCA (HT243) averaged triplicate confusion-matrix heatmaps for Portal / Seurat(CCA) / scVI.
For each method, the row-normalised ATAC-label confusion matrix (actual vs KNN-predicted, k=10)
is averaged over its 3 replicate runs (coolwarm, vmin=-1, vmax=1, auto-contrast annotations),
plus a right-side bar of true per-cell-type test-ATAC cell counts.

This is the Fig 4C generator, ported from the original analysis notebook (macro_sub-Fig4C1.ipynb,
Jan 2026); the published Fig4C.pdf was verified md5-identical to that notebook's output.
The functions are kept verbatim -- only the two input paths are repointed to the repo copies:
  - predictions: datasets/brca/03_metrics/reproducibility/rep_knn_k10_pred_label_7type.csv
    (written by 03_fig4_reproduce_7type.ipynb; a column SUPERSET of the original
    rep_knn_k10_pred_label.csv -- same '<Method>-1/-2/-3' naming, same 7-type labels; the extra
    methods' columns are ignored because plot_avg_triplicate_heatmaps only plots the requested
    `methods`)
  - labels: results/brca/label.csv (the 7-cell-type truth)

Run:  python plot_fig4c_confusion.py     -> Fig4C.pdf in the current directory
"""
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix


def _avg_confusion_percent(df, actual_col, pred_cols, labels):
    """
    Compute confusion matrix (row-normalized %) for each pred_col, then average across cols.
    Returns: avg_percent [C x C], and per-replicate matrices list.
    """
    actual = df[actual_col].astype(str).values

    mats = []
    for c in pred_cols:
        pred = df[c].astype(str).values
        cm = confusion_matrix(actual, pred, labels=labels).astype(float)

        # row-normalize to percentages; protect against empty rows
        row_sums = cm.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        cm_pct = cm / row_sums
        mats.append(cm_pct)

    avg = np.mean(mats, axis=0)
    return avg, mats

def plot_avg_triplicate_heatmaps(
    KNN_atac,
    methods=("Portal", "Seurat(CCA)", "scVI"),
    actual_col="cell_type",
    rep_pattern=r"^(?P<method>.+)-(?P<rep>[123])$",
    label_map=None,             # dict: raw_label -> pretty_label (optional)
    figsize=(9, 4),
    dpi=300,
    save=None
):
    """
    Plot average confusion heatmap per method (averaged across -1/-2/-3 columns),
    plus a right-side bar plot of true-class counts.
    """
    df = KNN_atac.copy()

    # ---- define labels (sorted by frequency in actual) ----
    actual = df[actual_col].astype(str).values
    categories, counts = np.unique(actual, return_counts=True)
    sorted_idx = np.argsort(counts)[::-1]
    labels = categories[sorted_idx]
    label_counts = counts[sorted_idx]

    # pretty labels
    if label_map is None:
        pretty = {x: x for x in labels}
    else:
        pretty = {x: label_map.get(x, x) for x in labels}
    pretty_labels = [pretty[x] for x in labels]

    # ---- find replicate columns for each method ----
    # build mapping: method -> [col1,col2,col3]
    col_method = {}
    rx = re.compile(rep_pattern)

    for c in df.columns:
        m = rx.match(c)
        if m:
            meth = m.group("method")
            col_method.setdefault(meth, []).append(c)

    # sanity: only keep requested methods
    method_cols = {}
    for meth in methods:
        cols = sorted(col_method.get(meth, []))
        # keep only -1/-2/-3 if present
        cols = [c for c in cols if rx.match(c)]
        if len(cols) == 0:
            raise ValueError(f"No replicate columns found for method '{meth}'. Expected like '{meth}-1'.")
        method_cols[meth] = cols

    # ---- make figure (heatmaps + count bar) ----
    n_methods = len(methods)
    width_ratios = [2.7] * n_methods + [1.2]
    fig, axes = plt.subplots(
        1, n_methods + 1,
        figsize=figsize,
        gridspec_kw={"width_ratios": width_ratios},
        dpi=dpi
    )
    *heat_axes, ax_bar = axes

    # ---- plot each averaged heatmap ----
    for ax, meth in zip(heat_axes, methods):
        avg_pct, _ = _avg_confusion_percent(df, actual_col, method_cols[meth], labels)

        im = ax.imshow(avg_pct, aspect="auto", cmap="coolwarm", vmin = -1, vmax = 1)  # 0..1 proportions

        #ax.set_title(f"{meth}\n(mean of {len(method_cols[meth])} runs)", fontsize=12)
        ax.set_title(meth, fontsize=12)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))

        ax.set_xticklabels(pretty_labels, rotation=90, ha="right", fontsize=10)
        ax.set_yticklabels(pretty_labels, fontsize=10)

        # annotate values
 # annotate values with auto-contrast text color
        norm = im.norm
        cmap = im.get_cmap()

        for i in range(avg_pct.shape[0]):
            for j in range(avg_pct.shape[1]):
                val = avg_pct[i, j]

                # background RGBA from colormap
                r, g, b, _ = cmap(norm(val))

                # perceived luminance (0=dark, 1=bright)
                luminance = 0.299 * r + 0.587 * g + 0.114 * b

                txt_color = "white" if luminance < 0.5 else "black"

                ax.text(j, i, f"{val:.2f}",
                        ha="center", va="center",
                        fontsize=8, color=txt_color)


        # only left-most keeps y tick labels
        if ax is not heat_axes[0]:
            ax.set_yticklabels([])
            ax.set_yticks(np.arange(len(labels)))

    heat_axes[0].set_ylabel("Actual ATAC cell labels", fontsize=12)
    heat_axes[-2].set_xlabel("Predicted ATAC cell labels", fontsize=12)

    # colorbar for all heatmaps
    #cbar = fig.colorbar(im, ax=heat_axes, fraction=0.025, pad=0.02)
    #cbar.set_label("Row-normalized fraction", fontsize=11)

    # ---- right-side counts bar ----
    ax_bar.barh(pretty_labels[::-1], label_counts[::-1],color = "tab:red")
    ax_bar.set_xlabel("Cell Counts", fontsize=12)
    ax_bar.set_yticks([])  # hide y labels (already on heatmaps)
    ax_bar.set_xlim(0, label_counts.max() * 1.1)

    plt.tight_layout()

    if save:
        plt.savefig(save, bbox_inches="tight")
    return fig


if __name__ == "__main__":
    # ---- inputs (repo copies of the original cluster files) ----
    # triplicate KNN (k=10) ATAC-label predictions, columns '<Method>-1/-2/-3'
    KNN_annot = '/path/to/multiomeBench/datasets/brca/03_metrics/reproducibility/rep_knn_k10_pred_label_7type.csv'
    # 7-cell-type truth labels (cell_type + modality per barcode)
    cell_annot = '/path/to/multiomeBench/results/brca/label.csv'

    KNN_atac_label = pd.read_csv(KNN_annot, index_col = 0)

    # test ATAC barcodes = not 'train multiomics', suffix '_atac'; merge truth with predictions
    labels_annot= pd.read_csv(cell_annot,index_col=0)
    bc_test=labels_annot.index.values[~labels_annot['modality'].isin(['train multiomics'])]
    bc_test_atac =[bc for bc in bc_test if bc.endswith("_atac")]
    KNN_atac=pd.merge(labels_annot.loc[bc_test_atac], KNN_atac_label, left_index=True, right_index=True)
    print('KNN_atac:', KNN_atac.shape)

    fig = plot_avg_triplicate_heatmaps(
        KNN_atac,
        methods=("Portal", "Seurat(CCA)", "scVI"),
        actual_col="cell_type",
        save="Fig4C.pdf"
    )
    print('wrote Fig4C.pdf')
