#!/usr/bin/env python3
"""
FigS2 + FigS3A = BRCA (HT243) triplicate reproducibility confusion-matrix grid. For each method, the
row-normalised ATAC-label confusion matrix (actual vs KNN-predicted, k=10) for its 3 replicate runs
+ a (std) panel = per-cell SD across the 3 reps. Same recipe/style as the published figure
(coolwarm, vmin=-1, vmax=1, center=0; std in Purples). Adds the 3 new methods (MaxFuse/MIDAS/scButterfly);
their reps are already in rep_knn_k10_pred_label_7type.csv (the Fig4-selected triplets).

Data: rep_knn_k10_pred_label_7type.csv (cols <method>-1/-2/-3, rows = common ATAC cells) + label.csv (truth).
Split into two panels via FIGS2 / FIGS3A below (FigS3B distance heatmaps are a separate/later plot).
Run:  python plot_figS4_S5_confusion.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
PRED = pd.read_csv(os.path.join(HERE, "rep_knn_k10_pred_label_7type.csv"), index_col=0)
LABEL = os.path.join(HERE, "..", "old", "HT243B1-S1H4", "label.csv")
if not os.path.exists(LABEL):
    LABEL = os.path.join(HERE, "label.csv")
lab = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])

Actual = lab.loc[PRED.index, "cell_type"].values
CATS = pd.Series(Actual).value_counts().index.tolist()          # 7 types, by frequency (row/col order)
LBL = list(CATS)                                                # single-line labels (font 10 fits: y horizontal, x vertical)

# method -> its 3 rep columns
methods = sorted({c.rsplit("-", 1)[0] for c in PRED.columns})
reps_of = {m: [f"{m}-{r}" for r in (1, 2, 3) if f"{m}-{r}" in PRED.columns] for m in methods}

# method order = the BRCA composite-metrics table (figS4a) integration-performance ranking, best->worst,
# so the reproducibility grid reads in the same order as the metrics figure. Top 7 -> FigS2, next 7 -> FigS3A.
# (FigS3B distance heatmaps moved to a later figure.)
METRICS_ORDER = ["scJoint", "scglue", "scglue(multiome)", "scBridge", "Portal", "MaxFuse", "simba",
                 "Seurat(CCA)", "BindSC", "scVI", "MIDAS", "scButterfly", "Cobolt", "scDART"]
FIGS2  = METRICS_ORDER[:7]
FIGS3A = METRICS_ORDER[7:]

def cms_and_std(method):
    cms = []
    for rep in reps_of[method]:
        cm = pd.crosstab(pd.Series(Actual, name="a"), PRED[rep].values).reindex(index=CATS, columns=CATS).fillna(0).values.astype(float)
        rs = cm.sum(axis=1, keepdims=True); rs[rs == 0] = 1
        cms.append(cm / rs)
    return cms, np.std(np.stack(cms, 0), axis=0)

PAGE = (8.5, 11)   # US Letter, portrait -- each panel is one letter page

def draw_panel(method_list, out_stub, title):
    n = len(method_list)
    fig, axs = plt.subplots(n, 4, figsize=PAGE)
    if n == 1:
        axs = axs[None, :]
    for i, m in enumerate(method_list):
        cms, cm_std = cms_and_std(m)
        last = (i == n - 1)
        for j in range(3):
            ax = axs[i, j]
            data = cms[j] if j < len(cms) else np.zeros((len(CATS), len(CATS)))
            sns.heatmap(data, annot=True, fmt=".2f", cbar=False, cmap="coolwarm",
                        vmin=0, vmax=1, linewidths=.3, annot_kws={"size": 6},
                        xticklabels=(LBL if last else False),
                        yticklabels=(LBL if j == 0 else False), ax=ax)
            ax.set_title(f"{m}-{j+1}", fontsize=10)
            if last:
                ax.set_xticklabels(LBL, rotation=90, ha="center", fontsize=10)
            if j == 0:
                ax.set_yticklabels(LBL, rotation=0, fontsize=10)
        ax = axs[i, 3]
        sns.heatmap(cm_std, annot=True, fmt=".2f", cbar=False, cmap="Purples",
                    vmin=0, vmax=0.5, linewidths=.3, annot_kws={"size": 6},
                    xticklabels=(LBL if last else False), yticklabels=False, ax=ax)
        ax.set_title(f"{m} (std)", fontsize=10)
        if last:
            ax.set_xticklabels(LBL, rotation=90, ha="center", fontsize=10)
    fig.supylabel("Actual ATAC cell type label", fontsize=10)
    fig.supxlabel("Predicted ATAC cell type label", fontsize=10)
    fig.tight_layout(rect=(0.02, 0.01, 1, 1), h_pad=0.35, w_pad=0.2)
    fig.savefig(f"{HERE}/{out_stub}.pdf")                       # exact 8.5x11 letter page
    fig.savefig(f"{HERE}/{out_stub}.png", dpi=200)
    plt.close(fig)
    print(f"wrote {out_stub}.pdf/png  ({title}: {', '.join(method_list)})")

draw_panel(FIGS2,  "figS2",  "FigS2")
draw_panel(FIGS3A, "figS3a", "FigS3A")
