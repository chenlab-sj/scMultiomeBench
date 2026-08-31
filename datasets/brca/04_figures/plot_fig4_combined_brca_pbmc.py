#!/usr/bin/env python3
"""
Combined Fig4 reproducibility across BOTH datasets (BRCA + pbmc3k), one row of 4 panels:
  A  BRCA   pairwise-NMI boxplot        C  BRCA   per-cell-type SD dot+bar
  B  pbmc3k pairwise-NMI boxplot        D  pbmc3k per-cell-type SD dot+bar
NMI pair (A,B) share a y-axis; SD pair (C,D) share a y-axis. Each SD panel keeps its own cell-type legend
(the two datasets have different cell types). Manuscript style (lightsteelblue boxes; gray bars + Dark2 dots).
    python plot_fig4_combined_brca_pbmc.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

FIG = os.path.dirname(os.path.abspath(__file__))
PBMC = os.path.normpath(os.path.join(FIG, "../../../pbmc/pbmc3k/benchmark/fig4/curated"))

# ---- inputs (Method/Pair/NMI  and  Method/CellType/Std) ----
brca_nmi = pd.read_csv(os.path.join(FIG, "nmi_df_7type.csv"))
brca_sd  = pd.read_csv(os.path.join(FIG, "sd_df_7type.csv"))
pbmc_nmi = pd.read_csv(os.path.join(PBMC, "nmi_df.csv"))
pbmc_sd  = pd.read_csv(os.path.join(PBMC, "sd_df.csv"))

BRCA_CT = ["Tumor", "Macrophages", "Fibroblasts", "Luminal progenitor", "T-cells", "Basal progenitor", "B-cells"]
PBMC_CT = ["CD14 Monocytes", "CD16 Monocytes", "B cells", "Memory T cells",
           "Naive CD4 T cells", "Naive CD8 T cells", "NK/Effector T cells"]

NMI_YLIM = (0.2, 1.08)          # shared across A,B (both fit; pbmc3k sits higher = more reproducible)
# SD panels use INDEPENDENT y (pbmc3k is much more reproducible -> a shared 0-0.2 axis flattens it to nothing)

# ---- font sizes (bumped up) ----
FS_TITLE, FS_YLAB, FS_XTICK, FS_YTICK, FS_LEG, FS_LEGT = 21, 19, 17, 16, 14, 15

def nmi_panel(ax, df, title, ylabel=None, order=None):
    if order is None:
        order = df.groupby("Method")["NMI"].mean().sort_values(ascending=False).index.tolist()
    sns.boxplot(data=df, x="Method", y="NMI", order=order, color="lightsteelblue",
                fliersize=2, linewidth=0.8, ax=ax)
    ax.set_ylim(*NMI_YLIM); ax.set_xlabel("")
    ax.set_ylabel(ylabel if ylabel else "", fontsize=FS_YLAB)
    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold")
    ax.set_xticklabels(order, rotation=90, ha="center", fontsize=FS_XTICK)
    ax.tick_params(axis="y", labelsize=FS_YTICK)

def sd_panel(ax, df, ct_order, title, ylabel=None, order=None, ylim=None):
    if order is None:
        order = df.groupby("Method")["Std"].mean().sort_values().index.tolist()   # most->least reproducible
    d = df.copy()
    d["CellType"] = pd.Categorical(d["CellType"], categories=ct_order, ordered=True)
    d["Method"] = pd.Categorical(d["Method"], categories=order, ordered=True)
    sns.barplot(data=d, x="Method", y="Std", color="lightgray", ci=None, ax=ax)   # seaborn 0.11.x
    sns.stripplot(data=d, x="Method", y="Std", hue="CellType", hue_order=ct_order,
                  order=order, dodge=True, size=4, alpha=0.85, palette="Dark2", ax=ax)
    ax.set_ylim(*(ylim if ylim else (0, df["Std"].max() * 1.18))); ax.set_xlabel("")
    ax.set_ylabel(ylabel if ylabel else "", fontsize=FS_YLAB)
    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold")
    ax.set_xticklabels(order, rotation=90, ha="center", fontsize=FS_XTICK)
    ax.tick_params(axis="y", labelsize=FS_YTICK)
    ax.legend(loc="upper left", frameon=True, fontsize=FS_LEG, title="Cell type", title_fontsize=FS_LEGT)

# pbmc3k first, then BRCA; each panel ordered by its OWN value. NMI pair shares a y-axis and the SD pair
# shares a y-axis (facet style: y-ticks only on the left panel of each pair).
SD_YLABEL = "Std. dev. of ATAC cell type\nprediction across triplicates"
sd_shared = (0, max(brca_sd["Std"].max(), pbmc_sd["Std"].max()) * 1.10)

fig = plt.figure(figsize=(20, 5.8))
# 5 columns: [NMI-pbmc | NMI-BRCA | spacer | SD-pbmc | SD-BRCA]. Small wspace -> each shared-y pair sits
# close (facet look); the empty spacer column keeps the two metric groups apart.
gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.28, 1, 1], wspace=0.06)
axes = [fig.add_subplot(gs[0, c]) for c in (0, 1, 3, 4)]
nmi_panel(axes[0], pbmc_nmi, "pbmc3k", ylabel="Pairwise NMI among triplicates")
nmi_panel(axes[1], brca_nmi, "BRCA")
sd_panel(axes[2], pbmc_sd, PBMC_CT, "pbmc3k", ylabel=SD_YLABEL, ylim=sd_shared)
sd_panel(axes[3], brca_sd, BRCA_CT, "BRCA", ylim=sd_shared)

for ax in (axes[1], axes[3]):        # facet style: hide inner y tick labels (shared scale within each pair)
    ax.tick_params(axis="y", labelleft=False)
plt.savefig(os.path.join(FIG, "fig4_combined.pdf"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(FIG, "fig4_combined.png"), dpi=140, bbox_inches="tight")
print("wrote fig4_combined.pdf/png")
print("NMI order BRCA :", brca_nmi.groupby('Method')['NMI'].mean().sort_values(ascending=False).round(3).to_dict())
