#!/usr/bin/env python3
"""
Combined RMS Fig4 (single dataset): pairwise-NMI boxplot (left) + per-cell-type SD dot+bar (right).
Same manuscript style as the BRCA+pbmc3k fig4_combined -- large fonts, vertical (90deg) method labels,
lightsteelblue boxes, gray bars + Dark2 per-cell-type dots. Reads the rebuilt nmi_df.csv / sd_df.csv.
    python plot_fig4_combined_brca_pbmc.py
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

FIG = os.path.dirname(os.path.abspath(__file__))
nmi_df = pd.read_csv(os.path.join(FIG, "nmi_df.csv"))
sd_df  = pd.read_csv(os.path.join(FIG, "sd_df.csv"))
RMS_CT = ["Mesoderm", "Myoblast", "Myocyte"]                       # RMS major cell types (Dark2 order)

# ---- font sizes (match the BRCA/pbmc3k fig4_combined) ----
FS_YLAB, FS_XTICK, FS_YTICK, FS_LEG, FS_LEGT = 19, 17, 16, 14, 15

def nmi_panel(ax, df):
    order = df.groupby("Method")["NMI"].mean().sort_values(ascending=False).index.tolist()
    sns.boxplot(data=df, x="Method", y="NMI", order=order, color="lightsteelblue",
                fliersize=2, linewidth=0.8, ax=ax)
    ax.set_ylim(min(0.2, df["NMI"].min() * 0.95), 1.1); ax.set_xlabel("")
    ax.set_ylabel("Pairwise NMI among triplicates", fontsize=FS_YLAB)
    ax.set_xticklabels(order, rotation=90, ha="center", fontsize=FS_XTICK)
    ax.tick_params(axis="y", labelsize=FS_YTICK)

def sd_panel(ax, df, ct_order):
    order = df.groupby("Method")["Std"].mean().sort_values().index.tolist()   # most->least reproducible
    d = df.copy()
    d["CellType"] = pd.Categorical(d["CellType"], categories=ct_order, ordered=True)
    d["Method"] = pd.Categorical(d["Method"], categories=order, ordered=True)
    sns.barplot(data=d, x="Method", y="Std", color="lightgray", ci=None, ax=ax)   # seaborn 0.11.x
    sns.stripplot(data=d, x="Method", y="Std", hue="CellType", hue_order=ct_order,
                  order=order, dodge=True, size=5, alpha=0.85, palette="Dark2", ax=ax)
    ax.set_ylim(0, max(0.22, df["Std"].max() * 1.15)); ax.set_xlabel("")
    ax.set_ylabel("Std. dev. of ATAC cell type\nprediction across triplicates", fontsize=FS_YLAB)
    ax.set_xticklabels(order, rotation=90, ha="center", fontsize=FS_XTICK)
    ax.tick_params(axis="y", labelsize=FS_YTICK)
    ax.legend(loc="upper left", frameon=True, fontsize=FS_LEG, title="Cell type", title_fontsize=FS_LEGT)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.8))
nmi_panel(ax1, nmi_df)
sd_panel(ax2, sd_df, RMS_CT)
plt.tight_layout(w_pad=2.0)
fig.savefig(os.path.join(FIG, "fig4_combined.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig4_combined.png"), dpi=150, bbox_inches="tight")
print("wrote fig4_combined.pdf/png")
