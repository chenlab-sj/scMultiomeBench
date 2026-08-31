#!/usr/bin/env python3
"""
Combined k-sensitivity figure: FigS1A = pbmc3k, FigS1B = pbmc10k (both = KNN ATAC-label-prediction accuracy vs k,
box = per-cell-type accuracy). Faceted 2 rows x 1 col, ONE shared method-color legend on top, ONE shared x-axis
(k value) at the bottom. 23 methods each (18 published + 5 new). Same palette as plot_figS1a.py / plot_ksens.py.
  python plot_figS1_ksens_combined.py
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

G = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
OUT = os.path.dirname(os.path.abspath(__file__))
KS = [5, 10, 20, 40, 80]

METHOD_COLORS = {
    "scglue": "#1f77b4", "scglue(multiome)": "#aec7e8", "Conos": "#ff7f0e", "Portal": "#ffbb78",
    "scBridge": "#2ca02c", "simba": "#98df8a", "Seurat(CCA)": "#d62728", "scJoint": "#ff9896",
    "scVI": "#9467bd", "MinNet": "#c5b0d5", "Cobolt": "#c49c94", "BindSC": "#e377c2",
    "LIGER": "#f7b6d2", "MultiMAP": "#7f7f7f", "scMoMaT": "#c7c7c7", "Unioncom": "#dbdb8d",
    "scDART": "#393b79", "Seurat(WNN)": "#9edae5",
    "MaxFuse": "#17becf", "MIDAS": "#bcbd22", "scButterfly": "#8c564b",
    "MIRA": "#843c39", "Multigrate": "#7b4173",
}
ORDER = ["BindSC", "Cobolt", "Conos", "LIGER", "MinNet", "MultiMAP", "Portal", "Seurat(CCA)",
         "Seurat(WNN)", "Unioncom", "scBridge", "scDART", "scJoint", "scMoMaT", "scVI", "scglue",
         "scglue(multiome)", "simba", "MaxFuse", "MIDAS", "scButterfly", "MIRA", "Multigrate"]

def load_long(accu_sum, new_long):
    df = pd.read_csv(accu_sum)
    cts = [c for c in df.columns if c not in ("pipeline", "overall", "accu_average", "k value")]
    long = df.melt(id_vars=["pipeline", "k value"], value_vars=cts, var_name="celltype", value_name="accuracy")
    long = long[long["k value"].isin(KS)]
    if os.path.exists(new_long):
        nl = pd.read_csv(new_long)
        nl = nl[nl["type"] != "overall"].rename(columns={"type": "celltype"})
        long = pd.concat([long, nl[["pipeline", "k value", "celltype", "accuracy"]]], ignore_index=True)
    return long

pbmc3k = load_long(f"{G}/pbmc/pbmc3k/benchmark/ksensitivity/kNN_celltype_accu_sum.csv",
                   f"{G}/pbmc/pbmc3k/benchmark/ksensitivity/new_methods_ktest_long.csv")
pbmc10k = load_long(f"{G}/pbmc/pbmc10k/benchmark/old/knn_test/kNN_celltype_accu_sum.csv",
                    f"{G}/pbmc/pbmc10k/benchmark/new_methods_ktest_long.csv")

present = [m for m in ORDER if m in (set(pbmc3k["pipeline"]) | set(pbmc10k["pipeline"]))]
palette = {m: METHOD_COLORS[m] for m in present}

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(0.42 * len(present) + 2, 4.4))
for ax, data, name in [(axes[0], pbmc3k, "PBMC 3k"), (axes[1], pbmc10k, "PBMC 10k")]:
    d = data.copy()
    d["pipeline"] = pd.Categorical(d["pipeline"], categories=present, ordered=True)
    d["k value"] = pd.Categorical(d["k value"], categories=KS, ordered=True)
    sns.boxplot(data=d, x="k value", y="accuracy", hue="pipeline", hue_order=present,
                palette=palette, ax=ax, linewidth=0.5, fliersize=2.2,
                flierprops=dict(marker="o", markerfacecolor="none", markeredgewidth=0.5))
    ax.set_ylabel("Accuracy", fontsize=13)
    ax.set_ylim(-0.03, 1.03)
    ax.tick_params(labelsize=12)
    ax.legend_.remove()
    ax.text(-0.14, 0.5, name, transform=ax.transAxes, rotation=90, va="center", ha="center",
            fontsize=13, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="0.9", edgecolor="0.7"))  # left facet strip, vertical
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

axes[0].set_xlabel("")                             # seaborn auto-sets it; drop -> single shared x-axis
axes[1].set_xlabel("k value", fontsize=13)         # only the bottom panel keeps the x-axis label

# one shared legend, on top
handles = [mpatches.Patch(facecolor=palette[m], edgecolor="0.3", label=m) for m in present]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0),
           ncol=8, fontsize=11, frameon=True, handlelength=1.1, columnspacing=1.0, handletextpad=0.4)
fig.tight_layout(rect=(0, 0, 1, 0.77))     # reserve top ~20% for the (now 3-row) legend so it clears the panels

fig.savefig(f"{OUT}/figS1_ksens_combined.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/figS1_ksens_combined.png", dpi=200, bbox_inches="tight")
print(f"methods ({len(present)}):", present)
print("wrote figS1_ksens_combined.pdf/png")
