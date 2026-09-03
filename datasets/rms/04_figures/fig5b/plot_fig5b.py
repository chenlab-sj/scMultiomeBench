#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
RMS Fig5B heatmap = Euclidean distance of each method's predicted-cell-type ATAC coverage track vs the
ground-truth (annotation) track, over the MYOD1 + FOXO1 regions, per cell type. Reconstructs the published
Fig5B from the per-position value CSVs and ADDS the 3 new methods (MaxFuse/MIDAS/scButterfly).

Old methods -> old/Mast607A/knn_test/predicted_ataclabel_{myod1,foxo1}value.csv (published tracks).
New methods -> fig5b/new_methods_{myod1,foxo1}value.csv (04_compute_new_pileups.R, same pipeline/params).
Fig5B set (matches Fig5A): scglue(multiome), scglue, Seurat(CCA), simba, scVI, scBridge, Portal, scJoint, BindSC
+ Random (scDART/Cobolt excluded). New 3 inserted by performance; Random last.  python plot_fig5b.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

G = "/path/to/multiomeBench"
OLD = f"{G}/results/rms/knn_test"
NEW = f"{G}/RMS/benchmark/fig5b"
CTS = ["Myocyte", "Myoblast", "Mesoderm"]                 # heatmap row order (top->bottom), matches published
REMAP = {"scglue.multiome.": "scglue(multiome)", "Seurat.CCA.": "Seurat(CCA)", "Bindsc": "BindSC", "random": "Random"}
# Fig5B method set (old) -- match Fig5A: BindSC kept, scDART + Cobolt excluded (Random still kept as baseline)
OLD_SET = ["scglue(multiome)", "scglue", "Seurat(CCA)", "simba", "scVI", "scBridge", "Portal", "scJoint", "BindSC"]
NEW_SET = ["MaxFuse", "MIDAS", "scButterfly"]

def dists(tf):
    old = pd.read_csv(f"{OLD}/predicted_ataclabel_{tf}value.csv"); old["pipeline"] = old["pipeline"].replace(REMAP)
    new = pd.read_csv(f"{NEW}/new_methods_{tf}value.csv");        new["pipeline"] = new["pipeline"].replace(REMAP)
    ann = old[old.pipeline == "annotation"].set_index("position")[CTS]
    rows = []
    for src in (old, new):
        for pl, g in src.groupby("pipeline"):
            if pl == "annotation":
                continue
            gg = g.set_index("position")[CTS].reindex(ann.index)
            for ct in CTS:
                rows.append({"pipeline": pl, "celltype": ct, "dist": float(np.sqrt(np.nansum((gg[ct] - ann[ct]) ** 2)))})
    return pd.DataFrame(rows)

d = {tf: dists(tf) for tf in ("myod1", "foxo1")}
both = pd.concat(d.values())

# order methods best->worst by mean distance over both TFs + all celltypes; Random pinned last
present = [m for m in OLD_SET + NEW_SET if m in set(both.pipeline)]
mean_d = both[both.pipeline.isin(present)].groupby("pipeline")["dist"].mean().sort_values()
order = mean_d.index.tolist() + (["Random"] if "Random" in set(both.pipeline) else [])

vmin, vmax = 10, 45
fig = plt.figure(figsize=(0.72 * len(order) + 1.6, 6.4))
gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.12)
for k, (tf, lab) in enumerate([("myod1", "MYOD1"), ("foxo1", "FOXO1")]):
    ax = fig.add_subplot(gs[k])
    mat = d[tf].pivot(index="celltype", columns="pipeline", values="dist").reindex(index=CTS, columns=order)
    im = ax.imshow(mat.values, cmap="Reds_r", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_yticks(range(len(CTS))); ax.set_yticklabels(CTS, fontsize=17)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order if k == 1 else [], rotation=45, ha="right", fontsize=17)
    ax.set_ylabel(lab, fontsize=19, rotation=270, labelpad=24); ax.yaxis.set_label_position("right")
    for s in ax.spines.values():
        s.set_visible(True)
    ax.set_xticks(np.arange(-.5, len(order), 1), minor=True); ax.set_yticks(np.arange(-.5, len(CTS), 1), minor=True)
    ax.grid(which="minor", color="black", linewidth=0.6); ax.tick_params(which="minor", length=0)

cax = fig.add_axes([0.30, 0.97, 0.4, 0.025])
cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[20, 40])
cb.ax.tick_params(labelsize=14); cb.set_label("Euclidean distance to ground truth", fontsize=16, labelpad=8)
cb.ax.xaxis.set_label_position("top"); cb.ax.xaxis.set_ticks_position("bottom")

fig.savefig(f"{NEW}/fig5b_heatmap.pdf", bbox_inches="tight")
fig.savefig(f"{NEW}/fig5b_heatmap.png", dpi=200, bbox_inches="tight")
print("method order:", order)
print("wrote fig5b_heatmap.pdf/png")
