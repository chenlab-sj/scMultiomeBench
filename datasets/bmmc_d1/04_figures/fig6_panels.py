#!/usr/bin/env python3
"""
Export Fig6 panels A, B, C as SEPARATE self-contained vector PDFs, for manual assembly/adjustment.
Reuses the exact functions behind the combined figure (make_fig6_combined + make_r1_brca_pbmc_composite),
so the panels are byte-for-byte the same plots, just split out. A carries its Setting legend + centred
y-label; B and C each carry their own rank-change colorbar. Also drops a .png preview beside each PDF.

    python fig6_panels.py
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_fig6_combined as fc                       # shared: bmmc_scores, draw_bmmc, cmap, FS_*, FACETS
sys.path.insert(0, os.path.join(fc.G, "BRCA", "benchmark"))
from make_r1_brca_pbmc_composite import brca_data, pbmc_data, plot_brca, plot_pbmc  # noqa: E402

cmap = fc.cmap
CBAR_LABEL = "Rank change vs sc-multiomics\nred: rank decrease (better)\nblue: rank increase (worse)"

# shared rank-change norm (identical to the combined figure)
brca_conds, brca_scores, brca_change = brca_data()
pbmc_methods, pbmc_x, pbmc_y, pbmc_change = pbmc_data()
dmax = max(1, int(max(brca_change.abs().values.max(), pbmc_change.abs().max())))   # .values -> scalar (brca_change is a DataFrame)
norm = plt.Normalize(-dmax, dmax)


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"{stem}.{ext}"), dpi=250, bbox_inches="tight")
    print("wrote", stem + ".pdf/.png")


def add_colorbar(fig, rect):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
    cb.set_label(CBAR_LABEL, fontsize=fc.FS_CBAR)
    cb.ax.tick_params(labelsize=fc.FS_TICK)


# ---------- Panel A: BMMC 4-facet benchmark ----------
tab, order = fc.bmmc_scores()
figA = plt.figure(figsize=(8.5, 10.5))
gsA = figA.add_gridspec(4, 1, hspace=0.30)
axesA = [figA.add_subplot(gsA[i, 0]) for i in range(4)]
fc.draw_bmmc(axesA, tab, order)
p0, p3 = axesA[0].get_position(), axesA[3].get_position()
figA.text(p0.x0 - 0.095, p0.y1 + 0.012, "A", ha="left", va="bottom",
          fontsize=fc.FS_LETTER, fontweight="bold", fontfamily="Arial")
figA.text(p0.x0 - 0.072, (p0.y1 + p3.y0) / 2, "Benchmark metrics value",
          rotation=90, va="center", ha="center", fontsize=fc.FS_YLAB)
save(figA, "fig6_panelA")

# ---------- Panel B: BRCA cross library/donor slope ----------
figB = plt.figure(figsize=(10.5, 6.2))                       # wide enough to separate the 3 condition labels
axB = figB.add_axes([0.12, 0.19, 0.71, 0.73])
plot_brca(axB, brca_conds, brca_scores, brca_change, cmap, norm, panel_label=None,
          label_fs=fc.FS_LABEL, tick_fs=fc.FS_TICK, axis_fs=fc.FS_AXIS, minsep_frac=0.058, fit_labels=True)
add_colorbar(figB, [0.87, 0.32, 0.018, 0.44])
figB.text(0.02, 0.95, "B", fontsize=fc.FS_LETTER, fontweight="bold", fontfamily="Arial", va="bottom")
save(figB, "fig6_panelB")

# ---------- Panel C: PBMC cross-platform scatter ----------
figC = plt.figure(figsize=(7.6, 6.6))
axC = figC.add_axes([0.15, 0.14, 0.64, 0.78])
plot_pbmc(axC, pbmc_methods, pbmc_x, pbmc_y, pbmc_change, cmap, norm, panel_label=None,
          label_fs=15, tick_fs=fc.FS_TICK, axis_fs=fc.FS_AXIS)   # larger labels; adjustText avoids overlap
add_colorbar(figC, [0.85, 0.28, 0.022, 0.46])
figC.text(0.02, 0.95, "C", fontsize=fc.FS_LETTER, fontweight="bold", fontfamily="Arial", va="bottom")
save(figC, "fig6_panelC")
