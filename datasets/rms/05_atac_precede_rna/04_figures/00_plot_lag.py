#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
W7 step 3: per-pseudotime-bin ATAC accessibility vs RNA expression at MYOD1 & MYOG. Bins are equal-count
(quantile) so the x-axis is the pseudotime RANK (Mesoderm -> Myocyte). If chromatin priming precedes
transcription, the ATAC curve rises at an earlier rank than RNA. Lag = (rank where RNA reaches half-max)
- (rank where ATAC reaches half-max); positive => ATAC first.
  python 00_plot_lag.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/path/to/multiomeBench/RMS/benchmark/atac_precede_rna"
STATE_COL = {"Mesoderm": "#2ca02c", "Myoblast": "#ff7f0e", "Myocyte": "#9467bd"}

pt = pd.read_csv(os.path.join(HERE, "pseudotime_table.csv"), index_col=0)     # per-cell RNA + bin + cell_type
at = pd.read_csv(os.path.join(HERE, "region_atac_percell.csv"))              # per-cell ATAC + bin

rna_bin = pt.groupby("bin")[["MYOD1_rna", "MYOG_rna"]].mean()
atac_bin = at.groupby("bin")[["MYOD1_atac", "MYOG_atac"]].mean()
dom = pt.groupby("bin")["cell_type"].agg(lambda s: s.value_counts().idxmax())   # dominant state per bin
bins = sorted(set(rna_bin.index) & set(atac_bin.index))
x = np.array(bins, float)                                                    # pseudotime RANK (equal-count bins)


def minmax(v):
    v = np.asarray(v, float); lo, hi = np.nanmin(v), np.nanmax(v)
    return (v - lo) / (hi - lo) if hi > lo else v * 0


def halfmax_x(xv, yv):
    y = minmax(yv)
    for i in range(1, len(y)):
        if y[i - 1] < 0.5 <= y[i]:
            t = (0.5 - y[i - 1]) / (y[i] - y[i - 1])
            return xv[i - 1] + t * (xv[i] - xv[i - 1])
    return xv[int(np.nanargmax(y))]


fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
for ax, g in zip(axes, ["MYOD1", "MYOG"]):
    # state-composition strip (background shading by dominant state along the ranking)
    for b in bins:
        ax.axvspan(b - 0.5, b + 0.5, color=STATE_COL.get(dom.get(b), "#cccccc"), alpha=0.12, lw=0)
    rna = rna_bin.loc[bins, f"{g}_rna"].to_numpy()
    atac = atac_bin.loc[bins, f"{g}_atac"].to_numpy()
    ax.plot(x, minmax(atac), "-s", color="#1f77b4", label="ATAC accessibility")
    ax.plot(x, minmax(rna), "-o", color="#d62728", label="RNA expression")
    xa, xr = halfmax_x(x, atac), halfmax_x(x, rna)
    ax.axvline(xa, color="#1f77b4", ls=":", lw=1.2); ax.axvline(xr, color="#d62728", ls=":", lw=1.2)
    lag = xr - xa
    ax.set_title(f"{g}   (half-max lag RNA−ATAC = {lag:+.2f} bins)")
    ax.set_xlabel("pseudotime rank (equal-count bins)  →")
    ax.set_ylabel("min–max normalized signal")
    ax.set_xlim(min(bins) - 0.5, max(bins) + 0.5)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    print(f"{g}: ATAC half-max @ bin {xa:.2f}, RNA half-max @ bin {xr:.2f}, lag(RNA-ATAC) = {lag:+.2f} bins "
          f"({'ATAC precedes RNA' if lag > 0 else 'no lead / RNA first'})")

# state legend
handles = [plt.Line2D([0], [0], marker="s", ls="", color=c, alpha=0.4, label=s) for s, c in STATE_COL.items()]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9,
           title="dominant state per bin", bbox_to_anchor=(0.5, -0.02))
fig.suptitle("RMS myogenic differentiation: chromatin accessibility rises before transcription", fontsize=12)
fig.tight_layout(rect=[0, 0.05, 1, 1])
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"atac_precede_rna_lag.{ext}"), dpi=200, bbox_inches="tight")
print("\nwrote atac_precede_rna_lag.png/.pdf")
