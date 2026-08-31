#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
Per-cell view of MYOD1/MYOG ATAC vs RNA along the pseudotime RANK (equal spacing, since raw DPT pseudotime is
skewed). Shows every cell + a rolling-mean trend (window of cells, far finer than 15 bins), and quantifies the
lead with a binning-free CROSS-CORRELATION lag: the pseudotime-rank shift that best aligns the ATAC and RNA
trends. Positive => ATAC leads RNA (chromatin precedes transcription).
  python 01_plot_percell.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = "/path/to/multiomeBench/RMS/benchmark/atac_precede_rna"
STATE_COL = {"Mesoderm": "#2ca02c", "Myoblast": "#ff7f0e", "Myocyte": "#9467bd"}
W = 400   # rolling window in cells

pt = pd.read_csv(os.path.join(HERE, "pseudotime_table.csv"))          # barcode_base, cell_type, dpt_pseudotime, MYOD1/MYOG_rna
at = pd.read_csv(os.path.join(HERE, "region_atac_percell.csv"))       # barcode_base, MYOD1/MYOG_atac
df = pt.merge(at[["barcode_base", "MYOD1_atac", "MYOG_atac"]], on="barcode_base", how="inner")
df = df.sort_values("dpt_pseudotime").reset_index(drop=True)
N = len(df)
df["rank"] = np.arange(N) / (N - 1)                                   # even spacing along the ordering
x = df["rank"].to_numpy()


def roll(y):
    return pd.Series(y).rolling(W, center=True, min_periods=W // 4).mean().to_numpy()


def minmax(v):
    v = np.asarray(v, float); lo, hi = np.nanmin(v), np.nanmax(v)
    return (v - lo) / (hi - lo) if hi > lo else v * 0


def xcorr_lag(a, b):
    """rank-shift s (fraction of trajectory) that maximizes corr(ATAC(x), RNA(x+s)); positive => ATAC leads."""
    a, b = minmax(a), minmax(b)
    good = ~np.isnan(a) & ~np.isnan(b)
    a, b = a[good], b[good]
    n = len(a); best_s, best_r = 0, -2
    for shift in range(-n // 3, n // 3):
        if shift >= 0:
            aa, bb = a[:n - shift], b[shift:]          # compare ATAC(t) vs RNA(t+shift)
        else:
            aa, bb = a[-shift:], b[:n + shift]
        if len(aa) > n // 2:
            r = np.corrcoef(aa, bb)[0, 1]
            if r > best_r:
                best_r, best_s = r, shift
    return best_s / n, best_r                                          # lag as fraction of trajectory


fig, axes = plt.subplots(3, 2, figsize=(12, 11))
for j, g in enumerate(["MYOD1", "MYOG"]):
    atac_s, rna_s = roll(df[f"{g}_atac"]), roll(df[f"{g}_rna"])
    # row 0: ATAC per-cell + trend
    ax = axes[0, j]
    for s, c in STATE_COL.items():
        m = df["cell_type"].to_numpy() == s
        ax.scatter(x[m], df[f"{g}_atac"].to_numpy()[m], s=4, c=c, alpha=0.12, edgecolors="none")
    ax.plot(x, atac_s, color="#1f77b4", lw=2.6)
    ax.set_title(f"{g} — ATAC accessibility (per cell + rolling mean)"); ax.set_ylabel("frags / 10k")
    # row 1: RNA per-cell + trend
    ax = axes[1, j]
    for s, c in STATE_COL.items():
        m = df["cell_type"].to_numpy() == s
        ax.scatter(x[m], df[f"{g}_rna"].to_numpy()[m], s=4, c=c, alpha=0.12, edgecolors="none")
    ax.plot(x, rna_s, color="#d62728", lw=2.6)
    ax.set_title(f"{g} — RNA expression (per cell + rolling mean)"); ax.set_ylabel("counts")
    # row 2: normalized overlay + cross-correlation lag
    ax = axes[2, j]
    ax.plot(x, minmax(atac_s), color="#1f77b4", lw=2.6, label="ATAC (norm.)")
    ax.plot(x, minmax(rna_s),  color="#d62728", lw=2.6, label="RNA (norm.)")
    lag, r = xcorr_lag(atac_s, rna_s)
    ax.set_title(f"{g} — overlay | cross-corr lag = {lag:+.2f} of trajectory (peak r={r:.2f})")
    ax.set_ylabel("min–max normalized"); ax.legend(frameon=False, fontsize=9, loc="upper left")
    for ax in axes[:, j]:
        ax.set_xlim(0, 1); ax.set_xlabel("pseudotime rank (Mesoderm → Myocyte)")
    print(f"{g}: cross-correlation lag = {lag:+.3f} of trajectory (peak r={r:.2f})  "
          f"{'-> ATAC leads RNA' if lag > 0 else '-> no clear lead'}")

fig.suptitle("RMS: per-cell MYOD1/MYOG accessibility vs expression along pseudotime rank", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.98])
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"atac_precede_rna_percell.{ext}"), dpi=180)
print("\nwrote atac_precede_rna_percell.png/.pdf")
