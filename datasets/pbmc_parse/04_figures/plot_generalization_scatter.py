#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""R1#1 generalization -- pbmc cross-PLATFORM score-comparison scatter.
x = integration performance score on Parse RNA x 10X pbmc3k ATAC (cross-platform, truly independent)
y = integration performance score on 10X pbmc3k sc-multiome
Each retained method = one point (scDART + Cobolt excluded, matching the main analysis), coloured by
rank change. A y=x diagonal is drawn: ABOVE the line = multiome score not matched cross-platform =
weak generalizer (scButterfly, scJoint); on/below = generalizes at least as well. NB the two composites
differ (multiome is paired -> includes ari/ami/ks.statistic; cross-platform unpaired drops them), a mild
systematic offset -> read the diagonal as a relative guide, not a strict equal-difficulty claim."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

G = "/path/to/multiomeBench"
def load(p):
    d = pd.read_csv(p); d.columns = [c.strip('"') for c in d.columns]
    d["method"] = d["method"].astype(str).str.strip('"')
    return d.set_index("method")["score"]
mult  = load(f"{G}/pbmc/pbmc3k/benchmark/fig2b/fig2b_matrix.csv")
parse = load(f"{G}/pbmc_parse/benchmark/metrics/metrics_matrix.csv")

EXCLUDE = {"scDART", "Cobolt"}
methods = [m for m in mult.index if m in parse.index and m not in EXCLUDE]
x = parse[methods].astype(float); y = mult[methods].astype(float)   # x = cross-platform, y = multiome
rho = spearmanr(x, y).correlation

# rank change (1 = best) -> colour. drop = cross-platform rank - multiomics rank.
# Negative values indicate a rank increase relative to multiomics; positive values indicate a rank decrease.
rank_mult = mult[methods].rank(ascending=False); rank_parse = parse[methods].rank(ascending=False)
drop = rank_parse - rank_mult
# RdBu: red = negative drop (rank increase), blue = positive drop (rank decrease), white = unchanged.
dmax = max(1, int(drop.abs().max())); norm = plt.Normalize(-dmax, dmax); cmap = plt.cm.RdBu

fig, ax = plt.subplots(figsize=(6.6, 6.4))
lo, hi = 0.5, 0.85
# y=x reference: both axes are the [0,1] integration performance score, so equal-score is the natural
# guide. Points ABOVE the line = higher multiome than cross-platform score = weaker cross-platform
# generalization (scButterfly, scJoint); on/below = generalizes at least as well. NB the two composites
# differ slightly (multiome is paired -> includes ari/ami/ks.statistic; cross-platform drops them), so
# read the line as a relative guide, not a strict equal-difficulty claim.
ax.plot([lo, hi], [lo, hi], ls="--", color="#888888", lw=1.0, zorder=1)
ax.scatter(x, y, s=95, c=[cmap(norm(drop[m])) for m in methods], edgecolor="black", lw=0.5, zorder=3)
# manual label placement (dx, dy, ha) to de-overlap the x0.70-0.76 x y0.62-0.72 cluster
LP = {
    "scButterfly": (0.008, 0.0, "left"),  "scJoint": (0.008, -0.003, "left"),
    "MaxFuse": (-0.008, -0.004, "right"), "MIDAS": (-0.008, 0.002, "right"),
    "Seurat(CCA)": (0.0, 0.011, "center"), "simba": (0.008, 0.002, "left"),
    "BindSC": (0.0, -0.012, "center"),    "scVI": (-0.007, 0.009, "right"),
    "scBridge": (0.0, -0.012, "center"),  "Portal": (0.008, 0.001, "left"),
    "scglue": (0.008, -0.004, "left"),    "scglue(multiome)": (-0.008, 0.004, "right"),
}
for m in methods:
    dx, dy, ha = LP.get(m, (0.008, 0, "left"))
    ax.annotate(m, (x[m], y[m]), xytext=(x[m] + dx, y[m] + dy), fontsize=8, va="center", ha=ha, zorder=4)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel("Integration performance score\nParse pbmc × 10X pbmc3k scATAC (cross-platform)", fontsize=10)
ax.set_ylabel("Integration performance score\n10X pbmc3k sc-multiome", fontsize=10)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04); cb.set_label("rank change vs sc-multiomics\nred: increase; blue: decrease", fontsize=9)
fig.tight_layout(); fig.savefig(f"{G}/pbmc_parse/benchmark/r1_pbmc_scatter.png", dpi=170, bbox_inches="tight")
print(f"Spearman rho = {rho:.3f}  over {len(methods)} methods")
print("scores (cross-platform x -> multiome y):")
for m in x.sort_values(ascending=False).index:
    print(f"  {m:18s} {x[m]:.3f} -> {y[m]:.3f}  (rankΔ {int(drop[m]):+d})")
print("wrote r1_pbmc_scatter.png")
