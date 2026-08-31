#!/usr/bin/env python3
"""
R2.1 deliverable: does the choice of ATAC->gene-activity estimator change the benchmark?
Compares each gene-activity-dependent method's pbmc3k composite score under the two estimators:
  Signac GeneActivity (the main fig2b baseline)  vs  ArchR GeneScoreMatrix (this experiment).
Both scored with the IDENTICAL paired composite (plot_metrics_matrix.R), so scores are directly comparable.

HEADLINE claim (Spearman + scatter) is computed over the 5 methods where ONLY the gene-activity
estimator changed -- a clean, like-for-like swap:
    Seurat(CCA), BindSC, MaxFuse, scBridge, Portal.
scJoint is kept in the CSV/table but is NOT plotted and is excluded from rho, for two confounds that the
others don't have:
  (1) its original Signac prep was overwritten, so the ArchR run is a pipeline RECONSTRUCTION, not the
      same pipeline with only the estimator swapped; and
  (2) scJoint uniquely requires COUNT-scale input, so the continuous ArchR scores had to be rounded to
      integer counts -- a method-specific sensitivity, not a property of the benchmark.

RESULT (do not mis-state it): the estimator DOES matter for some methods. rho=0.70 (p=0.19), mean|d|=0.045,
max|d|=0.110; Seurat(CCA)/BindSC/MaxFuse are ~unchanged (|d|<=0.017) while Portal (-0.110, rank 2->4) and
scBridge (-0.089) drop substantially. The benchmark's HEADLINE is nonetheless unaffected, because the
top-ranked methods (scGLUE, scVI) are peak-based and never consume gene activity.

Reads ../fig2b/sum_metrics_clean.csv (Signac) + metrics/sum_metrics_clean.csv (ArchR). Run AFTER
metrics/run_all.sh.  python3.9 compare_geneactivity.py   (python3 = 3.11 here has a broken scipy ABI)
Outputs: compare_geneactivity.csv + compare_geneactivity.png/pdf
"""
import math
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

HERE = os.path.dirname(os.path.abspath(__file__))
CLEAN_METHODS   = ["Seurat(CCA)", "BindSC", "MaxFuse", "scBridge", "Portal"]   # clean estimator swap -> headline
FLAGGED_METHODS = ["scJoint"]                                                  # confounded -> in table, not plotted

# Reference points ONLY (drawn on y=x, grey, never in rho and never written to the CSV): methods VERIFIED
# against the repo code not to consume an ATAC->gene-activity matrix. For these the swap is a literal no-op
# -- the input they read is byte-identical -- so ArchR == Signac is a logical identity, not an assumption or
# a measurement. They were NOT re-run, and the figure says so.
REFERENCE_METHODS = ["scglue(multiome)", "scglue", "scVI", "scDART", "scButterfly",
                     "MIDAS", "simba", "Cobolt", "MIRA", "Seurat(WNN)"]
# DELIBERATELY NOT PLOTTED -- gene-activity dependence UNRESOLVED (run code lives outside the repo).
# LIGER is very likely AFFECTED (its canonical scRNA+scATAC mode consumes a gene-activity matrix) and Conos
# likely so. Drawing any of these on y=x would assert "no effect" for methods we never tested and have
# reason to think ARE affected -- the one claim here a reviewer could falsify. Do not add them.
UNRESOLVED = ["LIGER", "Conos", "MultiMAP", "Unioncom", "scMoMaT", "Multigrate", "MinNet"]
COL = "Overall Performance"


def load(path, tag):
    df = pd.read_csv(path)
    df.columns = [c.strip('"') for c in df.columns]
    df["method"] = df["method"].astype(str).str.strip('"')
    return df.set_index("method")[COL].rename(tag)


sig = load(os.path.join(HERE, "..", "fig2b", "sum_metrics_clean.csv"), "Signac")
arc = load(os.path.join(HERE, "metrics", "sum_metrics_clean.csv"), "ArchR")

cmp = pd.concat([sig, arc], axis=1).reindex(CLEAN_METHODS + FLAGGED_METHODS).dropna()
cmp["delta(ArchR-Signac)"] = (cmp["ArchR"] - cmp["Signac"]).round(3)
cmp["in_headline"] = [m in CLEAN_METHODS for m in cmp.index]

# ranks + Spearman over the CLEAN methods only (the like-for-like swap)
clean = cmp[cmp["in_headline"]].copy()
clean["Signac_rank"] = clean["Signac"].rank(ascending=False).astype(int)
clean["ArchR_rank"]  = clean["ArchR"].rank(ascending=False).astype(int)
cmp = cmp.join(clean[["Signac_rank", "ArchR_rank"]])
cmp = cmp.round(3).sort_values(["in_headline", "Signac"], ascending=[False, False])
cmp.to_csv(os.path.join(HERE, "compare_geneactivity.csv"))

rho = p = None
if len(clean) >= 3:
    try:
        from scipy.stats import spearmanr
        rho, p = spearmanr(clean["Signac"], clean["ArchR"])
    except Exception:
        rho = clean[["Signac", "ArchR"]].corr(method="spearman").iloc[0, 1]

print(f"headline methods (clean swap): {len(clean)}  ->  {list(clean.index)}")
if rho is not None:
    print(f"Spearman rank (Signac vs ArchR): rho={rho:.2f}" + (f" (p={p:.2g})" if p is not None else ""))
print(f"mean |delta| (clean): {clean['delta(ArchR-Signac)'].abs().mean():.3f}   "
      f"max |delta| (clean): {clean['delta(ArchR-Signac)'].abs().max():.3f}")
flagged = cmp[~cmp["in_headline"]]
if len(flagged):
    print("flagged (excluded from rho):")
    print(flagged[["Signac", "ArchR", "delta(ArchR-Signac)"]].to_string())
print("\nfull table:")
print(cmp.to_string())

# ---- figure: Signac vs ArchR composite; y=x = unchanged by the gene-activity estimator ----
# scJoint is NOT plotted: it is not a clean estimator swap (reconstructed pipeline + count-scale rounding),
# so it does not belong on a controlled-comparison scatter, and its collapse (0.676 -> 0.101) would compress
# the axes ~5x and hide the effect we are actually showing. It stays in the CSV / table / response text.
# Peak-based methods (scGLUE, scVI, ...) are likewise NOT plotted: they never consume a gene-activity matrix,
# so re-running them is a no-op. Plotting them at Signac==ArchR would mean drawing ~15 unmeasured points
# exactly on y=x -- an assertion, not a measurement -- and LIGER/Conos/MultiMap were never resolved (LIGER is
# very likely AFFECTED), so "all other methods" is not a safe set. That scope claim belongs in prose.
ref = sig.reindex(REFERENCE_METHODS).dropna()                     # Signac composite; ArchR == Signac by identity
missing_ref = [m for m in REFERENCE_METHODS if m not in ref.index]
if missing_ref:
    print(f"NOTE reference methods not found in the baseline (skipped): {missing_ref}")
print(f"reference (not re-run, no gene-activity input): {len(ref)} -> {list(ref.index)}")
print(f"NOT plotted (dependence unresolved; LIGER/Conos likely affected): {UNRESOLVED}")

# Type sizes for the manuscript figure, kept in one place so they scale together. If you raise these
# further, raise `figsize` with them -- the label placer needs blank space to work with, and it will
# print a NOTE for any label it cannot fit.
FS_MEASURED = 18   # labels on the re-run (blue) points
FS_REF      = 17   # labels on the not-re-run (open) points
FS_AXIS      = 21
FS_TICK      = 17
FS_TITLE     = 22
FS_LEGEND    = 14

# Label offsets are (radius, angle) pairs tried NEAREST-FIRST, so a label stays as close to its own
# marker as it can -- a label shoved far away is ambiguous, which is its own kind of unreadable.
# Angles are ordered perpendicular-to-the-diagonal first (up-left, down-right): every reference point
# sits ON y=x, so that is where the free space is. Anything pushed past _LEADER_MIN gets a thin
# leader line back to its marker, so a displaced label is never ambiguous.
# This replaces the old fixed `-8 - step*9.5` stack, which compared only scores and so still let
# simba/MIDAS and scDART/scVI/scglue collide while pushing scglue(multiome) off the right edge.
_ANGLES = [135, 315, 180, 0, 90, 270, 225, 45]
_RADII = [15, 23, 33, 45, 59, 75, 93]
_LEADER_MIN = 30
# scDART/scVI/scglue are 0.0017-0.0028 apart, i.e. their markers physically overlap on the page
# (~11 and ~19 px). For a point with a neighbour that close, ANY nearby label is ambiguous, so those
# start further out, fan out over more angles, and always get a leader line -- the reader can then
# trace each name to its own circle. Keep _CROWD_PX tight: at 45 px it caught pairs that are merely
# near each other (simba/MIDAS, Seurat(CCA)/scButterfly) and pushed every label out until the
# placement ran out of room.
_CROWD_PX = 26
_RADII_CROWDED = [34, 44, 56, 70, 86]
_ANGLES_ALL = list(range(0, 360, 30))


def _ang_dist(a, b):
    d = abs((a - b) % 360)
    return min(d, 360 - d)


def place_labels(ax, items, occupied, fontsize, color="black", markers=()):
    """Put each label at the nearest offset that collides with nothing and stays inside the axes.

    Measures the ACTUAL rendered text box, so labels keep clear of the markers, the legend, the y=x
    line and each other. Every accepted box is appended to `occupied`, so later labels avoid earlier
    ones -- call with the points that matter most first, they get first pick of the close-in spots.
    `markers` is the display-space position of every plotted point, used to spot crowded neighbours.
    """
    renderer = ax.figure.canvas.get_renderer()
    for text, x, y in items:
        px, py = ax.transData.transform((x, y))
        crowded = any(0.5 < math.hypot(mx - px, my - py) < _CROWD_PX for mx, my in markers)
        cand = [(r_, a_) for r_ in _RADII for a_ in _ANGLES]
        if crowded:
            # Send each label RADIALLY OUTWARD from the centre of its cluster. Points on one side of
            # the cluster get labels on that same side, so the leader lines fan out instead of
            # crossing each other -- crossed leaders are worse than none, they mis-assign the name.
            near = [(mx, my) for mx, my in markers
                    if math.hypot(mx - px, my - py) < _CROWD_PX * 2] or [(px, py)]
            cx = sum(m[0] for m in near) / len(near)
            cy = sum(m[1] for m in near) / len(near)
            # A point sitting at the cluster centre has no meaningful outward direction -> send it
            # perpendicular to y=x, where there is free space.
            base = (math.degrees(math.atan2(py - cy, px - cx))
                    if math.hypot(px - cx, py - cy) > 3 else 135.0)
            angles = sorted(_ANGLES_ALL, key=lambda a_: _ang_dist(a_, base))
            # fanned-out spots first, falling back to the ordinary close-in ones
            cand = [(r_, a_) for r_ in _RADII_CROWDED for a_ in angles] + cand
        for r, a in cand:
            dx, dy = r * math.cos(math.radians(a)), r * math.sin(math.radians(a))
            ha = "left" if dx > 2 else "right" if dx < -2 else "center"
            va = "bottom" if dy > 2 else "top" if dy < -2 else "center"
            # Measure the TEXT only. Attaching arrowprops here would make get_window_extent return
            # the text box unioned with the leader line, i.e. the whole diagonal span back to the
            # marker -- which then blocks every neighbouring label. The leader is added afterwards.
            ann = ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points",
                              fontsize=fontsize, color=color, ha=ha, va=va, zorder=5)
            bb = ann.get_window_extent(renderer=renderer).expanded(1.05, 1.25)
            inside = (ax.bbox.containsx(bb.x0) and ax.bbox.containsx(bb.x1)
                      and ax.bbox.containsy(bb.y0) and ax.bbox.containsy(bb.y1))
            if inside and not any(bb.overlaps(o) for o in occupied):
                occupied.append(bb)
                if crowded or r >= _LEADER_MIN:
                    ax.annotate("", xy=(x, y), xycoords="data",
                                xytext=(dx, dy), textcoords="offset points", zorder=4,
                                arrowprops=dict(arrowstyle="-", color=color, lw=0.7,
                                                shrinkA=3, shrinkB=6))
                break
            ann.remove()
        else:   # nothing fit -- keep the label rather than drop it, and say so
            ann = ax.annotate(text, (x, y), xytext=(15, 0), textcoords="offset points",
                              fontsize=fontsize, color=color, ha="left", va="center", zorder=5)
            occupied.append(ann.get_window_extent(renderer=renderer))
            print(f"NOTE: no clear spot for '{text}' -- placed at the default offset")


fig, ax = plt.subplots(figsize=(11.5, 11.0))
allv = pd.concat([clean["Signac"], clean["ArchR"], ref])          # limits over everything drawn
lo, hi = allv.min(), allv.max()
pad = (hi - lo) * 0.12 + 1e-6; lims = (lo - pad, hi + pad)
ax.plot(lims, lims, "--", color="black", lw=1.6, zorder=1)

# reference: NOT re-run -- drawn on y=x because the swap cannot change them (they never read the matrix).
# Open markers (not colour) are what separate these from the measured points, so they read correctly in
# black and in greyscale print.
ax.scatter(ref, ref, s=95, facecolors="none", edgecolors="black", linewidths=1.7, zorder=2,
           label=f"not re-run — no gene-activity input, unchanged by construction (n={len(ref)})")
# NB the legend is now the ONLY on-figure statement that the open points were not re-run (the explanatory
# footnote was removed for the manuscript version) -- do not shorten that label. Everything else (which
# methods are excluded and why) must live in the figure CAPTION; see R2.1_response_draft.md.

# measured: the clean estimator swap
ax.scatter(clean["Signac"], clean["ArchR"], s=185, color="#4C72B0", zorder=3,
           label=f"re-run with ArchR — gene-activity–dependent (n={len(clean)})")

ax.set_xlim(lims); ax.set_ylim(lims); ax.set_aspect("equal", "box")
# "Signac vs ArchR" is already on the axes -> not repeated in the title; rho lives in the caption/text.
ax.set_xlabel("Integration performance score — Signac GeneActivity", fontsize=FS_AXIS)
ax.set_ylabel("Integration performance score — ArchR GeneScoreMatrix", fontsize=FS_AXIS)
ax.set_title("Gene-activity estimator sensitivity (pbmc3k)", fontsize=FS_TITLE)
ax.tick_params(labelsize=FS_TICK)
leg = ax.legend(loc="upper left", fontsize=FS_LEGEND, frameon=False)
fig.tight_layout()

# Labels go on last, once the axes geometry is final, so the measured text boxes are exact.
fig.canvas.draw()
_renderer = fig.canvas.get_renderer()
occupied = [leg.get_window_extent(renderer=_renderer)]
for _x, _y in list(zip(ref.values, ref.values)) + list(zip(clean["Signac"], clean["ArchR"])):
    _px, _py = ax.transData.transform((_x, _y))
    occupied.append(Bbox.from_bounds(_px - 11, _py - 11, 22, 22))   # keep text off the markers
# ...and off the y=x line itself: sample it and block a thin band, otherwise labels for the reference
# points (which all sit on the line) get struck through by it.
for _i in range(161):
    _t = lims[0] + (lims[1] - lims[0]) * _i / 160.0
    _px, _py = ax.transData.transform((_t, _t))
    occupied.append(Bbox.from_bounds(_px - 4, _py - 4, 8, 8))
_markers = [ax.transData.transform((_x, _y))
            for _x, _y in list(zip(ref.values, ref.values)) + list(zip(clean["Signac"], clean["ArchR"]))]
place_labels(ax, [(m, r["Signac"], r["ArchR"]) for m, r in clean.iterrows()], occupied,
             FS_MEASURED, markers=_markers)
place_labels(ax, [(m, v, v) for m, v in ref.sort_values().items()], occupied,
             FS_REF, markers=_markers)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"compare_geneactivity.{ext}"), dpi=300)
print("\nwrote compare_geneactivity.csv + compare_geneactivity.png/pdf")
