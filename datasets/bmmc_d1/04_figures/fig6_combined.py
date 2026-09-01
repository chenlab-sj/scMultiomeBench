#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
NEW Figure 6 = BMMC_d1 benchmark (panel A) + two metric-fair generalization scatters that SHARE one
rank-change colour scale:
  B: BMMC same-donor cross-site  (s1d1 scRNA + s4d1 scATAC  vs  s1d1 sc-multiome, technical batch)
  C: PBMC cross-platform         (Parse scRNA + 10X scATAC  vs  10X sc-multiome)

Layout:
    +----------------------+   +--------------------+  +--+
    |                      |   |  B  BMMC crosssite |  |  |
    |  A  BMMC benchmark   |   +--------------------+  |cb|   (single shared colourbar)
    |     (4 metric facets)|   |  C  PBMC scatter   |  |  |
    +----------------------+   +--------------------+  +--+

Both scatters are scored on the SAME 6 metrics (load_score recomputes the unpaired composite, so the
sc-multiome baseline is NOT credited for the same-cell metrics ks.statistic/ari/ami it uniquely has).
The BRCA slope was moved to the supplement. B and C reuse plot_pbmc; the shared norm spans the larger of
the two panels' |rank change|. A's y-label is centred over the 4 facets; the diverging colormap has no
white middle.  python fig6_combined.py
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 14,
})

G = "/path/to/multiomeBench"
sys.path.insert(0, os.path.join(G, "BRCA", "benchmark"))
from make_r1_brca_pbmc_composite import bmmc_crosssite_data, pbmc_data, plot_pbmc  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
D    = os.path.normpath(os.path.join(HERE, ".."))
PUB  = "/path/to/multiomeBench/results/bmmc_d1/published_reference/sum_metrics_clean.csv"

# ---- font sizes (enlarged) ----
FS_TICK, FS_FACET, FS_XMETH, FS_YLAB, FS_LEG, FS_AXIS, FS_LABEL, FS_CBAR, FS_LETTER = 17, 17, 17, 20, 17, 19, 16, 17, 23

# ---------- BMMC panel config ----------
BATCH = {"scglue(multiome)": "scglue(multiome,batch)", "scVI": "scVI(batch)", "scglue": "scglue(batch)",
         "BindSC": "BindSC(batch)", "MIDAS": "MIDAS(batch)"}          # scJoint(batch) excluded
SELECTED = ["scglue(multiome)", "scVI", "scglue", "scBridge", "simba", "Seurat(CCA)", "BindSC",
            "Portal", "scJoint", "MaxFuse", "MIDAS", "scButterfly"]
COLORS = {"scglue": "#1f77b4", "scglue(multiome)": "#aec7e8", "Portal": "#ffbb78", "scBridge": "#2ca02c",
          "simba": "#98df8a", "Seurat(CCA)": "#d62728", "scJoint": "#ff9896", "scVI": "#9467bd",
          "BindSC": "#e377c2", "MaxFuse": "#17becf", "MIDAS": "#bcbd22", "scButterfly": "#7f7f7f"}
FACETS = ["Biological diversity\npreservation", "Batch effects\ncorrection",
          "Omics gaps\nreduction", "Optimal alignment\nbetween modalities"]


def bmmc_scores():
    # Single-file mode (default): load the shipped final score frame and plot directly.
    # Set REBUILD_MATRIX=1 to rebuild it from the per-metric tables (+ the published_reference splice).
    scores_csv = os.path.join(D, "fig6_scores.csv")
    if os.environ.get("REBUILD_MATRIX", "0") != "1" and os.path.exists(scores_csv):
        tab = pd.read_csv(scores_csv)
    else:
        sm = pd.read_csv(os.path.join(D, "sum_metrics.csv")); sm.columns = [c.strip('"') for c in sm.columns]
        sm["method"] = sm["method"].astype(str).str.strip('"')
        ct = pd.read_csv(os.path.join(D, "celltype_metrics.csv"))
        ctm = ct.groupby("method").agg(ks_celltype_mean=("ks_celltype", "mean"),
                                       ks_sample_mean=("ks_sample", "mean"),
                                       ks_inter_celltype_mean=("ks_inter_celltype", "mean")).reset_index()
        ac = pd.read_csv(os.path.join(D, "adj_atac_predaccu.csv")); ac = ac.rename(columns={ac.columns[0]: "method"})
        ac = ac[["method", "average_accu"]]
        pk = pd.read_csv(os.path.join(D, "peakdist.csv")); pk.columns = [c.strip('"') for c in pk.columns]
        pk["method"] = pk["method"].astype(str).str.strip('"'); pk = pk[["method", "peakdist_adj"]]
        raw = sm.merge(ctm, on="method").merge(ac, on="method").merge(pk, on="method")
        raw[FACETS[0]] = (raw.ks_celltype_mean + raw.asw) / 2
        raw[FACETS[1]] = (raw.ks_sample_mean + raw.sample_asw) / 2
        raw[FACETS[2]] = (raw["ks.statistic"] + raw.omics_asw) / 2
        raw[FACETS[3]] = ((raw.ari + raw.ami) / 2 + raw.ks_inter_celltype_mean
                          + (raw.average_accu + raw.peakdist_adj) / 2) / 3
        raw = raw.set_index("method")[FACETS]
        pub = pd.read_csv(PUB).rename(columns={"Bio-conservation": FACETS[0], "Sample batch correction": FACETS[1],
                                               "Omics gap reduction": FACETS[2], "Alignment accuracy": FACETS[3]})
        pub = pub.set_index("method")[FACETS]

        def sc(m):
            if m in raw.index: return raw.loc[m].to_dict()
            if m in pub.index: return pub.loc[m].to_dict()
            return None
        rows = []
        for base in SELECTED:
            for setting, m in [("Default", base)] + ([("Batch-aware", BATCH[base])] if base in BATCH else []):
                s = sc(m)
                if s is None:
                    print("  skip missing:", m); continue
                rows.append({"pipeline": base, "setting": setting, "overall": np.mean([s[f] for f in FACETS]), **s})
        tab = pd.DataFrame(rows)
    order = (tab[tab.setting == "Default"].set_index("pipeline")["overall"]
             .sort_values(ascending=False).index.tolist())
    return tab, order


def draw_bmmc(axes, tab, order):
    xpos = {p: i for i, p in enumerate(order)}
    BW = 0.38
    for k, (ax, facet) in enumerate(zip(axes, FACETS)):
        for _, r in tab.iterrows():
            x = xpos[r.pipeline]
            has_batch = r.pipeline in BATCH and (tab.pipeline == r.pipeline).sum() > 1
            off = (-BW/2 if r.setting == "Default" else BW/2) if has_batch else 0.0
            ax.bar(x + off, r[facet], width=BW if has_batch else BW*1.4,
                   color=COLORS.get(r.pipeline, "#cccccc"), edgecolor="black", linewidth=0.7,
                   hatch=None if r.setting == "Default" else "////")
        ax.set_ylim(0, 0.9); ax.set_yticks([0, 0.25, 0.5, 0.75]); ax.tick_params(labelsize=FS_TICK)
        ax.text(1.008, 0.5, facet, transform=ax.transAxes, rotation=270, va="center", ha="left", fontsize=FS_FACET)
        ax.margins(x=0.01)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order if k == len(FACETS) - 1 else [], rotation=45, ha="right", fontsize=FS_XMETH)
    # Setting legend: single horizontal row ABOVE panel A, no title
    axes[0].legend(handles=[Patch(facecolor="white", edgecolor="black", label="Default"),
                            Patch(facecolor="white", edgecolor="black", hatch="////", label="Batch-aware")],
                   loc="lower center", bbox_to_anchor=(0.5, 1.04), ncol=2, frameon=False, fontsize=FS_LEG,
                   handlelength=1.6, columnspacing=1.6)


# ---------- diverging colormap WITHOUT the white centre ----------
_base = plt.cm.RdBu
cmap = LinearSegmentedColormap.from_list(
    "RdBu_nowhite", np.vstack([_base(np.linspace(0.0, 0.40, 128)), _base(np.linspace(0.60, 1.0, 128))]))

# ---------- shared scatter config ----------
# B and C share the SAME x/y range. Start at 0.52 (so Panel C's scButterfly at x=0.54 is not clipped)
# and extend to 0.95 to give the near-tied top methods (scglue/scglue(multiome)/scVI) label headroom.
SCAT_LIM   = (0.52, 0.95)
SCAT_TICKS = [0.6, 0.7, 0.8, 0.9]   # identical numeric ticks on x AND y, shared by B and C
EXPAND     = (2.0, 2.6)     # adjustText repulsion -- larger to separate the dense top-right cluster
BMMC_XLAB  = "Integration performance score\nBMMC s1d1 scRNA + s4d1 scATAC"
BMMC_YLAB  = "Integration performance score\nBMMC s1d1 sc-multiome"


def _load():
    tab, order = bmmc_scores()
    bmmc = bmmc_crosssite_data()          # B: (methods, x, y, rank_change) same-donor cross-site
    pbmc = pbmc_data()                    # C: (methods, x, y, rank_change) cross-platform
    dmax = max(1, int(max(bmmc[3].abs().max(), pbmc[3].abs().max())))   # ONE shared colour scale
    return tab, order, bmmc, pbmc, plt.Normalize(-dmax, dmax), dmax


def _scatter_B(ax, bmmc, norm):
    plot_pbmc(ax, *bmmc, cmap, norm, panel_label=None, label_fs=17, tick_fs=FS_TICK, axis_fs=FS_AXIS,
              xlabel=BMMC_XLAB, ylabel=BMMC_YLAB, lim=SCAT_LIM, expand=EXPAND, ticks=SCAT_TICKS)


def _scatter_C(ax, pbmc, norm):
    plot_pbmc(ax, *pbmc, cmap, norm, panel_label=None, label_fs=17, tick_fs=FS_TICK, axis_fs=FS_AXIS,
              lim=SCAT_LIM, expand=EXPAND, ticks=SCAT_TICKS)   # default PBMC axis labels


def _cbar(fig, cax, norm):
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax)
    cb.set_label("Rank change vs sc-multiome", fontsize=FS_CBAR)   # red/blue meaning -> figure legend
    cb.ax.tick_params(labelsize=FS_TICK)


def _save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{stem}.{ext}"), dpi=250, bbox_inches="tight")
    print("wrote", stem + ".png/.pdf")


# ---------- 1) combined A | B/C ----------
def render_combined(tab, order, bmmc, pbmc, norm, dmax):
    fig = plt.figure(figsize=(19, 15))
    gs = fig.add_gridspec(4, 2, width_ratios=[1.15, 1.05], wspace=0.33, hspace=0.30)
    bmmc_axes = [fig.add_subplot(gs[i, 0]) for i in range(4)]
    right = gs[0:4, 1].subgridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.30)
    ax_B = fig.add_subplot(right[0]); ax_C = fig.add_subplot(right[1])

    draw_bmmc(bmmc_axes, tab, order)
    _scatter_B(ax_B, bmmc, norm); _scatter_C(ax_C, pbmc, norm)

    pB, pC = ax_B.get_position(), ax_C.get_position()
    cb_h = (pB.y1 - pC.y0) * 0.62
    cax = fig.add_axes([pB.x1 + 0.012, (pB.y1 + pC.y0) / 2 - cb_h / 2, 0.012, cb_h])
    _cbar(fig, cax, norm)

    pA, pB, pC = bmmc_axes[0].get_position(), ax_B.get_position(), ax_C.get_position()
    top_AB, left_BC = max(pA.y1, pB.y1) + 0.012, min(pB.x0, pC.x0) - 0.052
    let = dict(fontsize=FS_LETTER, fontweight="bold", fontfamily="Arial", va="bottom")
    fig.text(pA.x0 - 0.052, top_AB, "A", ha="left", **let)
    fig.text(left_BC,       top_AB, "B", ha="left", **let)
    fig.text(left_BC,       pC.y1 + 0.012, "C", ha="left", **let)
    y_mid = (bmmc_axes[0].get_position().y1 + bmmc_axes[3].get_position().y0) / 2
    fig.text(pA.x0 - 0.040, y_mid, "Benchmark metrics value", rotation=90, va="center", ha="center", fontsize=FS_YLAB)
    _save(fig, "fig6_combined")


# ---------- 2) standalone Panel A (BMMC benchmark) ----------
def render_panelA(tab, order):
    fig = plt.figure(figsize=(10, 14))
    axes = [fig.add_subplot(4, 1, i + 1) for i in range(4)]
    fig.subplots_adjust(left=0.20, right=0.90, top=0.94, bottom=0.14, hspace=0.30)
    draw_bmmc(axes, tab, order)
    pA = axes[0].get_position()
    y_mid = (axes[0].get_position().y1 + axes[3].get_position().y0) / 2
    fig.text(pA.x0 - 0.11, y_mid, "Benchmark metrics value", rotation=90, va="center", ha="center", fontsize=FS_YLAB)
    fig.text(pA.x0 - 0.11, pA.y1 + 0.012, "A", fontsize=FS_LETTER, fontweight="bold", fontfamily="Arial", va="bottom")
    _save(fig, "fig6_panelA")


# ---------- 3) standalone Panels B & C (the two shared-scale scatters) ----------
def render_BC(bmmc, pbmc, norm):
    fig = plt.figure(figsize=(8.8, 15))
    gs = fig.add_gridspec(2, 1, hspace=0.24)
    fig.subplots_adjust(left=0.16, right=0.82, top=0.96, bottom=0.06)
    ax_B = fig.add_subplot(gs[0]); ax_C = fig.add_subplot(gs[1])
    _scatter_B(ax_B, bmmc, norm); _scatter_C(ax_C, pbmc, norm)
    fig.canvas.draw()                                       # finalise aspect='equal' boxes before reading positions
    pB, pC = ax_B.get_position(), ax_C.get_position()
    cb_h = (pB.y1 - pC.y0) * 0.5
    cax = fig.add_axes([pB.x1 + 0.03, (pB.y1 + pC.y0) / 2 - cb_h / 2, 0.022, cb_h])
    _cbar(fig, cax, norm)
    let = dict(fontsize=FS_LETTER, fontweight="bold", fontfamily="Arial", va="bottom")
    fig.text(pB.x0 - 0.10, pB.y1 + 0.006, "B", **let)
    fig.text(pC.x0 - 0.10, pC.y1 + 0.006, "C", **let)
    _save(fig, "fig6_BC")


def main():
    tab, order, bmmc, pbmc, norm, dmax = _load()
    render_combined(tab, order, bmmc, pbmc, norm, dmax)
    render_panelA(tab, order)
    render_BC(bmmc, pbmc, norm)
    print("BMMC order:", order)
    print(f"shared rank-change scale: +/-{dmax}")


if __name__ == "__main__":
    main()
