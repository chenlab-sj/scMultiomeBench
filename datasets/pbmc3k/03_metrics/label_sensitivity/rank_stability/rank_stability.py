#!/usr/bin/env python3
"""
R4 label-circularity rank stability: does the method ranking change when the cell-type labels are derived by
an INDEPENDENT RNA-only method (SingleR) instead of the 10X-provided (multimodal, RNA+ATAC) labels?

Assembles the SAME composite (plot_metrics_matrix.R formula) twice -- once from the 10X-labelled fig2b metrics, once
from the SingleR-labelled recompute (this dir) -- and Spearmans the two rankings. Peak similarity is held at
the 10X value for BOTH (peak was not re-run), so the comparison isolates the label effect on the cell-type
metrics. High rho => the ranking is NOT an artifact of how the labels were derived, and (since SingleR never
sees ATAC) not driven by ATAC-shaped labels either.

NB "10X" here = the 10X Genomics R&D annotation used as ground truth (NOT Seurat WNN, which is a benchmarked
*method*). Run AFTER recompute_sub.sh.  python3.9 rank_stability.py
Outputs: rank_stability.csv + rank_stability.png/pdf
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE  = os.path.dirname(os.path.abspath(__file__))
FIG2B = os.path.normpath(os.path.join(HERE, "..", "..", "fig2b"))     # 10X-labelled metrics + peakdist


def _clean(df):
    df.columns = [c.strip('"') for c in df.columns]
    if "method" in df.columns:
        df["method"] = df["method"].astype(str).str.strip('"')
    return df


def assemble(sum_f, ct_f, accu_f, peak_f):
    """Reproduce the plot_metrics_matrix.R composite from the 4 source CSVs -> DataFrame[method, score]."""
    sm = _clean(pd.read_csv(sum_f))[["method", "asw", "omics_asw", "ami", "ari", "ks.statistic"]]
    ct = _clean(pd.read_csv(ct_f)).groupby("method", as_index=False).agg(
        ks_celltype_mean=("ks_celltype", "mean"),
        ks_inter_celltype_mean=("ks_inter_celltype", "mean"))
    ac = _clean(pd.read_csv(accu_f))[["method", "average_accu"]]
    pk = _clean(pd.read_csv(peak_f))
    pk = pk[pk["method"] != "random"][["method", "peakdist_adj"]]
    d = sm.merge(ct, on="method").merge(ac, on="method").merge(pk, on="method")
    # exact plot_metrics_matrix.R composite (paired pbmc3k: includes ks.statistic + (ari+ami)/2)
    d["score"] = ((d["ks_celltype_mean"] + d["asw"]) / 2                                   # bio-conservation
                  + (d["ks.statistic"] + d["omics_asw"]) / 2                               # omics gap reduction
                  + ((d["ari"] + d["ami"]) / 2 + d["ks_inter_celltype_mean"]
                     + (d["average_accu"] + d["peakdist_adj"]) / 2) / 3) / 3               # alignment accuracy
    return d[["method", "score"]]


def main():
    tenx = assemble(f"{FIG2B}/sum_metrics.csv", f"{FIG2B}/celltype_metrics.csv",
                    f"{FIG2B}/adj_atac_predaccu.csv", f"{FIG2B}/peakdist.csv").rename(columns={"score": "label_10X"})
    sgl  = assemble(f"{HERE}/sum_metrics.csv", f"{HERE}/celltype_metrics.csv",
                    f"{HERE}/adj_atac_predaccu.csv", f"{FIG2B}/peakdist.csv").rename(columns={"score": "label_SingleR"})

    cmp = tenx.merge(sgl, on="method").set_index("method").dropna()
    cmp["rank_10X"]     = cmp["label_10X"].rank(ascending=False).astype(int)
    cmp["rank_SingleR"] = cmp["label_SingleR"].rank(ascending=False).astype(int)
    cmp["rank_shift"]   = cmp["rank_SingleR"] - cmp["rank_10X"]
    cmp = cmp.round(3).sort_values("label_10X", ascending=False)
    cmp.to_csv(os.path.join(HERE, "rank_stability.csv"))

    from scipy.stats import spearmanr
    rho, p = spearmanr(cmp["label_10X"], cmp["label_SingleR"])
    print(f"methods compared: {len(cmp)}")
    print(f"Spearman rank (10X vs SingleR labels): rho={rho:.3f}  (p={p:.2g})")
    print(f"max |rank shift|: {cmp['rank_shift'].abs().max()}   mean |rank shift|: {cmp['rank_shift'].abs().mean():.2f}")
    print(cmp.to_string())

    n = len(cmp)
    shift_max = int(cmp["rank_shift"].abs().max())

    # ---- PRIMARY (panel D): score scatter with repelled labels. On y=x -> ranking independent of the
    #      label-derivation method. adjustText repels the 23 method labels and draws leader lines. ----
    from adjustText import adjust_text
    fig, ax = plt.subplots(figsize=(7.9, 7.7))
    allv = pd.concat([cmp["label_10X"], cmp["label_SingleR"]]); lo, hi = allv.min(), allv.max()
    pad = (hi - lo) * 0.07 + 1e-6; lims = (lo - pad, hi + pad)
    ax.plot(lims, lims, "--", color="grey", lw=1, zorder=1)
    ax.scatter(cmp["label_10X"], cmp["label_SingleR"], s=48, color="#4C72B0", zorder=3)
    ax.set_xlim(lims); ax.set_ylim(lims); ax.set_aspect("equal", "box")
    texts = [ax.text(r["label_10X"], r["label_SingleR"], m, fontsize=11.5, color="#222")
             for m, r in cmp.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.5),
                expand=(1.3, 1.7), force_text=(0.5, 0.8))
    ax.set_xlabel("Integration performance score — 10X-annotated (multimodal) labels", fontsize=14)
    ax.set_ylabel("Integration performance score — SingleR (RNA-only) labels", fontsize=14)
    ax.tick_params(labelsize=13)
    ax.set_title("PBMC 3k label source rank stability", fontsize=16)
    ax.text(0.03, 0.97, f"Spearman $\\rho$ = {rho:.3f}", transform=ax.transAxes,
            fontsize=13, ha="left", va="top")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"rank_stability.{ext}"), dpi=300)

    # ---- SECONDARY (reference): bump plot of the ranks (red = rank changed, grey = held) ----
    fig2, ax2 = plt.subplots(figsize=(6.4, 8.6))
    for m, r in cmp.iterrows():
        shifted = r["rank_shift"] != 0
        col = "#d1495b" if shifted else "#c9ccd1"
        ax2.plot([0, 1], [r["rank_10X"], r["rank_SingleR"]], "-",
                 color=col, lw=1.8 if shifted else 1.0, zorder=3 if shifted else 1, solid_capstyle="round")
        ax2.scatter([0, 1], [r["rank_10X"], r["rank_SingleR"]], s=22, color=col, zorder=4 if shifted else 2)
        ax2.text(-0.03, r["rank_10X"], m, ha="right", va="center", fontsize=7.5, color="#111" if shifted else "#444")
        ax2.text( 1.03, r["rank_SingleR"], m, ha="left", va="center", fontsize=7.5, color="#111" if shifted else "#444")
    ax2.set_xlim(-0.62, 1.62); ax2.set_ylim(n + 0.6, 0.4)
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(["10X-annotated\n(multimodal) labels", "SingleR\n(RNA-only) labels"], fontsize=10)
    ax2.set_yticks([]); ax2.tick_params(length=0)
    for s in ("top", "right", "left", "bottom"): ax2.spines[s].set_visible(False)
    ax2.set_title(f"PBMC 3k label source rank stability (bump)\nSpearman rho = {rho:.3f}, max shift {shift_max}", fontsize=10)
    fig2.tight_layout()
    for ext in ("png", "pdf"):
        fig2.savefig(os.path.join(HERE, f"rank_stability_bump.{ext}"), dpi=300)
    print("\nwrote rank_stability.csv + rank_stability.png/pdf (scatter) + rank_stability_bump.png/pdf")


if __name__ == "__main__":
    main()
