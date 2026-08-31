#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""Composite R1#1 BRCA/PBMC generalization figure with a shared rank-change legend."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr


G = Path("/path/to/multiomeBench")
OUT = G / "BRCA" / "benchmark"
EXCLUDE = {"scDART", "Cobolt"}


# The 6 metrics shared by the paired (multiome) AND unpaired matrices. The paired matrices additionally
# carry ks.statistic/ari/ami (SAME-CELL metrics), which the unpaired data can't have. To compare paired
# vs unpaired fairly (Panels B/C are unpaired-vs-multiome scatters), we RECOMPUTE every score from these
# 6 columns with the unpaired composite formula (plot_metrics_matrix.R, SPLICE_PUBLISHED=0) rather than reading
# the stored `score` -- so the multiome baseline is not credited for metrics the unpaired point lacks.
SIX = ["ks_celltype_mean", "asw", "omics_asw", "ks_inter_celltype_mean", "average_accu", "peakdist_adj"]


def load_score(path):
    d = pd.read_csv(path)
    d.columns = [c.strip('"') for c in d.columns]
    d["method"] = d["method"].astype(str).str.strip('"')
    d = d.set_index("method")
    if all(c in d.columns for c in SIX):
        f = d[SIX].astype(float)
        # unpaired 6-metric composite: bio-conservation + omics-gap (ks.statistic dropped) + alignment
        # (ari/ami dropped). Identical to plot_metrics_matrix.R's SPLICE_PUBLISHED=0 score.
        score = ((f.ks_celltype_mean + f.asw) / 2
                 + f.omics_asw
                 + (f.ks_inter_celltype_mean + (f.average_accu + f.peakdist_adj) / 2) / 2) / 3
        return score
    return d["score"].astype(float)   # fallback: matrix lacks the 6 columns (shouldn't happen)


def spread(vals, minsep):
    out = list(vals)
    for i in range(1, len(out)):
        if out[i] - out[i - 1] < minsep:
            out[i] = out[i - 1] + minsep
    # recenter the spread block on the dots' mean so labels sit ABOVE and BELOW their points
    # (not all pushed upward) -> each label ends up closer to its own dot
    shift = (sum(vals) - sum(out)) / len(out)
    return [o + shift for o in out]


def brca_data():
    cond = {
        # figS4a: plot_metrics_matrix.R now splices the 7-major figS2a source (old) + major new -> fig2b_matrix.csv is correct.
        "HT243\nsc-multiomics": load_score(OUT / "figS4a" / "fig2b_matrix.csv"),
        # figS4b + crossdonor: user's plot-script outputs are already correct (old = published figS2c / all-major new)
        "HT137 scRNA\n+\nHT137 scATAC": load_score(OUT / "figS4b" / "metrics_matrix.csv"),
        "HT137 scRNA\n+\nHT243 scATAC": load_score(OUT / "crossdonor" / "metrics_matrix.csv"),
    }
    conds = list(cond)
    scores = pd.DataFrame(cond).dropna()
    scores = scores.loc[[m for m in scores.index if m not in EXCLUDE]]
    ranks = scores.rank(ascending=False, method="first").astype(int)
    # per-condition rank change vs the sc-multiomics baseline (conds[0]); DataFrame, baseline column = 0
    rank_change = ranks.sub(ranks[conds[0]], axis=0)
    return conds, scores, rank_change


def pbmc_data():
    multiome = load_score(G / "pbmc" / "pbmc3k" / "benchmark" / "fig2b" / "fig2b_matrix.csv")
    parse = load_score(G / "pbmc_parse" / "benchmark" / "metrics" / "metrics_matrix.csv")
    methods = [m for m in multiome.index if m in parse.index and m not in EXCLUDE]
    x = parse[methods].astype(float)
    y = multiome[methods].astype(float)
    rank_multiome = multiome[methods].rank(ascending=False)
    rank_parse = parse[methods].rank(ascending=False)
    rank_change = rank_parse - rank_multiome
    return methods, x, y, rank_change


def bmmc_crosssite_data():
    """BMMC same-donor (d1) technical-batch scatter: x = s1d1 RNA + s4d1 ATAC (unpaired cross-site),
    y = s1d1 sc-multiome (paired baseline). Both scored on the SAME 6 metrics (each matrix built with
    SPLICE_PUBLISHED=0), so this is the metric-fair analogue of the PBMC cross-platform panel."""
    crosssite = load_score(G / "BMMC_d1" / "benchmark" / "crosssite" / "metrics_matrix.csv")
    paired    = load_score(G / "BMMC_d1" / "benchmark" / "s1d1_paired" / "metrics_matrix.csv")
    methods = [m for m in paired.index if m in crosssite.index and m not in EXCLUDE]
    x = crosssite[methods].astype(float)                 # unpaired cross-site
    y = paired[methods].astype(float)                    # paired multiome baseline
    rank_change = crosssite[methods].rank(ascending=False) - paired[methods].rank(ascending=False)
    return methods, x, y, rank_change


def brca_indep_data():
    """BRCA supplement scatter (same design as Fig 6 B/C): x = HT137 scRNA + HT137 scATAC
    (independently sequenced), y = HT243 sc-multiome baseline. Both scored on the 6 shared metrics
    (load_score recomputes the unpaired composite), so paired vs independent is compared fairly."""
    indep    = load_score(OUT / "figS4b" / "metrics_matrix.csv")     # HT137 scRNA + HT137 scATAC
    multiome = load_score(OUT / "figS4a" / "fig2b_matrix.csv")       # HT243 sc-multiome baseline
    methods = [m for m in multiome.index if m in indep.index and m not in EXCLUDE]
    x = indep[methods].astype(float)                                 # independently sequenced
    y = multiome[methods].astype(float)                              # paired sc-multiome baseline
    rank_change = indep[methods].rank(ascending=False) - multiome[methods].rank(ascending=False)
    return methods, x, y, rank_change


def plot_brca(ax, conds, scores, rank_change, cmap, norm, panel_label="A",
              label_fs=7.2, tick_fs=8.5, axis_fs=10, letter_fs=13, minsep_frac=0.045,
              fit_labels=False):
    nC = len(conds)
    for method in scores.index:
        ys = [scores.loc[method, c] for c in conds]
        # each segment is coloured by the rank change (vs the baseline) of the condition it arrives at,
        # so HT137+HT137 and HT137+HT243 each show their own rank-change colour
        for i in range(nC - 1):
            ax.plot([i, i + 1], [ys[i], ys[i + 1]], "-",
                    color=cmap(norm(rank_change.loc[method, conds[i + 1]])), lw=2, alpha=0.9, zorder=3)
        for i, c in enumerate(conds):
            ax.plot(i, ys[i], "o", ms=6, zorder=3,
                    markerfacecolor=cmap(norm(rank_change.loc[method, c])),
                    markeredgecolor="0.35", markeredgewidth=0.4)

    span = scores.max().max() - scores.min().min()
    sep = span * minsep_frac
    spread_y = []                                            # track spread label positions
    for side, col, ha, dx in [
        (0, conds[0], "right", -0.07),
        (nC - 1, conds[-1], "left", 0.07),
    ]:
        order = scores[col].sort_values().index
        yvals = spread(scores.loc[order, col].values, sep)
        spread_y += list(yvals)
        for method, yy in zip(order, yvals):
            y_dot = scores.loc[method, col]
            ax.plot([side + dx, side], [yy, y_dot], color="0.7", lw=0.5, zorder=2)   # leader: label -> its dot
            ax.text(side + dx, yy, method, ha=ha, va="center", fontsize=label_fs, zorder=4)

    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels(conds, fontsize=tick_fs)
    ax.tick_params(axis="y", labelsize=tick_fs)
    ax.set_ylabel("Integration performance score", fontsize=axis_fs)
    ax.set_xlim(-1.15, len(conds) - 1 + 1.15)
    if fit_labels:                                           # expand y so every spread label AND data point is inside
        lo = min(min(spread_y), scores.min().min())          # include the middle-column low point (else it clips)
        hi = max(max(spread_y), scores.max().max())
        pad = (hi - lo) * 0.05
        ax.set_ylim(lo - pad, hi + pad)
    else:
        ax.set_ylim(scores.min().min() - span * 0.08, scores.max().max() + span * 0.12)
    if panel_label:
        ax.text(-0.12, 1.02, panel_label, transform=ax.transAxes, fontsize=letter_fs, fontweight="bold")


def plot_pbmc(ax, methods, x, y, rank_change, cmap, norm, panel_label="B",
              label_fs=7.2, tick_fs=8.5, axis_fs=9.2, letter_fs=13,
              xlabel="Integration performance score\nParse PBMC scRNA + 10X PBMC scATAC",
              ylabel="Integration performance score\n10X PBMC sc-multiomics",
              lim=(0.5, 0.85), expand=(1.4, 1.7), ticks=None):
    lo, hi = lim
    ax.plot([lo, hi], [lo, hi], ls="--", color="#888888", lw=1.0, zorder=1)
    ax.scatter(
        x,
        y,
        s=82,
        c=[cmap(norm(rank_change[m])) for m in methods],
        edgecolor="black",
        lw=0.5,
        zorder=3,
    )

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    if ticks is not None:                 # identical numeric ticks on BOTH x and y
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
    ax.tick_params(labelsize=tick_fs)
    ax.set_xlabel(xlabel, fontsize=axis_fs)
    ax.set_ylabel(ylabel, fontsize=axis_fs)

    # method labels: auto-repelled with adjustText so none overlap, thin leader line back to each point.
    # (replaces the old hand-tuned label_pos offsets, which collided once the font was enlarged.)
    from adjustText import adjust_text
    xs = [float(x[m]) for m in methods]
    ys = [float(y[m]) for m in methods]
    texts = [ax.text(xs[i], ys[i], m, fontsize=label_fs, zorder=5) for i, m in enumerate(methods)]
    adjust_text(
        texts, x=xs, y=ys, ax=ax,
        expand=expand,
        arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5),
    )

    if panel_label:
        ax.text(-0.14, 1.02, panel_label, transform=ax.transAxes, fontsize=letter_fs, fontweight="bold")


def main():
    brca_conds, brca_scores, brca_change = brca_data()
    pbmc_methods, pbmc_x, pbmc_y, pbmc_change = pbmc_data()
    dmax = max(1, int(max(brca_change.abs().values.max(), pbmc_change.abs().max())))
    cmap = plt.cm.RdBu
    norm = plt.Normalize(-dmax, dmax)

    fig = plt.figure(figsize=(12.2, 5.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.24, 1.0, 0.045], wspace=0.34)
    ax_brca = fig.add_subplot(gs[0, 0])
    ax_pbmc = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    plot_brca(ax_brca, brca_conds, brca_scores, brca_change, cmap, norm)
    plot_pbmc(ax_pbmc, pbmc_methods, pbmc_x, pbmc_y, pbmc_change, cmap, norm)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(
        "Rank change vs sc-multiomics\nred: rank decrease (better)\nblue: rank increase (worse)",
        fontsize=9,
    )

    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"r1_brca_pbmc_composite.{ext}", dpi=300, bbox_inches="tight")
    print(f"wrote {OUT / 'r1_brca_pbmc_composite.png'}")
    print(f"wrote {OUT / 'r1_brca_pbmc_composite.pdf'}")
    print(f"BRCA rho consecutive: {[round(spearmanr(brca_scores[brca_conds[i]], brca_scores[brca_conds[i+1]]).correlation, 2) for i in range(len(brca_conds)-1)]}")
    print(f"PBMC rho: {spearmanr(pbmc_x, pbmc_y).correlation:.2f}")


if __name__ == "__main__":
    main()
