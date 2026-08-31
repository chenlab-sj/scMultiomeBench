#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
Fig4D (macrophage subsampling) -- UPDATED to add MaxFuse, MIDAS, scButterfly.

Faithful port of BRCA/benchmark/old/HT243B1-S1H4/macophage_sub_fig4.ipynb, with two changes:
  1. the 3 new methods' Macrophages accuracy (new_methods_macro_accu.csv, from
     compute_new_methods_accu.py) is concatenated onto the old-method accuracy;
  2. nothing about the old methods is recomputed -- their per-(rep,sub) accuracy CSVs are
     read exactly as the old notebook read them:
       rep1  sub0..5 : BRCA/benchmark/old/HT243B1-S1H4/knn_k10sub{0..5}_pred_accu.csv
       rep2-5 sub1..5: BRCA/HT243-S1H4_subsample/HT243_S1H4_macrosub{rep}/knn_k10sub{1..5}_pred_accu.csv

Left panel  = Macrophage ATAC-label prediction accuracy vs subset fraction (mean +/- SD over
              the 5 replicate subsampling series).
Right panel = minimum macrophage count needed to still reach accuracy > 0.9 (852 * smallest
              fraction that stays > 0.9, per replicate; lower = more robust).
scDART and Cobolt stay excluded (matching the main analysis).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

G = "/path/to/multiomeBench"
OLD_REP1  = f"{G}/BRCA/benchmark/old/HT243B1-S1H4"          # rep1: knn_k10sub{0..5}_pred_accu.csv
SUBSAMPLE = f"{G}/BRCA/HT243-S1H4_subsample"                # rep2-5: HT243_S1H4_macrosub{rep}/...
HERE      = f"{G}/BRCA/benchmark/fig4d_macrophage_subset"
MACRO_BASE_COUNT = 852                                      # full test-ATAC macrophage count

# lineplot method set: old 9 (scDART + Cobolt excluded, as in the old notebook) + 3 new
PIPELINES_SEL = ["scglue(multiome)", "scVI", "scglue", "scJoint", "Seurat(CCA)",
                 "Portal", "simba", "scBridge", "BindSC", "MaxFuse", "MIDAS", "scButterfly"]

METHOD_COLORS = {
    "scglue": "#1f77b4", "scglue(multiome)": "#aec7e8", "BindSC": "#e377c2",
    "Portal": "#ffbb78", "scBridge": "#2ca02c", "simba": "#98df8a",
    "Seurat(CCA)": "#d62728", "scJoint": "#ff9896", "scVI": "#9467bd",
    "MaxFuse": "#17becf", "MIDAS": "#bcbd22", "scButterfly": "#8c564b",
}

# ---------- load OLD-method Macrophage accuracy (read the existing CSVs; do not recompute) ----------
def load_macro(path, sub, rep):
    df = pd.read_csv(path, index_col=0)
    m = df[df["type"] == "Macrophages"][["pipeline", "accuracy", "type"]].copy()
    m["sub"] = sub; m["rep"] = rep; m["percentage"] = 1 / 2 ** sub
    return m

old_macro = []
for sub in range(0, 6):                                     # rep1 has the full sub0 point
    old_macro.append(load_macro(f"{OLD_REP1}/knn_k10sub{sub}_pred_accu.csv", sub, 1))
for rep in range(2, 6):
    for sub in range(1, 6):
        old_macro.append(load_macro(f"{SUBSAMPLE}/HT243_S1H4_macrosub{rep}/knn_k10sub{sub}_pred_accu.csv", sub, rep))
old_macro = pd.concat(old_macro, ignore_index=True)

# ---------- new-method Macrophage accuracy (computed by compute_new_methods_accu.py) ----------
new_macro = pd.read_csv(f"{HERE}/new_methods_macro_accu.csv")

macro_accu = pd.concat([old_macro, new_macro], ignore_index=True)
macro_accu = macro_accu[macro_accu["pipeline"].isin(PIPELINES_SEL)].copy()
print("methods in figure:", sorted(macro_accu["pipeline"].unique()))
print("rows:", macro_accu.shape[0])

# ---------- right panel: minimum macrophage count needed for accuracy > 0.9 ----------
high = macro_accu[macro_accu["accuracy"] > 0.9].copy()
high = high.loc[high.groupby(["pipeline", "rep"])["sub"].idxmax()]     # smallest fraction still > 0.9
high["count"] = (1 / 2 ** high["sub"]) * MACRO_BASE_COUNT
high["value"] = high["count"].apply(lambda x: 500 if x > 800 else x)   # cap the >0.9-only-at-full case
mean_count = high.groupby("pipeline")["count"].mean().reset_index()
order = mean_count.sort_values("count")["pipeline"].tolist()
high["pipeline"] = pd.Categorical(high["pipeline"], categories=order, ordered=True)
never_high = sorted(set(PIPELINES_SEL) - set(high["pipeline"].dropna().unique()))
if never_high:
    print("NOTE: never reach >0.9 (absent from right panel):", never_high)

# ---------- plot ----------
plt.rc("font", size=16)
fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={"width_ratios": [2.5, 2]}, figsize=(15, 6))

sns.lineplot(data=macro_accu, x="sub", y="accuracy", hue="pipeline", marker="o",
             palette=METHOD_COLORS, ax=ax1, ci="sd", err_style="bars")   # ci="sd" for seaborn 0.11
ax1.set_xticks([0, 1, 2, 3, 4, 5])
ax1.set_xticklabels(["1", "1/2", "1/4", "1/8", "1/16", "1/32"])   # 1 on the left -> 1/32 on the right
ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1]); ax1.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1])
ax1.set_xlabel("Macrophage subset fraction")
ax1.set_ylabel("ATAC label prediction \n accuracy for macrophage")
ax1.axhline(y=0.9, color="tab:red", linestyle="--", linewidth=1.5)
ax1.legend([], [], frameon=False)

sns.barplot(data=high, y="pipeline", x="value", hue="pipeline", dodge=False,
            ci="sd", capsize=0.1, palette=METHOD_COLORS, ax=ax2)   # ci="sd" for seaborn 0.11
patches = ax2.patches
for i, line in enumerate(ax2.get_lines()):
    if i < len(patches):                                   # recolor each error bar to match its bar
        line.set_color(patches[i].get_facecolor())
ax2.set_ylabel("")
ax2.set_xticks([0, 100, 200, 300, 400, 500])
ax2.set_xticklabels(["0", "100", "200", "300", "400", ">500"])
ax2.set_xlabel("Minimum macrophage counts required for \n ATAC label prediction accuracy > 0.9")
ax2.legend([], [], frameon=False)

handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=6)
plt.tight_layout(rect=(0.02, 0.01, 0.98, 0.82))
fig.savefig(f"{HERE}/fig4d_macrophage_subset.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{HERE}/fig4d_macrophage_subset.pdf", bbox_inches="tight")
macro_accu.to_csv(f"{HERE}/fig4d_macro_accu_combined.csv", index=False)

# summary table (mean over reps at each fraction) for the record
summ = (macro_accu.groupby(["pipeline", "sub"])["accuracy"].mean().unstack("sub")
        .reindex(order[::-1] if order else None))
print("\nmean Macrophage accuracy by sub (0=full ... 5=1/32):")
print(summ.round(3).to_string())
print("\nmin macrophage count for >0.9 (mean over reps, lower=better):")
print(mean_count.sort_values("count").to_string(index=False))
print(f"\nwrote fig4d_macrophage_subset.png/.pdf and fig4d_macro_accu_combined.csv to {HERE}")
