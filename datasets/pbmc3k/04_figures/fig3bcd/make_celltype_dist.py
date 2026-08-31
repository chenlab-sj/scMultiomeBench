#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""Fig3 B / C / D -- pbmc3k per-cell-type inter-omics distance panels.

Clean standalone extraction of the distance-plot cells (27,28,30,33,34,36) of
benchmark/benchmarkmatrix_sel.ipynb. Curated 4-method subset (unchanged, per request):
scglue(multiome), scJoint, scBridge, cobolt.

Produces four files (same names as the notebook):
  celltype_samemodality_dist0303.pdf  (cell 30) intra-omics intra/inter-celltype distance hists
  CD16_inter_dist0303.pdf             (cell 33) Fig3B  CD16 Monocytes inter-omics violins
  T_cell_dist0303.pdf                 (cell 34) Fig3C  Naive CD8 T cells inter-omics violins
  Memory_T_inter_dist.png             (cell 36) Fig3D  Memory T cells inter-omics distance hists

The notebook called benchmark_fun.lat2pdist / .get_celltype_dist (a cluster-only module). Those are
REPRODUCED inline below so this script is self-contained; to use the exact cluster module instead,
set USE_BENCHMARK_FUN=True (needs /path/to/multiomeBench/common on sys.path).

INPUTS -- per-cell latents live on the cluster; edit LATENTS / LABEL to repoint. Run on a node.
Usage:  python make_celltype_dist.py [outdir=.]
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances

# ---- config -------------------------------------------------------------------------------------
BASE = "/path/to/tools"
LATENTS = {                                                       # curated 4-method subset (order = columns)
    "scglue\n(multiome)": f"{BASE}/scglue/scglue_withpair/pbmc3k/scglue_latent.csv",
    "scJoint":            f"{BASE}/scJoint/pbmc3k/output/scJoint_latent.csv",
    "scBridge":           f"{BASE}/scBridge/pbmc3k/latent.csv",
    "cobolt":             f"{BASE}/cobolt/pbmc3k/res_pbmc3k/cobolt_latent.csv",
}
LABEL = f"{BASE}/benchmark/pbmc3k/label.csv"                      # also in repo: ../fig2b/label.csv
OUT   = sys.argv[1] if len(sys.argv) > 1 else "."
option = "cosine"
rna_pipelines = []                                               # (notebook: empty -> use both modalities)
USE_BENCHMARK_FUN = False

# ---- benchmark_fun reproductions (self-contained) ----------------------------------------------
if USE_BENCHMARK_FUN:
    sys.path.append(f"{BASE}/benchmark")
    import benchmark_fun
    lat2pdist = benchmark_fun.lat2pdist
    get_celltype_dist = benchmark_fun.get_celltype_dist
else:
    def lat2pdist(lat_df, bc, option):
        """pairwise distance matrix over the barcodes in bc (index == columns == bc order)."""
        df = lat_df.loc[lat_df.index.isin(bc)]
        D = cosine_distances(df.values) if option.lower() == "cosine" else euclidean_distances(df.values)
        return pd.DataFrame(D, index=df.index, columns=df.index)

    def _col_typedist(col, label_dic):
        wi, bt = [], []
        col_type = label_dic.get(col.name)
        for index, value in col.items():
            if label_dic.get(index) == col_type:
                wi.append(value)
            else:
                bt.append(value)
        return wi, bt

    def _combine(row):
        out = []
        for sub in row:
            out.extend(sub)
        return out

    def get_celltype_dist(pdist_df, label_dic):
        """-> [within-celltype distances, between-celltype distances] pooled over all columns."""
        dist_out = pdist_df.apply(lambda column: _col_typedist(column, label_dic))
        results = dist_out.apply(lambda row: _combine(row), axis=1)
        return [results.loc[0], results.loc[1]]

# ---- labels -------------------------------------------------------------------------------------
labels_annot = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
labels_df = labels_annot
bc_test = labels_annot[labels_annot["modality"] != "train multiomics"].index.values
bc_test_rna = [bc for bc in bc_test if bc.endswith("_rna")]
bc_test_atac = [bc for bc in bc_test if bc.endswith("_atac")]
labels_test = labels_annot.loc[bc_test]
label_dic = dict(zip(labels_test.index, labels_test["cell_type"]))
cell_test = [s.split("_")[0] for s in bc_test]

pipelines = list(LATENTS.keys())
lat_dfs = [pd.read_csv(p, index_col=0) for p in LATENTS.values()]


def _inter_pdist(lat_df):
    """rna-rows x atac-cols inter-omics distance block (shared by the three inter-omics panels)."""
    pdist_df = lat2pdist(lat_df, bc_test, option)
    d = pdist_df.loc[pdist_df.index.isin(bc_test_rna)]
    return d[[c for c in d.columns if c in bc_test_atac]]


def _atac_cols(celltype):
    return [k for k in label_dic if label_dic[k] == celltype and k.endswith("_atac")]


# ================================ cell 30: same-modality hists ===================================
def fig_samemodality():
    celltype_sel = ["CD16 Monocytes", "Naive CD8 T cells"]
    celltype_sel_label = ["CD16 \n Monocytes", "Naive \n CD8 T cells"]
    fig2, ax2s = plt.subplots(2, len(lat_dfs), figsize=(12, 5)); plt.rc("font", size=15)
    for i, lat_df in enumerate(lat_dfs):
        pipeline = pipelines[i]; print("samemodality:", pipeline)
        pdist_df = lat2pdist(lat_df, bc_test, option)
        pd_rna = pdist_df.loc[pdist_df.index.isin(bc_test_rna)]; pd_rna = pd_rna[[c for c in pd_rna.columns if c in bc_test_rna]]
        pd_atac = pdist_df.loc[pdist_df.index.isin(bc_test_atac)]; pd_atac = pd_atac[[c for c in pd_atac.columns if c in bc_test_atac]]
        for j, type in enumerate(celltype_sel):
            catac = _atac_cols(type); crna = [k for k in label_dic if label_dic[k] == type and k.endswith("_rna")]
            rna_d = get_celltype_dist(pd_rna[[c for c in pd_rna.columns if c in crna]], label_dic)
            atac_d = get_celltype_dist(pd_atac[[c for c in pd_atac.columns if c in catac]], label_dic)
            if pipeline in rna_pipelines:
                wi, bt = atac_d[0], atac_d[1]
            else:
                wi, bt = atac_d[0] + rna_d[0], atac_d[1] + rna_d[1]
            ax2 = ax2s[j, i]
            if i == 0: ax2.set_ylabel(celltype_sel_label[j], fontsize=18)
            if j == 0: ax2.set_title(pipeline, fontsize=18)
            ax2.hist(np.array(wi), bins=50, density=True, alpha=0.5, color="tab:blue", label="intra-celltype same modalities")
            ax2.hist(np.array(bt), bins=50, density=True, alpha=0.5, color="tab:orange", label="inter-celltype same modalities")
    fig2.text(0.5, 0.01, "intra-omics cosine distance ", ha="center", fontsize=18)
    fig2.text(0.01, 0.5, "Density", va="center", rotation="vertical", fontsize=18)
    fig2.legend(labels=["intra-cell type distances", "inter-cell type distances"], loc="upper center",
                bbox_to_anchor=(0.5, 0.98), fontsize=16, ncol=2)
    plt.tight_layout(rect=(0.02, 0.02, 0.95, 0.90)); plt.subplots_adjust(hspace=0.35)
    fig2.savefig(os.path.join(OUT, "celltype_samemodality_dist0303.pdf"), dpi=400); plt.close(fig2)


# ================================ cell 33: Fig3B CD16 Monocytes ==================================
def fig_cd16():
    fig1, ax1s = plt.subplots(1, 4, figsize=(18.5, 3)); plt.rc("font", size=12)
    for i, lat_df in enumerate(lat_dfs):
        pipeline = pipelines[i]; print("CD16:", pipeline)
        inter = _inter_pdist(lat_df)
        type = "CD16 Monocytes"
        cell = inter[[c for c in inter.columns if c in _atac_cols(type)]]
        wi, bt = get_celltype_dist(cell, label_dic)                          # intra-celltype (blue)
        # CD14 Monocytes rna rows -> CD16 atac cols
        sel_rna = [k for k in label_dic if label_dic[k] == "CD14 Monocytes" and k.endswith("_rna")]
        bt_sel = get_celltype_dist(cell.loc[sel_rna], label_dic)[1]           # CD16-CD14 distance
        # all other cell types (not CD16, not CD14) rna rows -> CD16 atac cols
        sel2_rna = [k for k in label_dic if label_dic[k] not in (type, "CD14 Monocytes") and k.endswith("_rna")]
        bt_sel2 = get_celltype_dist(cell.loc[sel2_rna], label_dic)[1]         # CD16-other distance
        ax1 = ax1s[i]
        sns.violinplot(data=[np.array(wi), np.array(bt_sel), np.array(bt_sel2)], ax=ax1)
        ax1.set(xlabel=None); ax1.tick_params(bottom=False); ax1.set_xticklabels([]); ax1.set_title(pipeline, fontsize=16)
    fig1.text(0.5, 0.02, " Inter-omics cosine distance ", ha="center", fontsize=18)
    fig1.text(0.02, 0.5, "Density", va="center", rotation="vertical", fontsize=18)
    fig1.legend(labels=["CD16 Monocyte \n intra-cell type distance", "CD16-CD14 Monocyte distance",
                        "CD16 Monocyte - \n other cell type distance"], loc="right",
               bbox_to_anchor=(0.95, 0.4), fontsize=14, ncol=1)
    plt.tight_layout(rect=(0.04, 0.05, 0.7, 0.98))
    fig1.savefig(os.path.join(OUT, "CD16_inter_dist0303.pdf"), dpi=300); plt.close(fig1)


# ================================ cell 34: Fig3C Naive CD8 T cells ===============================
def fig_tcell():
    fig1, ax1s = plt.subplots(1, 4, figsize=(18.5, 3)); plt.rc("font", size=12)
    for i, lat_df in enumerate(lat_dfs):
        pipeline = pipelines[i]; print("Tcell:", pipeline)
        inter = _inter_pdist(lat_df)
        type = "Naive CD8 T cells"
        cell = inter[[c for c in inter.columns if c in _atac_cols(type)]]
        wi, bt = get_celltype_dist(cell, label_dic)
        sel_rna = [k for k in label_dic if label_dic[k] == "Naive CD4 T cells" and k.endswith("_rna")]
        bt_sel = get_celltype_dist(cell.loc[sel_rna], label_dic)[1]           # CD8-CD4
        sel2_rna = [k for k in label_dic if label_dic[k] == "Memory T cells" and k.endswith("_rna")]
        bt_sel2 = get_celltype_dist(cell.loc[sel2_rna], label_dic)[1]         # CD8-Memory
        sel3_rna = [k for k in label_dic if label_dic[k] not in (type, "Naive CD4 T cells", "Memory T cells") and k.endswith("_rna")]
        bt_sel3 = get_celltype_dist(cell.loc[sel3_rna], label_dic)[1]         # CD8-other
        ax1 = ax1s[i]
        sns.violinplot(data=[np.array(wi), np.array(bt_sel), np.array(bt_sel2), np.array(bt_sel3)], ax=ax1)
        ax1.set(xlabel=None); ax1.tick_params(bottom=False); ax1.set_xticklabels([]); ax1.set_title(pipeline, fontsize=16)
    fig1.text(0.5, 0.02, "Inter-omics cosine distance", ha="center", fontsize=18)
    fig1.text(0.02, 0.5, "Density", va="center", rotation="vertical", fontsize=18)
    fig1.legend(labels=["Naive CD8 T cells \n intra-cell type distance", "Naive CD8 - Naive CD4 T cells distance",
                        "Naive CD8 - Memory T cells distance", "Naive CD8 T cell - \n other cell types distance"],
               loc="right", bbox_to_anchor=(0.95, 0.4), fontsize=14, ncol=1)
    plt.tight_layout(rect=(0.04, 0.05, 0.7, 0.98))
    fig1.savefig(os.path.join(OUT, "T_cell_dist0303.pdf"), dpi=400); plt.close(fig1)


# ================================ cell 36: Fig3D Memory T cells ==================================
def fig_memoryT():
    fig1, ax1s = plt.subplots(3, 6, figsize=(16, 10)); plt.rc("font", size=12)
    for i, lat_df in enumerate(lat_dfs):
        row1, col1 = i // 6, i % 6
        pipeline = pipelines[i]; print("MemoryT:", pipeline)
        inter = _inter_pdist(lat_df)
        type = "Memory T cells"
        cell = inter[[c for c in inter.columns if c in _atac_cols(type)]]
        wi, bt = get_celltype_dist(cell, label_dic)
        sel = ["Naive CD4 T cells", "Naive CD8 T cells"]
        sel_rna = [k for k in label_dic if label_dic[k] in sel and k.endswith("_rna")]
        bt_sel = get_celltype_dist(cell.loc[sel_rna], label_dic)[1]           # Memory-Naive CD4/8
        ax1 = ax1s[row1, col1]
        ax1.hist(np.array(bt), bins=50, density=False, alpha=0.5, color="tab:orange", label="Memory T cells inter-cell type inter-omics distance")
        ax1.hist(np.array(bt_sel), bins=50, density=False, alpha=0.5, color="tab:red", label="Memory - Naive CD4/8 T cells inter-omics distance")
        ax1.hist(np.array(wi), bins=50, density=False, alpha=0.5, color="tab:blue", label="Memory T cells inter-omics distance")
        ax1.set_title(pipeline, fontsize=12)
    fig1.text(0.5, 0.02, option + " distance ", ha="center", fontsize=12)
    fig1.text(0.02, 0.5, "Density", va="center", rotation="vertical", fontsize=12)
    fig1.legend(labels=["Memory T cells inter-cell type inter-omics distance", "Memory - Naive CD4/8 T cells inter-omics distance",
                        "Memory T cells inter-omics distance"], loc="upper center", bbox_to_anchor=(0.5, 0.99), fontsize=14, ncol=2)
    plt.tight_layout(rect=(0.01, 0.02, 0.99, 0.89))
    fig1.savefig(os.path.join(OUT, "Memory_T_inter_dist.png"), dpi=300); plt.close(fig1)


if __name__ == "__main__":
    fig_samemodality()
    fig_cd16()      # Fig3B
    fig_tcell()     # Fig3C
    fig_memoryT()   # Fig3D
    print("wrote 4 files to", OUT)
