#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""Fig2A -- pbmc3k UMAP grid.

Clean standalone extraction of the UMAP section (cells 1-14) of benchmark/benchmarkmatrix_sel.ipynb.
Grid = 5 methods (rows) x 5 panels (cols): modality masks multiomics / scRNA / scATAC, true Cell type,
Predicted ATAC cell type. The 5 rows are ordered top->bottom by Fig2b rank (BindSC dropped to fit the
figure): scglue(multiome), scVI, Seurat(CCA), scJoint, Conos.

INPUTS -- per-cell latent embeddings live on the cluster, so run this on a node with them mounted.
Edit the LATENTS dict / LABEL / KNN below to repoint (e.g. to repo latents under scripts/<method>/).
  * LATENTS[method] : coembed latent csv, index = <barcode>_rna / <barcode>_atac, cols = latent dims
  * LABEL           : label.csv (Barcode index, cell_type, modality)
  * KNN             : KNN_res_label_pred.csv (per-method predicted ATAC cell type, cols named by method)
Needs scanpy (computes neighbors + UMAP per method).

Usage:  python make_umap_grid.py [outdir=.]
"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ---- config: inputs (default = the original cluster paths from the notebook) --------------------
BASE = "/path/to/tools"
LATENTS = {                                                       # 5 methods, rows ordered by Fig2b rank
    "scglue(multiome)": f"{BASE}/scglue/scglue_withpair/pbmc3k/scglue_latent.csv",  # Fig2b rank 1
    "scVI":             f"{BASE}/scvi/pbmc3k/res_pbmc3k/scvi_latent.csv",           # Fig2b rank 3
    "Seurat(CCA)":      f"{BASE}/Seuratv3/pbmc3k/res_pbmc3k_testall/lat_df.csv",    # Fig2b rank 5
    "scJoint":          f"{BASE}/scJoint/pbmc3k/output/scJoint_latent.csv",         # Fig2b rank 8
    "Conos":            f"{BASE}/Concos/pbmc3k/res_pbmc3k/coembed_coor.csv",        # Fig2b rank 18
}                                                                  # BindSC removed to fit the figure
LABEL    = f"{BASE}/benchmark/pbmc3k/label.csv"                   # also in repo: ../fig2b/label.csv
# per-method KNN predicted ATAC labels from the CURRENT fig2b outputs (knn_pred_label__<method>.csv:
# index = <bc>-1_atac, single column named by the method). The old single wide KNN_res_label_pred.csv
# used legacy names ("Seurat v3", "Concos", no scVI) that don't match the LATENTS keys -> KeyError.
KNN_DIR  = os.environ.get("KNN_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fig2b"))
OUT      = sys.argv[1] if len(sys.argv) > 1 else "."
FILENAME = "umap.png"

# ---- labels + KNN predictions ------------------------------------------------------------------
labels_annot = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_sel       = labels_annot.index.values
labels_test  = labels_annot[labels_annot["modality"] != "train multiomics"]
celltype     = labels_test["cell_type"].unique().tolist()        # category order for Cell type panels

pipelines = list(LATENTS.keys())
# one column per method (col name = method); concat aligns on the <bc>-1_atac index
KNN_atac_label = pd.concat(
    [pd.read_csv(os.path.join(KNN_DIR, f"knn_pred_label__{m}.csv"), index_col=0) for m in pipelines], axis=1)
lat_dfs   = [pd.read_csv(LATENTS[m], index_col=0) for m in pipelines]

cell_annot, data_annot = "cell_type", "modality"
omics = ["multiomics", "scRNA", "scATAC"]
panel = ["multiomics", "scRNA", "scATAC", "Cell type", "Predicted ATAC cell type"]
n_cols, n_rows = 5, len(lat_dfs)

# ---- UMAP grid (faithful to notebook cell 14) --------------------------------------------------
fig, axs = plt.subplots(n_rows, n_cols, figsize=(13, 7.3))   # taller so the vertical method names fit at font 12; legends on top
plt.rc("font", size=14)
for i in range(len(lat_dfs)):
    row = i
    lat_df = lat_dfs[i]
    lat_df = lat_df.loc[lat_df.index.isin(bc_sel)]
    pipeline = pipelines[i]
    print("work on " + pipeline)
    adata = sc.AnnData(lat_df)
    KNN_atac = KNN_atac_label[pipeline]
    KNN_atac.name = "Predicted ATAC cell type"
    adata_pp = adata.copy()
    sc.pp.neighbors(adata_pp)
    adata_pp.obs = pd.merge(adata_pp.obs, KNN_atac, how="left", left_index=True, right_index=True)
    sc.tl.umap(adata_pp)

    annot_df_reorder = labels_annot.reindex(adata_pp.obs.index)
    adata_pp.obs["Cell type"] = annot_df_reorder[cell_annot]
    adata_pp.obs["Modality"]  = annot_df_reorder[data_annot]
    adata_pp.obs["Cell type"] = pd.Categorical(adata_pp.obs["Cell type"], categories=celltype, ordered=True)
    adata_pp.obs["Predicted ATAC cell type"] = pd.Categorical(
        adata_pp.obs["Predicted ATAC cell type"], categories=celltype, ordered=True)
    adata_pp.obs["multiomics"] = adata_pp.obs["Modality"].apply(lambda x: "multiomics" if x == "train multiomics" else "")
    adata_pp.obs["scRNA"]      = adata_pp.obs["Modality"].apply(lambda x: "scRNA" if x == "test scRNA" else "")
    adata_pp.obs["scATAC"]     = adata_pp.obs["Modality"].apply(lambda x: "scATAC" if x == "test scATAC" else "")
    for col in ("multiomics", "scRNA", "scATAC"):
        adata_pp.obs[col] = pd.Categorical(adata_pp.obs[col], categories=omics, ordered=True)

    for j in range(len(panel)):
        ax = axs[row, j]
        if j <= 2:
            sc.pl.umap(adata_pp, color=panel[j], palette=sns.color_palette("Set2"), show=False, ax=ax, title=None, frameon=False)
        else:
            sc.pl.umap(adata_pp, color=panel[j], show=False, ax=ax, title=None, frameon=False)
        if j == 0:
            ax.axis("on")
            ax.tick_params(top="off", bottom="off", left="off", right="off", labelleft="on", labelbottom="off")
            ax.set_ylabel(pipeline, rotation=90, va="center", fontsize=12)   # vertical method name, bigger font
            ax.set_xlabel("")
            ax.set(frame_on=False)
        ax.set_title(None)
        if row == 0:
            ax.set_title(panel[j], fontsize=14)
            if j == 0:
                handles1, labels1 = ax.get_legend_handles_labels()
            if j == 3:
                handles2, labels2 = ax.get_legend_handles_labels()
        ax.legend_.remove()

fig.text(0.5, 0.02, "UMAP1", ha="center", fontsize=14)
fig.text(0.01, 0.5, "UMAP2", va="center", rotation="vertical", fontsize=14)

# drop the "NA" entry (empty modality mask / unlabelled cells) that the categorical legends pick up
def _drop_na(hs, ls):
    keep = [(h, l) for h, l in zip(hs, ls) if str(l) not in ("NA", "nan", "")]
    return [h for h, _ in keep], [l for _, l in keep]
handles1, labels1 = _drop_na(handles1, labels1)
handles2, labels2 = _drop_na(handles2, labels2)

# both legends on TOP, 2 rows each (modality: 3 items -> ncol 2; cell type: 7 items -> ncol 4)
fig.legend(handles1, labels1, loc="upper left", bbox_to_anchor=(0.03, 0.99), ncol=2, fontsize=11, title=panel[0])
fig.legend(handles2, labels2, loc="upper left", bbox_to_anchor=(0.27, 0.99), ncol=4, fontsize=11, title=panel[3])
fig.tight_layout(rect=(0.02, 0.03, 0.98, 0.85))   # top room for the 2-row legends
plt.subplots_adjust(hspace=0.1)
plt.savefig(os.path.join(OUT, FILENAME), dpi=300)
print("wrote", os.path.join(OUT, FILENAME))
