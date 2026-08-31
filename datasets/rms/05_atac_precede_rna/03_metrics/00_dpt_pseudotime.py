#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
W7 step 1 (scanpy / benchmark_env):  order RMS Mast607A cells along myogenic pseudotime (DPT) and pull
MYOD1/MYOG RNA. Produces pseudotime_table.csv (per-cell pseudotime + bin + cell_type + MYOD1/MYOG RNA)
which the R step (01_region_atac_by_bin.R) then joins with per-bin ATAC accessibility.

Data: rna_embed.h5ad = 8291 x 50 PCs, obs_names '<bc>-1_rna' (RNA-only PCA embedding of the paired multiome).
Labels: RMS/benchmark/figS4a/label.csv, 'test scRNA' rows carry Mesoderm/Myoblast/Myocyte, same '<bc>-1_rna' index.
RNA counts: the multiome filtered_feature_bc_matrix.h5 (GEX), barcodes '<bc>-1'.
  run:  python 00_dpt_pseudotime.py
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

HERE  = "/path/to/multiomeBench/RMS/benchmark/atac_precede_rna"
EMBED = "/path/to/data/RMS/Mast607A_TB19_22652/rna_embed.h5ad"
LABEL = "/path/to/multiomeBench/RMS/benchmark/figS4a/label.csv"
H5    = "/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5"
NBIN  = 15
ORDER = ["Mesoderm", "Myoblast", "Myocyte"]   # expected differentiation direction
# myogenic panel: MYOD1/MYOG (MRFs) -> MEF2C (myocyte TF) -> compact terminal structural genes
# (MYL1/TNNT3/TNNT2/ACTN2) picked for good promoter-ATAC signal. Dropped MYH3/TTN (huge genes -> weak ATAC).
GENES = ["MYOD1", "MYOG", "MYH3", "MEF2C"]

# ---- embedding + labels ----
a = ad.read_h5ad(EMBED)                                   # (8291, 50) PCs, obs '<bc>-1_rna'
lab = pd.read_csv(LABEL, index_col=0)
lab = lab[lab["modality"] == "test scRNA"]                # '<bc>-1_rna' -> Mesoderm/Myoblast/Myocyte
common = [b for b in a.obs_names if b in set(lab.index)]
a = a[common].copy()
a.obs["cell_type"] = pd.Categorical(lab.loc[common, "cell_type"].values, categories=ORDER)
print(f"cells with embedding + label: {a.n_obs}")
print(a.obs["cell_type"].value_counts().to_string())

# ---- MYOD1 / MYOG RNA (match on base barcode '<bc>-1') ----
gex = sc.read_10x_h5(H5)                                  # multiome h5 -> GEX by default
gex.var_names_make_unique()
sc.pp.normalize_total(gex, target_sum=1e4)               # library-size normalize (per 10k UMIs) to match the
                                                         # ATAC depth normalization; kept LINEAR (no log) so RNA
                                                         # and ATAC are on comparable per-10k scales before binning.
base = [b[:-4] if b.endswith("_rna") else b for b in a.obs_names]   # '<bc>-1_rna' -> '<bc>-1'
for g in GENES:
    if g in gex.var_names:
        col = np.asarray(gex[:, g].X.todense()).ravel()
        s = pd.Series(col, index=gex.obs_names)
        a.obs[f"{g}_rna"] = s.reindex(base).to_numpy()
        print(f"{g} RNA: matched {a.obs[f'{g}_rna'].notna().sum()}/{a.n_obs} cells, "
              f"mean={np.nanmean(a.obs[f'{g}_rna']):.3f}")
    else:
        print(f"WARNING: {g} not in GEX var_names")

# ---- diffusion pseudotime, rooted in Mesoderm ----
a.obsm["X_pca"] = a.X
sc.pp.neighbors(a, use_rep="X_pca", n_neighbors=15)
sc.tl.diffmap(a)
meso = np.where(a.obs["cell_type"].to_numpy() == "Mesoderm")[0]
dc1 = a.obsm["X_diffmap"][:, 1]
a.uns["iroot"] = int(meso[np.argmin(dc1[meso])])          # extreme-DC1 Mesoderm cell as root
sc.tl.dpt(a)

# orient so pseudotime increases Mesoderm -> Myoblast -> Myocyte
mpt = a.obs.groupby("cell_type")["dpt_pseudotime"].mean()
print("mean pseudotime per state (pre-orient):\n" + mpt.to_string())
if mpt["Mesoderm"] > mpt["Myocyte"]:
    a.obs["dpt_pseudotime"] = a.obs["dpt_pseudotime"].max() - a.obs["dpt_pseudotime"]
    print("-> reversed pseudotime so Mesoderm is earliest")
print("mean pseudotime per state (final):\n" +
      a.obs.groupby("cell_type")["dpt_pseudotime"].mean().to_string())

# equal-COUNT (quantile) bins: the DPT pseudotime is heavily skewed (most cells at high pt), so equal-width
# bins leave the early bins nearly empty. qcut gives ~equal cells/bin -> robust per-bin means along the ranking.
a.obs["bin"] = pd.qcut(a.obs["dpt_pseudotime"], q=NBIN, labels=False, duplicates="drop")
print("\ncells per quantile bin:\n" + a.obs["bin"].value_counts().sort_index().to_string())
print("dominant state per bin:\n" +
      a.obs.groupby("bin")["cell_type"].agg(lambda s: s.value_counts().idxmax()).to_string())

# ---- write ----
rna_cols = [f"{g}_rna" for g in GENES if f"{g}_rna" in a.obs.columns]
out = a.obs[["cell_type", "dpt_pseudotime", "bin"] + rna_cols].copy()
out.insert(0, "barcode_base", base)
out.to_csv(os.path.join(HERE, "pseudotime_table.csv"))
pd.Series(base).to_csv(os.path.join(HERE, "cells_base.txt"), index=False, header=False)
print(f"\nwrote pseudotime_table.csv ({out.shape[0]} cells, {NBIN} bins) + cells_base.txt")
