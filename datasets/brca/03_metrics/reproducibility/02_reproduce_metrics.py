# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
BRCA HT243B1-S1H4 reproducibility (Fig4A/B), 14 methods x 3 seed-replicate runs.
Parameterized port of HT243_S1H4_reproduce.ipynb (the deep-dive top-performer set + the 3 new methods):
  Fig4A = per-cell-type std of the KNN ATAC-label prediction across the 3 reps
          -> reproducibility = 1 - Avg_SD
  Fig4B = pairwise NMI of louvain clusterings across the 3 reps
scVI rep1 path is FIXED to HT243B1-S1H4 (the notebook had HT514B1-S1H3 = a different patient).
Outputs (consumed by plot_fig4.py): rep_knn_k10_pred_label.csv, louvain_cluster_reproduce.csv,
sd_df.csv, nmi_df.csv, reproducibility.csv.
Run in benchmark_env with BENCHMARK_FUN_DIR set (CPU; the 42 louvain clusterings are the slow part).
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from functools import reduce
from itertools import combinations
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, normalized_mutual_info_score

sys.path.append(os.environ.get("BENCHMARK_FUN_DIR", "/path/to/multiomeBench/common"))
import benchmark_fun

HTAN = "/path/to/data/HTAN/HT243B1-S1H4"
REPO = "/path/to/multiomeBench/BRCA/HT243B1-S1H4"
OUT   = os.environ.get("REPRO_OUT", ".")
LABEL = os.environ.get("LABEL", "label.csv")
K, SEED = 10, 420

# method -> [rep1, rep2, rep3] latent paths (3 seed-replicate runs each); paths from the reproduce notebook
METHODS = {
    "scDART":           [f"{HTAN}/scDART/latent.csv",             f"{HTAN}/scDART/rep2/latent.csv",             f"{HTAN}/scDART/rep3/latent.csv"],
    "Seurat(CCA)":      [f"{HTAN}/Seurat3/coembed_coor.csv",      f"{HTAN}/Seurat3/rep2/coembed_coor.csv",      f"{HTAN}/Seurat3/rep3/coembed_coor.csv"],
    "BindSC":           [f"{HTAN}/bindSC/coembed_coor.csv",       f"{HTAN}/bindSC/rep2/coembed_coor.csv",       f"{HTAN}/bindSC/rep3/coembed_coor.csv"],
    "scBridge":         [f"{HTAN}/scBridge/latent.csv",           f"{HTAN}/scBridge/rep2/latent.csv",           f"{HTAN}/scBridge/rep3/latent.csv"],
    "scglue":           [f"{HTAN}/scglue/latent.csv",             f"{HTAN}/scglue/rep2/latent.csv",             f"{HTAN}/scglue/rep3/latent.csv"],
    "scglue(multiome)": [f"{HTAN}/scglue_paired/scglue_latent.csv", f"{HTAN}/scglue_paired/rep2/scglue_latent.csv", f"{HTAN}/scglue_paired/rep3/scglue_latent.csv"],
    "scJoint":          [f"{HTAN}/scJoint/latent.csv",            f"{HTAN}/scJoint/rep2/latent.csv",            f"{HTAN}/scJoint/rep3/latent.csv"],
    "scVI":             [f"{HTAN}/scVI/scvi_latent.csv",          f"{HTAN}/scVI/rep2/scvi_latent.csv",          f"{HTAN}/scVI/rep3/scvi_latent.csv"],  # rep1 FIXED (was HT514B1-S1H3)
    "Cobolt":           [f"{HTAN}/cobolt/cobolt_latent.csv",      f"{HTAN}/cobolt/rep2/cobolt_latent.csv",      f"{HTAN}/cobolt/rep3/cobolt_latent.csv"],
    "simba":            [f"{HTAN}/simba/latent.csv",              f"{HTAN}/simba/rep2/latent.csv",              f"{HTAN}/simba/rep3/latent.csv"],
    "Portal":           [f"{HTAN}/portal/lat_df.csv",             f"{HTAN}/portal/rep2/lat_df.csv",             f"{HTAN}/portal/rep3/lat_df.csv"],
    "MaxFuse":          [f"{REPO}/maxfuse/latent.csv",            f"{REPO}/maxfuse/rep2/latent.csv",            f"{REPO}/maxfuse/rep3/latent.csv"],
    "MIDAS":            [f"{REPO}/midas/latent.csv",              f"{REPO}/midas/rep2/latent.csv",              f"{REPO}/midas/rep3/latent.csv",              f"{REPO}/midas/rep4/latent.csv",              f"{REPO}/midas/rep5/latent.csv"],
    "scButterfly":      [f"{REPO}/scbutterfly/latent.csv",        f"{REPO}/scbutterfly/rep4/latent.csv",        f"{REPO}/scbutterfly/rep5/latent.csv"],  # forced-seed triplet (rep2/3 byte-identical -> artificial 1.0; rep4/5 FORCE_SEED 100/200)
}

labels_annot = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test = labels_annot.index.values[~labels_annot["modality"].isin(["train multiomics"])]
bc_test_rna  = [b for b in bc_test if b.endswith("_rna")]
bc_test_atac = [b for b in bc_test if b.endswith("_atac")]

# load the 42 latents; name them <method>-1/2/3 and remember the per-method grouping
lat_dfs, pipelines, pipeline_groups = [], [], {}
for method, paths in METHODS.items():
    reps = []
    for r, p in enumerate(paths, 1):
        name = f"{method}-{r}"
        lat_dfs.append(pd.read_csv(p, index_col=0))
        pipelines.append(name); reps.append(name)
    pipeline_groups[method] = reps
print(f"loaded {len(lat_dfs)} latents ({len(pipeline_groups)} methods x 3 reps)")

# ---- KNN ATAC-label transfer (k=10, cosine, distance-weighted), one column per rep ----
df_preds = []
for lat_df, name in zip(lat_dfs, pipelines):
    rna  = [b for b in bc_test_rna  if b in lat_df.index]
    atac = [b for b in bc_test_atac if b in lat_df.index]
    clf = KNeighborsClassifier(n_neighbors=K, metric="cosine", weights="distance")
    clf.fit(lat_df.loc[rna], labels_annot.loc[rna, "cell_type"])
    pred = clf.predict(lat_df.loc[atac])
    df_preds.append(pd.DataFrame({name: pred}, index=atac))
KNN_atac = reduce(lambda l, r: l.join(r, how="inner"), df_preds)   # atac cells common to all reps
KNN_atac.to_csv(os.path.join(OUT, "rep_knn_k10_pred_label.csv"))

Actual = labels_annot.loc[KNN_atac.index, "cell_type"].values
counts = pd.Series(Actual).value_counts()
sorted_categories = counts.index.tolist()                # all types, by frequency (confusion-matrix columns)
# MATCH the original Fig4: aggregate std ONLY over MAJOR cell types (>100 cells). The confusion matrix
# still spans all types (so a rare type can appear as a predicted column), but rare-type ROWS do not
# enter the reproducibility mean -- exactly what the notebook's sorted_categories_label (7) did.
MIN_CELLS = int(os.environ.get("MIN_CELLS", "100"))
major = [c for c in sorted_categories if counts[c] > MIN_CELLS]
print(f"cell types: {len(sorted_categories)} total; {len(major)} major (>{MIN_CELLS}) used for Fig4A: {major}")

# ---- Fig4A: per-cell-type std of the row-normalized confusion across the 3 reps (major types only) ----
sd_data = []
for method, reps in pipeline_groups.items():
    cms = []
    for rep in reps:
        cm = confusion_matrix(Actual, KNN_atac[rep].values, labels=sorted_categories).astype(float)
        cms.append(cm / cm.sum(axis=1, keepdims=True))
    cm_std = np.std(np.stack(cms, 0), axis=0)
    for i, ct in enumerate(sorted_categories):
        if ct not in major:                              # rare types excluded from the reproducibility mean
            continue
        sd_data.append({"Method": method, "CellType": ct, "Std": float(np.mean(cm_std[i, :]))})
sd_df = pd.DataFrame(sd_data)
sd_df.to_csv(os.path.join(OUT, "sd_df.csv"), index=False)
repro = sd_df.groupby("Method")["Std"].mean().reset_index().rename(columns={"Std": "Avg_SD"})
repro["reproducibility"] = 1 - repro["Avg_SD"]
repro.sort_values("Avg_SD").to_csv(os.path.join(OUT, "reproducibility.csv"), index=False)

# ---- louvain cluster each rep -> Fig4B pairwise NMI ----
# FIX_RES=1: derive ONE louvain resolution per method from its rep1 and reuse it for every rep, so NMI
#   reflects EMBEDDING stability rather than per-rep resolution-search noise (the main 4A/4B-discordance
#   driver on continuous / few-cluster data). FIX_RES=0 = original per-rep search.
FIX_RES = os.environ.get("FIX_RES", "0") == "1"
nclust = len(major)
bc_test_nmi = labels_annot.index.values[
    labels_annot.index.isin(bc_test) & labels_annot["cell_type"].isin(major)
]
lat_by_name = dict(zip(pipelines, lat_dfs))

def _emb(lat_df):
    lt = lat_df.loc[lat_df.index.isin(bc_test_nmi)]
    a = sc.AnnData(lt); a.obsm["X_emb"] = lt
    return lt, a

df_clusters = []
for method, reps in pipeline_groups.items():
    method_res = benchmark_fun.find_louvain_res(_emb(lat_by_name[reps[0]])[1], nclust) if FIX_RES else None
    if FIX_RES:
        print(f"{method}: fixed louvain res {method_res:.4f} from {reps[0]} -> all {len(reps)} reps", flush=True)
    for name in reps:
        lt, a = _emb(lat_by_name[name])
        res = method_res if FIX_RES else benchmark_fun.find_louvain_res(a, nclust)
        app = a.copy(); sc.pp.neighbors(app); sc.tl.louvain(app, resolution=res, random_state=SEED)
        df_clusters.append(pd.DataFrame({name: app.obs["louvain"].astype(str).values}, index=lt.index))
clusters = reduce(lambda l, r: l.join(r, how="inner"), df_clusters)   # cells common to all reps
clusters.to_csv(os.path.join(OUT, "louvain_cluster_reproduce.csv"))

nmi_rows = []
for method, reps in pipeline_groups.items():
    for a, b in combinations(reps, 2):
        nmi_rows.append({"Method": method, "Pair": f"{a} vs {b}",
                         "NMI": float(normalized_mutual_info_score(clusters[a], clusters[b]))})
pd.DataFrame(nmi_rows).to_csv(os.path.join(OUT, "nmi_df.csv"), index=False)

print("DONE: reproducibility.csv, sd_df.csv, nmi_df.csv, rep_knn_k10_pred_label.csv, louvain_cluster_reproduce.csv")
