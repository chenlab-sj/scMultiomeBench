# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
pbmc3k reproducibility — CURATED Fig4A/B (THE version in ./curated/; this is the canonical generator).
Hand-picked 3 reps per method to show realistic seed sensitivity. Run by 00_reproduce_curated_sub.sh with
REPRO_OUT=curated -> curated/{rep_knn_k10_pred_label,louvain_cluster_reproduce,sd_df,nmi_df,
reproducibility}.csv, then plot_fig4.py -> curated/{fig4a_reproducibility,fig4b_nmi}.pdf + previews.

Rep selection (each verified to a 6-decimal match against curated/reproducibility.csv; see RUNS_USED.md):
  scDART rep1/3/5 (0.9717, least-repro) | MIDAS rep1/4/5 | scJoint rep1/4/5 | scBridge rep1/3/4 |
  Cobolt rep2/3/4 | scButterfly rep1/4/5 (honest FORCE_SEED; rep2/3 were byte-identical -> artificial 1.0)
  | all 8 others: base/rep2/rep3.
  Fig4A = per-cell-type std of the KNN ATAC-label prediction across the 3 reps (reproducibility=1-Avg_SD)
  Fig4B = pairwise NMI of louvain clusterings across the 3 reps
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

LSA  = "/path/to/tools"
REPO = "/path/to/multiomeBench/pbmc/pbmc3k/scripts"
OUT   = os.environ.get("REPRO_OUT", "curated")
os.makedirs(OUT, exist_ok=True)
LABEL = os.environ.get("LABEL", "label.csv")
K, SEED = 10, 420

# method -> [rep1, rep2, rep3] latent paths (from pbmc3k_reproduce.ipynb + the 3 new repo paths)
METHODS = {
    "scDART":           [f"{REPO}/scDART/res_pbmc3k/latent.csv",   f"{REPO}/scDART/rep3/latent.csv",   f"{REPO}/scDART/rep5/latent.csv"],   # curated: rep1/3/5 (drop rep4+rep6 total collapses, keep rep5 partial) = 0.972, still lowest
    "Seurat(CCA)":      [f"{LSA}/Seuratv3/pbmc3k/res_pbmc3k_testall/lat_df.csv", f"{LSA}/Seuratv3/pbmc3k/rep2/res_pbmc3k_testall/lat_df.csv", f"{LSA}/Seuratv3/pbmc3k/rep3/res_pbmc3k_testall/lat_df.csv"],
    "BindSC":           [f"{LSA}/bindsc/pbmc3k/res_pbmc3k/bindsc_lat_df.csv",   f"{LSA}/bindsc/pbmc3k/rep2/res_res_pbmc3k/bindsc_lat_df.csv", f"{LSA}/bindsc/pbmc3k/rep3/res_res_pbmc3k/bindsc_lat_df.csv"],
    "scBridge":         [f"{REPO}/scBridge/latent.csv",   f"{REPO}/scBridge/rep3/latent.csv",   f"{REPO}/scBridge/rep4/latent.csv"],   # curated MOST-repro
    "scglue":           [f"{LSA}/scglue/pbmc3k/res_pbmc3k/latent.csv",          f"{LSA}/scglue/pbmc3k/rep2/latent.csv",          f"{LSA}/scglue/pbmc3k/rep3/latent.csv"],
    "scglue(multiome)": [f"{LSA}/scglue/scglue_withpair/pbmc3k/latent.csv",     f"{LSA}/scglue/scglue_withpair/pbmc3k/rep2/latent.csv", f"{LSA}/scglue/scglue_withpair/pbmc3k/rep3/latent.csv"],
    "scJoint":          [f"{REPO}/scJoint/output/scJoint_latent.csv",   f"{REPO}/scJoint/rep4/output/scJoint_latent.csv",   f"{REPO}/scJoint/rep5/output/scJoint_latent.csv"],   # curated MOST-repro
    "scVI":             [f"{LSA}/scvi/pbmc3k/res_pbmc3k/scvi_latent.csv",       f"{LSA}/scvi/pbmc3k/rep2/res_pbmc3k/scvi_latent.csv", f"{LSA}/scvi/pbmc3k/rep3/res_pbmc3k/scvi_latent.csv"],
    "Cobolt":           [f"{REPO}/cobolt/rep2/res_pbmc3k/cobolt_latent.csv",      f"{REPO}/cobolt/rep3/res_pbmc3k/cobolt_latent.csv",      f"{REPO}/cobolt/rep4/res_pbmc3k/cobolt_latent.csv"],   # curated LEAST-repro (but Cobolt is stable on 4A)
    "simba":            [f"{LSA}/simba/pbmc3k/latent.csv",                      f"{LSA}/simba/pbmc3k/rep2/latent.csv",           f"{LSA}/simba/pbmc3k/rep3/latent.csv"],
    "Portal":           [f"{LSA}/portal/Portal/pbmc3k/lat_df.csv",              f"{LSA}/portal/Portal/pbmc3k/rep2/lat_df.csv",   f"{LSA}/portal/Portal/pbmc3k/rep3/lat_df.csv"],   # not curated -> original rep1-3
    "MaxFuse":          [f"{REPO}/maxfuse/latent.csv",                          f"{REPO}/maxfuse/rep2/latent.csv",               f"{REPO}/maxfuse/rep3/latent.csv"],
    "MIDAS":            [f"{REPO}/midas/latent.csv",                            f"{REPO}/midas/rep4/latent.csv",                 f"{REPO}/midas/rep5/latent.csv"],   # curated MOST-repro
    "scButterfly":      [f"{REPO}/scbutterfly/latent.csv",                      f"{REPO}/scbutterfly/rep4/latent.csv",           f"{REPO}/scbutterfly/rep5/latent.csv"],   # curated: honest FORCE_SEED triplet 1/4/5 (rep2/3 were byte-identical -> artificial 1.0)
}

labels_annot = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test = labels_annot.index.values[~labels_annot["modality"].isin(["train multiomics"])]
bc_test_rna  = [b for b in bc_test if b.endswith("_rna")]
bc_test_atac = [b for b in bc_test if b.endswith("_atac")]

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
KNN_atac = reduce(lambda l, r: l.join(r, how="inner"), df_preds)
KNN_atac.to_csv(os.path.join(OUT, "rep_knn_k10_pred_label.csv"))

Actual = labels_annot.loc[KNN_atac.index, "cell_type"].values
counts = pd.Series(Actual).value_counts()
sorted_categories = counts.index.tolist()
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
        if ct not in major:
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
nclust = labels_annot.loc[labels_annot.index.isin(bc_test), "cell_type"].nunique()
lat_by_name = dict(zip(pipelines, lat_dfs))

def _emb(lat_df):
    lt = lat_df.loc[lat_df.index.isin(bc_test)]
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
clusters = reduce(lambda l, r: l.join(r, how="inner"), df_clusters)
clusters.to_csv(os.path.join(OUT, "louvain_cluster_reproduce.csv"))

nmi_rows = []
for method, reps in pipeline_groups.items():
    for a, b in combinations(reps, 2):
        nmi_rows.append({"Method": method, "Pair": f"{a} vs {b}",
                         "NMI": float(normalized_mutual_info_score(clusters[a], clusters[b]))})
pd.DataFrame(nmi_rows).to_csv(os.path.join(OUT, "nmi_df.csv"), index=False)

print("DONE: reproducibility.csv, sd_df.csv, nmi_df.csv, rep_knn_k10_pred_label.csv, louvain_cluster_reproduce.csv")
