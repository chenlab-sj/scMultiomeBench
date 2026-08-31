# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
RMS Mast607A reproducibility (Fig4A/B), 14 methods x 3 seed-replicate runs.
Port of BRCA/benchmark/fig4/reproduce_metrics.py -- reads the STAGED latents
(RMS/benchmark/staged/<method>/rep{1,2,3}.csv, already normalized to <bc>-1_rna / <bc>-1_atac by
01_prep_latents.py), so no per-method barcode handling is needed here.
  Fig4A = per-cell-type std of the KNN ATAC-label prediction across the 3 reps -> reproducibility = 1 - Avg_SD
  Fig4B = pairwise NMI of louvain clusterings across the 3 reps
RMS label filter = modality.startswith("test") (train rows are "train multiomics (Mast39)"/"(Mast213F)").
Outputs (consumed by plot_fig4.py): rep_knn_k10_pred_label.csv, louvain_cluster_reproduce.csv,
sd_df.csv, nmi_df.csv, reproducibility.csv.
Run in benchmark_env with BENCHMARK_FUN_DIR set (CPU; the louvain clusterings are the slow part).
"""
import os
import re
import sys
import glob
import numpy as np
import pandas as pd
import scanpy as sc
from functools import reduce
from itertools import combinations
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, normalized_mutual_info_score

sys.path.append(os.environ.get("BENCHMARK_FUN_DIR", "/path/to/multiomeBench/common"))
import benchmark_fun

HERE = os.path.dirname(os.path.abspath(__file__))
STAGE = os.path.normpath(os.path.join(HERE, "..", "staged"))
OUT = os.environ.get("REPRO_OUT", ".")
LABEL = os.environ.get("LABEL", "label.csv")
K, SEED = 10, 420
MIN_CELLS = int(os.environ.get("MIN_CELLS", "100"))

METHODS = ["scVI", "Cobolt", "Portal", "scglue", "scglue(multiome)", "BindSC", "Seurat(CCA)",
           "MaxFuse", "MIDAS", "scDART", "scBridge", "simba", "scButterfly", "scJoint"]

labels_annot = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test = labels_annot.index.values[labels_annot["modality"].astype(str).str.startswith("test")]
bc_test_rna = [b for b in bc_test if b.endswith("_rna")]
bc_test_atac = [b for b in bc_test if b.endswith("_atac")]

# load the staged latents; name them <method>-1/2/3 and remember the per-method grouping.
# A method missing any rep (e.g. simba/scButterfly not finished) is SKIPPED so the rest still run.
lat_dfs, pipelines, pipeline_groups = [], [], {}
for method in METHODS:
    # use ALL staged seeds per method (some have 5, most 3); reproducibility is computed within each
    # method over ITS OWN reps, so unequal rep counts are valid (SD over N reps; NMI over C(N,2) pairs).
    paths = sorted(glob.glob(os.path.join(STAGE, method, "rep*.csv")),
                   key=lambda p: int(re.search(r"rep(\d+)\.csv$", p).group(1)))
    if len(paths) < 2:
        print(f"  WARNING: {method} has <2 staged reps -> skipped (re-run 01_prep_latents.py when ready)")
        continue
    reps = []
    for p in paths:
        r = int(re.search(r"rep(\d+)\.csv$", p).group(1))
        name = f"{method}-{r}"
        lat_dfs.append(pd.read_csv(p, index_col=0))
        pipelines.append(name)
        reps.append(name)
    pipeline_groups[method] = reps
print("loaded %d latents; reps/method: %s" % (len(lat_dfs), {m: len(v) for m, v in pipeline_groups.items()}))

# ---- KNN ATAC-label transfer (k=10, cosine, distance-weighted), one column per rep ----
df_preds = []
for lat_df, name in zip(lat_dfs, pipelines):
    rna = [b for b in bc_test_rna if b in lat_df.index]
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
# aggregate std ONLY over MAJOR cell types (>MIN_CELLS). The confusion matrix still spans all types
# (a rare type can appear as a predicted column), but rare-type ROWS do not enter the reproducibility mean.
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
#   reflects EMBEDDING stability rather than per-rep resolution-search noise (find_louvain_res binary-
#   searches a resolution separately per rep; on continuous / few-cluster data that is a big NMI-noise
#   source and is what makes 4B disagree with the 4A label-transfer panel). FIX_RES=0 (default) keeps
#   the original per-rep search.
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
clusters = reduce(lambda l, r: l.join(r, how="inner"), df_clusters)   # cells common to all reps
clusters.to_csv(os.path.join(OUT, "louvain_cluster_reproduce.csv"))

nmi_rows = []
for method, reps in pipeline_groups.items():
    for a, b in combinations(reps, 2):
        nmi_rows.append({"Method": method, "Pair": f"{a} vs {b}",
                         "NMI": float(normalized_mutual_info_score(clusters[a], clusters[b]))})
pd.DataFrame(nmi_rows).to_csv(os.path.join(OUT, "nmi_df.csv"), index=False)

print("DONE: reproducibility.csv, sd_df.csv, nmi_df.csv, rep_knn_k10_pred_label.csv, louvain_cluster_reproduce.csv")
