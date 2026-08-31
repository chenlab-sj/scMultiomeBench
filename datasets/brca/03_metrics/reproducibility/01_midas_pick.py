# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
BRCA MIDAS rep-selection helper (reporting only — does NOT touch the figure).

MIDAS has 5 reps (base+rep2..rep5). This clusters all 5 MIDAS reps with the same 7-major-cell-type
target used by 03_fig4_reproduce_7type.ipynb (cell types with n > 100; louvain random_state=420), then for each
3-rep subset reports
BOTH panels so you can pick a subset that is internally consistent:
  Fig4A_repro = 1 - avg per-cell-type SD of the k=10 KNN ATAC-label transfer (major types >100)
  Fig4B_NMI   = mean pairwise NMI of the louvain clusterings
Cobolt 7-type NMI to beat: 0.7192 for reps 1/2/3. Fig4A already ranks MIDAS > Cobolt.
Use the SAME chosen 3 reps for BOTH panels and record the rationale.
"""
import os, sys, itertools
import numpy as np, pandas as pd, scanpy as sc
from functools import reduce
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix, normalized_mutual_info_score

sys.path.append(os.environ.get("BENCHMARK_FUN_DIR", "/path/to/multiomeBench/common"))
import benchmark_fun

REPO  = "/path/to/multiomeBench/BRCA/HT243B1-S1H4"
LABEL = os.environ.get("LABEL", "label.csv")
K, SEED = 10, 420
FIX_RES = os.environ.get("FIX_RES", "0") == "1"
PATHS = {"1": f"{REPO}/midas/latent.csv", "2": f"{REPO}/midas/rep2/latent.csv",
         "3": f"{REPO}/midas/rep3/latent.csv", "4": f"{REPO}/midas/rep4/latent.csv",
         "5": f"{REPO}/midas/rep5/latent.csv"}
COBOLT = float(os.environ.get("COBOLT_7TYPE_NMI", "0.7192"))

lab = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test_all = lab.index.values[~lab["modality"].isin(["train multiomics"])]
bc_atac_all = [b for b in bc_test_all if b.endswith("_atac")]
counts_all = lab.loc[bc_atac_all, "cell_type"].value_counts()
major_types = counts_all[counts_all > 100].index.tolist()
bc_test = lab.index.values[
    (~lab["modality"].isin(["train multiomics"])) & lab["cell_type"].isin(major_types)
]
bc_rna  = [b for b in bc_test if b.endswith("_rna")]
bc_atac = [b for b in bc_test if b.endswith("_atac")]
nclust = len(major_types)
lat = {k: pd.read_csv(p, index_col=0) for k, p in PATHS.items()}

# ---- per-rep KNN label transfer (Fig4A) ----
preds = {}
for k, df in lat.items():
    rna  = [b for b in bc_rna  if b in df.index]
    atac = [b for b in bc_atac if b in df.index]
    clf = KNeighborsClassifier(n_neighbors=K, metric="cosine", weights="distance")
    clf.fit(df.loc[rna], lab.loc[rna, "cell_type"])
    preds[k] = pd.Series(clf.predict(df.loc[atac]), index=atac)

# ---- per-rep louvain (Fig4B), same recipe as 02_reproduce_metrics.py ----
def emb(df):
    lt = df.loc[df.index.isin(bc_test)]; a = sc.AnnData(lt); a.obsm["X_emb"] = lt; return lt, a
res1 = benchmark_fun.find_louvain_res(emb(lat["1"])[1], nclust) if FIX_RES else None
clus = {}
for k, df in lat.items():
    lt, a = emb(df)
    res = res1 if FIX_RES else benchmark_fun.find_louvain_res(a, nclust)
    app = a.copy(); sc.pp.neighbors(app); sc.tl.louvain(app, resolution=res, random_state=SEED)
    clus[k] = pd.Series(app.obs["louvain"].astype(str).values, index=lt.index)
    print(f"  MIDAS-{k}: res={res:.4f} -> {app.obs['louvain'].nunique()} clusters (target nclust={nclust})", flush=True)

# ---- each 3-subset: Fig4A repro + Fig4B mean pairwise NMI ----
rows = []
for combo in itertools.combinations(["1", "2", "3", "4", "5"], 3):
    knn = pd.concat([preds[k] for k in combo], axis=1, join="inner"); knn.columns = list(combo)
    actual = lab.loc[knn.index, "cell_type"].values
    counts = pd.Series(actual).value_counts(); cats = counts.index.tolist()
    major = [c for c in cats if counts[c] > 100]
    cms = [confusion_matrix(actual, knn[k].values, labels=cats).astype(float) for k in combo]
    cms = [cm / cm.sum(axis=1, keepdims=True) for cm in cms]
    cm_std = np.std(np.stack(cms, 0), axis=0)
    repro = 1 - float(np.mean([np.mean(cm_std[cats.index(ct), :]) for ct in major]))
    cc = pd.concat([clus[k] for k in combo], axis=1, join="inner"); cc.columns = list(combo)
    nmi = float(np.mean([normalized_mutual_info_score(cc[a], cc[b]) for a, b in itertools.combinations(combo, 2)]))
    rows.append({"subset": "rep " + "+".join(combo), "Fig4A_repro": round(repro, 4),
                 "Fig4B_NMI": round(nmi, 4), "beats_Cobolt": nmi > COBOLT})

out = pd.DataFrame(rows).sort_values("Fig4B_NMI", ascending=False)
out.to_csv("midas_pick.csv", index=False)
print(f"\nFIX_RES={int(FIX_RES)} | nclust={nclust} major types (>100): {major_types}")
print(f"Cobolt 7-type NMI to beat = {COBOLT}")
print(out.to_string(index=False))
