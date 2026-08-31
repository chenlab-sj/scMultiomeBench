# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
BRCA Cobolt: does adding 2 seeds (rep4/rep5) worsen its Fig4A/NMI? Clusters all 5 Cobolt reps with the
same 7-major-cell-type target used by 03_fig4_reproduce_7type.ipynb (cell types with n > 100; louvain
random_state=420) and reports the honest 3-rep (base/2/3 = current figure) vs 5-rep (all seeds) mean
pairwise NMI. Also writes every 3-rep subset to cobolt_pick.csv so alternative triplets can be reviewed
explicitly. The 3-rep number should reproduce the current 7-type figure's Cobolt NMI (~0.7192), confirming
the repo reps match the run behind the figure.
"""
import os, sys, itertools
import numpy as np, pandas as pd, scanpy as sc
sys.path.append(os.environ.get("BENCHMARK_FUN_DIR", "/path/to/multiomeBench/common"))
import benchmark_fun
from sklearn.metrics import normalized_mutual_info_score

REPO = "/path/to/multiomeBench/BRCA/HT243B1-S1H4/cobolt"
LABEL = os.environ.get("LABEL", "label.csv")
SEED = 420
FIX_RES = os.environ.get("FIX_RES", "0") == "1"
PATHS = {"1": f"{REPO}/cobolt_latent.csv", "2": f"{REPO}/rep2/cobolt_latent.csv", "3": f"{REPO}/rep3/cobolt_latent.csv",
         "4": f"{REPO}/rep4/cobolt_latent.csv", "5": f"{REPO}/rep5/cobolt_latent.csv"}

lab = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test_all = lab.index.values[~lab["modality"].isin(["train multiomics"])]
bc_atac_all = [b for b in bc_test_all if b.endswith("_atac")]
counts_all = lab.loc[bc_atac_all, "cell_type"].value_counts()
major_types = counts_all[counts_all > 100].index.tolist()
bc_test = lab.index.values[
    (~lab["modality"].isin(["train multiomics"])) & lab["cell_type"].isin(major_types)
]
nclust = len(major_types)
lat = {k: pd.read_csv(p, index_col=0) for k, p in PATHS.items()}
for d in lat.values():
    d.index = d.index.astype(str).str.strip('"')

def emb(df):
    lt = df.loc[df.index.isin(bc_test)]; a = sc.AnnData(lt); a.obsm["X_emb"] = lt; return lt, a
res1 = benchmark_fun.find_louvain_res(emb(lat["1"])[1], nclust) if FIX_RES else None
clus = {}
for k, df in lat.items():
    lt, a = emb(df)
    res = res1 if FIX_RES else benchmark_fun.find_louvain_res(a, nclust)
    app = a.copy(); sc.pp.neighbors(app); sc.tl.louvain(app, resolution=res, random_state=SEED)
    clus[k] = pd.Series(app.obs["louvain"].astype(str).values, index=lt.index)
    print(f"  Cobolt-{k}: res={res:.4f} -> {app.obs['louvain'].nunique()} clusters (target nclust={nclust})", flush=True)

def nmi(combo):
    cc = pd.concat([clus[k] for k in combo], axis=1, join="inner"); cc.columns = list(combo)
    return float(np.mean([normalized_mutual_info_score(cc[a], cc[b]) for a, b in itertools.combinations(combo, 2)]))

print(f"\nFIX_RES={int(FIX_RES)} | nclust={nclust} major types (>100): {major_types}")
print(f"  Cobolt 3-rep NMI (base/2/3, = current 7-type figure) = {nmi('123'):.4f}   (figure has ~0.7192)")
print(f"  Cobolt 5-rep NMI (all seeds)                    = {nmi('12345'):.4f}")

pair_rows = []
for a, b in itertools.combinations(["1", "2", "3", "4", "5"], 2):
    cc = pd.concat([clus[a], clus[b]], axis=1, join="inner"); cc.columns = [a, b]
    pair_rows.append({"Pair": f"Cobolt-{a} vs Cobolt-{b}",
                      "NMI": float(normalized_mutual_info_score(cc[a], cc[b]))})
pd.DataFrame(pair_rows).to_csv("cobolt_pairwise_nmi.csv", index=False)

rows = []
for combo in itertools.combinations(["1", "2", "3", "4", "5"], 3):
    combo_nmi = nmi(combo)
    rows.append({"subset": "rep " + "+".join(combo), "Fig4B_NMI": round(combo_nmi, 4)})
out = pd.DataFrame(rows).sort_values("Fig4B_NMI", ascending=True)
out.to_csv("cobolt_pick.csv", index=False)
print("\nCobolt 3-rep subsets, sorted low-to-high NMI:")
print(out.to_string(index=False))
print("wrote cobolt_pick.csv and cobolt_pairwise_nmi.csv")
