#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
FigS1B k-test = pbmc10k KNN RNA-label-prediction accuracy vs k (the REVERSE of FigS1A): fit KNN on the test
ATAC cells, predict the test RNA cells' cell type, per-cell-type accuracy, k in {5,10,20,40,80}. Same
KNeighborsClassifier(n_neighbors=k, metric='cosine', weights='distance') + class-wise-accuracy recipe as
00_compute_ktest_legacy18.py::knn_ataclabel, just with RNA<->ATAC swapped. No published reverse data exists, so this
recomputes ALL 23 methods (18 original latents on the cluster + 5 new in the repo) -> figS1b_ktest_long.csv.

VALIDATE: also runs the FORWARD direction for a few old methods and compares to the published FigS1A
(kNN_celltype_accu_sum.csv) so we know the recipe/barcodes reproduce it before trusting the reverse.
Run in an env with sklearn/pandas (benchmark_env on cluster, or Mac).  python 04_compute_ktest_reverse.py
"""
import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, confusion_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
G = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
LSA = "/path/to/tools" if os.path.exists("/path/to/tools") else "/path/to/tools"
LABEL = os.path.join(HERE, "label.csv")
KS = [5, 10, 20, 40, 80]

# 23 methods -> latent (18 original on the cluster; 5 new in the repo) -- same map as 00_compute_metrics_sub.sh
LAT = {
    "BindSC": f"{LSA}/bindsc/pbmc10k/res_pbmc10k_control/bindsc_lat_df.csv",
    "Conos": f"{LSA}/Concos/pbmc10k/res_pbmc10k_control/coembed_coor.csv",
    "LIGER": f"{LSA}/liger/pbmc10k/res_pbmc10k_control/coembed_coor.csv",
    "MultiMAP": f"{LSA}/multimap/pbmc10k/res_pbmc10k_control/coembed_coor.csv",
    "Seurat(CCA)": f"{LSA}/Seuratv3/pbmc10k/res_pbmc10k_control/lat_df.csv",
    "scDART": f"{LSA}/scDART/scDART/pbmc10k/res_pbmc10k/latent.csv",
    "scMoMaT": f"{LSA}/scMoMaT/pbmc10k/latent.csv",
    "simba": f"{LSA}/simba/pbmc10k/latent.csv",
    "scglue": f"{LSA}/scglue/pbmc10k/res_pbmc10K_control/latent.csv",
    "Unioncom": f"{LSA}/unioncom/pbmc10k/latent.csv",
    "Portal": f"{LSA}/portal/Portal/pbmc10k/lat_df.csv",
    "scJoint": f"{LSA}/scJoint/pbmc10k/res_control/scjoint_lat_df.csv",
    "scBridge": f"{LSA}/scBridge/pbmc10k/latent.csv",
    "Cobolt": f"{LSA}/cobolt/pbmc10k/res_pbmc10k_test_all/cobolt_latent.csv",
    "scglue(multiome)": f"{LSA}/scglue/pbmc10k/res_pbmc10K_withpair_control/scglue_latent.csv",
    "scVI": f"{LSA}/scvi/pbmc10k/res_pbmc10k_test_all/scvi_latent.csv",
    "Seurat(WNN)": f"{LSA}/Seuratv4/pbmc10k/res_pbmc10k_test_all/seurat4_latent.csv",
    "MinNet": f"{LSA}/MinNet/pbmc10k_control/MinNet_lat_df.csv",
    "MaxFuse": f"{G}/pbmc/pbmc10k/scripts/maxfuse/latent.csv",
    "MIDAS": f"{G}/pbmc/pbmc10k/scripts/midas/latent.csv",
    "scButterfly": f"{G}/pbmc/pbmc10k/scripts/scbutterfly/latent.csv",
    "MIRA": f"{G}/pbmc/pbmc10k/scripts/MIRA/latent.csv",
    "Multigrate": f"{G}/pbmc/pbmc10k/scripts/multigrate/latent.csv",
}

lab = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
bc_test = lab[lab["modality"] != "train multiomics"].index.values
bc_rna = [b for b in bc_test if b.endswith("_rna")]
bc_atac = [b for b in bc_test if b.endswith("_atac")]

def knn_accu(k, lat, fit_bc, pred_bc):
    """VERBATIM 00_compute_ktest_legacy18.py::knn_ataclabel recipe (fit_bc -> pred_bc). Returns [pl?,acc,type] rows."""
    fit = [i for i in fit_bc if i in lat.index]
    prd = [i for i in pred_bc if i in lat.index]
    clf = KNeighborsClassifier(n_neighbors=k, metric="cosine", weights="distance")
    clf.fit(lat.loc[fit], lab.loc[fit, "cell_type"])
    pred = clf.predict(lat.loc[prd])
    true = lab.loc[prd, "cell_type"].values
    out = [(accuracy_score(true, pred), "overall")]
    conf = confusion_matrix(true, pred)                 # default labels = sorted(unique(true u pred))
    cwa = np.diagonal(conf) / np.sum(conf, axis=1)
    types = np.unique(pred)
    for j in range(len(types)):
        out.append((cwa[j], types[j]))
    return out

records, fwd_check = [], []
for method, path in LAT.items():
    if not os.path.exists(path):
        print(f"SKIP {method}: no latent ({path})"); continue
    lat = pd.read_csv(path, index_col=0)
    for k in KS:
        for acc, typ in knn_accu(k, lat, bc_atac, bc_rna):        # REVERSE: fit ATAC -> predict RNA
            records.append({"pipeline": method, "accuracy": acc, "type": typ, "k value": k})
    # forward spot-check (k=10) for validating against published FigS1A
    for acc, typ in knn_accu(10, lat, bc_rna, bc_atac):
        if typ != "overall":
            fwd_check.append({"pipeline": method, "type": typ, "fwd_accu_k10": acc})
    print(f"OK   {method}")

out = pd.DataFrame(records)[["pipeline", "accuracy", "type", "k value"]]
out.to_csv(os.path.join(HERE, "figS1b_ktest_long.csv"), index=False)
print(f"\nwrote figS1b_ktest_long.csv  ({out.pipeline.nunique()} methods x {len(KS)} k)")

# ---- validate forward vs published FigS1A (kNN_celltype_accu_sum.csv) ----
try:
    pub = pd.read_csv(os.path.join(HERE, "old/knn_test/kNN_celltype_accu_sum.csv"))
    pub10 = pub[pub["k value"] == 10]
    fc = pd.DataFrame(fwd_check)
    diffs = []
    for m in fc.pipeline.unique():
        if m not in set(pub10.pipeline):
            continue
        pr = pub10[pub10.pipeline == m].iloc[0]
        for _, r in fc[fc.pipeline == m].iterrows():
            if r["type"] in pr.index:
                diffs.append(abs(float(pr[r["type"]]) - r["fwd_accu_k10"]))
    if diffs:
        print(f"FORWARD validation vs published FigS1A (k=10, per-celltype): mean|diff|={np.mean(diffs):.3f} "
              f"max={np.max(diffs):.3f}  (small -> recipe/barcodes reproduce FigS1A, reverse trustworthy)")
except Exception as e:
    print("forward validation skipped:", e)
