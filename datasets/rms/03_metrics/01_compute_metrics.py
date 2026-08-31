#!/usr/bin/env python
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
00_compute_metrics.py -- per-method integration metrics for the pbmc3k Fig2B matrix.

CLEAN, PARAMETERIZED replacement for benchmark_metrics.py + 00_compute_ktest_legacy18.py.
One method per call -> APPENDS a row to the shared output CSVs, so adding a method is:
    python 00_compute_metrics.py --latent scripts/MIRA/latent.csv --method MIRA

It keeps the EXACT metric definitions (the benchmark_fun.py helpers + the same scib/KS
calls) so new methods stay directly comparable to the original 18. Differences vs the old
scripts are only organizational:
  * parameterized (no hardcoded 18-method list)
  * cell_test ("paired" cells for ARI/AMI) defined directly as test stems present as BOTH
    <stem>_rna and <stem>_atac -- NOT via Conos's latent (the old legacy dependency)
  * KNN label transfer done ONCE at k=10 (the only k Fig2B uses), not the {5,10,20,40,80} loop
  * dropped: the Cosinedisttop KS test, histogram grids, distance-distribution PDFs, dead code

Outputs (appended; one row/method, deduped on method) into --out:
  sum_metrics.csv      method, asw, asw_random, omics_asw, ami, ari, ks.statistic
  celltype_metrics.csv method, celltype, ks_celltype, ks_inter_celltype, celltype_omics_ASW
  knn_pred_accu.csv    method, accuracy, type            (type='overall' or a cell type; k=10)
  knn_pred_label.csv   per-method predicted ATAC labels  (index = ATAC barcode)
NOTE: the random-ADJUSTED accuracy (average_accu) and the peak score are added downstream
(see assemble step); this script writes the raw per-cell-type KNN accuracy they build on.
"""
import argparse
import os
import sys

# Import the scientific stack -- including the REAL umap-learn that scanpy's sc.pp.neighbors
# imports (`from umap.umap_ import ...`) -- BEFORE putting BENCHMARK_FUN_DIR on sys.path.
# That dir (the old benchmark scripts) contains a umap.py that would otherwise SHADOW
# umap-learn -> "ModuleNotFoundError: No module named 'umap.umap_'; 'umap' is not a package".
# Importing it first caches the real package in sys.modules so the shadow is ignored.
import numpy as np
import pandas as pd
import scanpy as sc
import umap.umap_  # noqa: F401  (force-load the real umap-learn before the shadow dir is added)
from scipy import stats
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
import scib

# benchmark_fun.py holds the exact metric helpers (lat2pdist, find_louvain_res,
# cluster_cosistent, get_parallel_dist, get_celltype_dist) -- keep using it verbatim.
BENCHMARK_FUN_DIR = os.environ.get("BENCHMARK_FUN_DIR", "/path/to/multiomeBench/common")
if BENCHMARK_FUN_DIR not in sys.path:
    sys.path.append(BENCHMARK_FUN_DIR)
import benchmark_fun

SEED = 420
OPTION = "Cosine"


def load_label_context(label_csv):
    """Parse label.csv into the test-cell context the metrics need."""
    lab = pd.read_csv(label_csv, index_col=0)
    # RMS: test rows are "test scRNA"/"test scATAC"; train rows are "train multiomics (Mast39)" /
    # "(Mast213F)" (suffixed) -> filter by startswith("test"), NOT != "train multiomics".
    bc_test = lab[lab["modality"].astype(str).str.startswith("test")].index.values
    labels_test = lab.loc[bc_test]
    bc_test_rna = [b for b in bc_test if b.endswith("_rna")]
    bc_test_atac = [b for b in bc_test if b.endswith("_atac")]
    label_dic = dict(zip(labels_test.index, labels_test["cell_type"]))
    celltype = labels_test["cell_type"].unique().tolist()
    # "paired" cells = stems present as BOTH _rna and _atac (clean; no Conos dependency)
    rna_stems = {b[:-4] for b in bc_test_rna}
    atac_stems = {b[:-5] for b in bc_test_atac}
    cell_test = list(rna_stems & atac_stems)
    return dict(labels_df=lab, bc_test=bc_test, bc_test_rna=bc_test_rna,
                bc_test_atac=bc_test_atac, label_dic=label_dic, celltype=celltype,
                cell_test=cell_test, nclust=len(lab["cell_type"].unique()))


def integration_metrics(lat_df, method, ctx):
    """The exact per-method body of the old benchmark_metrics() loop, for one method."""
    labels_df, celltype = ctx["labels_df"], ctx["celltype"]
    bc_test, bc_test_rna, bc_test_atac = ctx["bc_test"], ctx["bc_test_rna"], ctx["bc_test_atac"]
    label_dic, nclust = ctx["label_dic"], ctx["nclust"]

    bc_test_sel = list(set(lat_df.index).intersection(bc_test))
    # cell_test = "paired" cells among the TEST cells present in this method's latent (stems with
    # both <stem>_rna and <stem>_atac in bc_test_sel). Derived from bc_test_sel -- NOT the full
    # latent -- so every barcode it references is guaranteed to be in adata_pp.obs below. This is
    # what makes methods whose latent carries extra/non-test barcodes (e.g. Conos) not KeyError.
    rna_stems = {b[:-4] for b in bc_test_sel if b.endswith("_rna")}
    atac_stems = {b[:-5] for b in bc_test_sel if b.endswith("_atac")}
    cell_test = list(rna_stems & atac_stems)

    lat_df_test = lat_df.loc[lat_df.index.isin(bc_test_sel)]
    labels_test = labels_df.loc[bc_test_sel]
    adata = sc.AnnData(lat_df_test)
    adata.obsm["X_emb"] = lat_df_test
    adata.obs = pd.merge(adata.obs, labels_test, how="left", left_index=True, right_index=True)
    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"])
    adata.obs["random_celltype"] = pd.Categorical(adata.obs["random_celltype"])

    asw = scib.me.silhouette(adata, label_key="cell_type", embed="X_emb")
    asw_random = scib.me.silhouette(adata, label_key="random_celltype", embed="X_emb")
    adata.obs["modality"] = pd.Categorical(adata.obs["modality"])
    omics_asw = scib.me.silhouette_batch(adata, batch_key="modality", label_key="cell_type",
                                         embed="X_emb", return_all=True)
    omics_asw_ave, omics_asw_celltype = omics_asw[0], omics_asw[1]

    # ari/ami measure SAME-CELL cross-omics cluster consistency -> need paired cells (a stem present
    # as both <stem>_rna and <stem>_atac). The Parse cross-platform data is UNPAIRED (Parse RNA cells
    # != pbmc3k ATAC cells) -> cell_test is empty and these metrics are undefined (set NaN).
    if len(cell_test) == 0:
        ami, ari = float("nan"), float("nan")
    else:
        louvain_res = benchmark_fun.find_louvain_res(adata, nclust)
        adata_pp = adata.copy()
        sc.pp.neighbors(adata_pp)
        sc.tl.louvain(adata_pp, resolution=louvain_res, random_state=SEED)
        # cluster_cosistent returns >2 values; original used [0]/[1] -> ami, ari
        modality_consis = benchmark_fun.cluster_cosistent(adata_pp.obs, cell_test, "louvain")
        ami, ari = modality_consis[0], modality_consis[1]

    pdist_df = benchmark_fun.lat2pdist(lat_df, bc_test, OPTION)
    # samecell_dist = RNA-vs-ATAC distance of the SAME cell -> needs paired cells; empty when unpaired
    samecell_dist = benchmark_fun.get_parallel_dist(pdist_df, cell_test) if len(cell_test) else []
    pdist_dfintera = pdist_df.loc[pdist_df.index.isin(bc_test_rna)]
    pdist_dfintera = pdist_dfintera[[c for c in pdist_dfintera.columns if c in bc_test_atac]]
    pdist_df_rna = pdist_df.loc[pdist_df.index.isin(bc_test_rna)]
    pdist_df_rna = pdist_df_rna[[c for c in pdist_df_rna.columns if c in bc_test_rna]]
    pdist_df_atac = pdist_df.loc[pdist_df.index.isin(bc_test_atac)]
    pdist_df_atac = pdist_df_atac[[c for c in pdist_df_atac.columns if c in bc_test_atac]]

    celltype_rows, intra_celltype_sameomics_mean = [], []
    for type_ in celltype:
        celltype_atac = [k for k, v in label_dic.items() if v == type_ and k.endswith("_atac")]
        celltype_rna = [k for k, v in label_dic.items() if v == type_ and k.endswith("_rna")]
        di = pdist_dfintera[[c for c in pdist_dfintera.columns if c in celltype_atac]]
        dr = pdist_df_rna[[c for c in pdist_df_rna.columns if c in celltype_rna]]
        da = pdist_df_atac[[c for c in pdist_df_atac.columns if c in celltype_atac]]

        # get_celltype_dist returns >=2 values; original indexed [0]/[1] (wiclust/btclust)
        d_inter = benchmark_fun.get_celltype_dist(di, label_dic)
        wi_inter, bt_inter = d_inter[0], d_inter[1]
        ks1 = stats.kstest(bt_inter, wi_inter, alternative="less")
        d_rna = benchmark_fun.get_celltype_dist(dr, label_dic)
        d_atac = benchmark_fun.get_celltype_dist(da, label_dic)
        wi_rna, bt_rna = d_rna[0], d_rna[1]
        wi_atac, bt_atac = d_atac[0], d_atac[1]
        ks2 = stats.kstest(bt_atac + bt_rna, wi_atac + wi_rna, alternative="less")

        intra_celltype_sameomics_mean.extend(
            dr.loc[dr.index.isin(celltype_rna)].mean(axis=1).tolist())
        intra_celltype_sameomics_mean.extend(
            da.loc[da.index.isin(celltype_atac)].mean(axis=1).tolist())
        celltype_rows.append([method, type_, ks2.statistic, ks1.statistic])

    celltype_df = pd.merge(pd.DataFrame(celltype_rows), omics_asw_celltype,
                           left_on=1, right_index=True)
    celltype_df.columns = ["method", "celltype", "ks_celltype", "ks_inter_celltype",
                           "celltype_omics_ASW"]
    # ks.statistic compares same-omics-celltype distances vs SAME-CELL distances -> needs paired cells
    ks_stat = (stats.kstest(intra_celltype_sameomics_mean, samecell_dist,
                            alternative="less").statistic if len(cell_test) else float("nan"))
    sum_row = dict(method=method, asw=asw, asw_random=asw_random, omics_asw=omics_asw_ave,
                   ami=ami, ari=ari)
    sum_row["ks.statistic"] = ks_stat
    return sum_row, celltype_df


def knn_accuracy(lat_df, method, ctx, k=10):
    """KNN (cosine, distance-weighted) ATAC-label transfer from RNA -> overall + per-class
    accuracy + the predicted labels. Same as the old knn_ataclabel at k=10."""
    labels_df = ctx["labels_df"]
    rna = [b for b in ctx["bc_test_rna"] if b in lat_df.index]
    atac = [b for b in ctx["bc_test_atac"] if b in lat_df.index]
    y_rna = labels_df.loc[rna, "cell_type"]
    y_atac = labels_df.loc[atac, "cell_type"]
    clf = KNeighborsClassifier(n_neighbors=k, metric="cosine", weights="distance")
    clf.fit(lat_df.loc[rna], y_rna)
    pred = clf.predict(lat_df.loc[atac])
    rows = [[method, accuracy_score(y_atac, pred), "overall"]]
    cm = confusion_matrix(y_atac, pred)
    acc_diag = np.diagonal(cm) / np.sum(cm, axis=1)
    for j, t in enumerate(np.unique(pred)):
        rows.append([method, acc_diag[j], t])
    accu_df = pd.DataFrame(rows, columns=["method", "accuracy", "type"])
    pred_df = pd.DataFrame({method: pred}, index=atac)
    return accu_df, pred_df


def append_dedup(df, path, key_cols):
    """Append df to the CSV at path, dropping any prior rows for the same method (idempotent)."""
    if os.path.exists(path):
        old = pd.read_csv(path)
        merged = pd.concat([old[~old[key_cols[0]].isin(df[key_cols[0]].unique())], df])
    else:
        merged = df
    merged.to_csv(path, index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--latent", required=True, help="method latent.csv (index=barcodes _rna/_atac)")
    ap.add_argument("--method", required=True, help="method display name (matrix row label)")
    ap.add_argument("--label", required=True, help="label.csv (modality, cell_type, random_celltype)")
    ap.add_argument("--out", default=".", help="output dir for the shared metric CSVs")
    ap.add_argument("--knn-label-out", default=None,
                    help="optional path to append the predicted-ATAC-label column (for peak step)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ctx = load_label_context(args.label)
    lat_df = pd.read_csv(args.latent, index_col=0)
    # Score only on barcodes that are in label.csv. A method may embed more cells than the
    # benchmark's labeled test set (e.g. MIRA runs all of TEST_H5 = 1642, but label.csv has the
    # 1531 labeled test cells); those extra cells have no ground-truth cell_type, so they can't be
    # scored and would KeyError in cluster_cosistent. Dropping them up front = every method scored
    # on the identical labeled cells.
    n_before = lat_df.shape[0]
    lat_df = lat_df[lat_df.index.isin(ctx["labels_df"].index)]
    dropped = n_before - lat_df.shape[0]
    print(f"[{args.method}] {lat_df.shape[0]} cells"
          + (f"  (dropped {dropped} not in label.csv)" if dropped else ""))

    sum_row, celltype_df = integration_metrics(lat_df, args.method, ctx)
    accu_df, pred_df = knn_accuracy(lat_df, args.method, ctx, k=10)

    append_dedup(pd.DataFrame([sum_row]), os.path.join(args.out, "sum_metrics.csv"), ["method"])
    append_dedup(celltype_df, os.path.join(args.out, "celltype_metrics.csv"), ["method"])
    append_dedup(accu_df, os.path.join(args.out, "knn_pred_accu.csv"), ["method"])
    pred_path = args.knn_label_out or os.path.join(args.out, f"knn_pred_label__{args.method}.csv")
    pred_df.to_csv(pred_path)
    print(f"[{args.method}] done -> {args.out}")


if __name__ == "__main__":
    main()
