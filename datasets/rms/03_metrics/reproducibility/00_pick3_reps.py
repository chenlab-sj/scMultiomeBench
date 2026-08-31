#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""RMS Fig4 helper: evaluate every 3-rep subset for selected methods.

Reads staged latents from RMS/benchmark/staged/<method>/rep*.csv, clusters each
replicate with the same test-cell label target used by RMS Fig4, and writes:
  pick3/<method>_pick3.csv
  pick3/<method>_pairwise_nmi.csv
  pick3/pick3_summary.csv

Run on the cluster benchmark_env because this needs scanpy + benchmark_fun:
  python 00_pick3_reps.py
  METHODS=scJoint,scVI,Portal,BindSC,Cobolt,scDART python 00_pick3_reps.py
"""

import glob
import itertools
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import confusion_matrix, normalized_mutual_info_score
from sklearn.neighbors import KNeighborsClassifier

sys.path.append(os.environ.get("BENCHMARK_FUN_DIR", "/path/to/multiomeBench/common"))
import benchmark_fun


HERE = Path(__file__).resolve().parent
STAGE = (HERE / ".." / "staged").resolve()
LABEL = Path(os.environ.get("LABEL", HERE / "label.csv"))
OUT = Path(os.environ.get("OUT", HERE / "pick3"))
METHODS = [
    m.strip()
    for m in os.environ.get("METHODS", "scJoint,scVI,Portal,BindSC,Cobolt,scDART").split(",")
    if m.strip()
]
K = int(os.environ.get("K", "10"))
SEED = int(os.environ.get("SEED", "420"))
MIN_CELLS = int(os.environ.get("MIN_CELLS", "100"))
FIX_RES = os.environ.get("FIX_RES", "0") == "1"


def rep_number(path):
    match = re.search(r"rep(\d+)\.csv$", str(path))
    if not match:
        raise ValueError(f"cannot parse replicate number from {path}")
    return int(match.group(1))


def load_labels():
    labels = pd.read_csv(LABEL, index_col=0).dropna(subset=["cell_type"])
    test_mask = labels["modality"].astype(str).str.startswith("test")
    test_all = labels.index[test_mask]
    test_rna = [b for b in test_all if str(b).endswith("_rna")]
    test_atac = [b for b in test_all if str(b).endswith("_atac")]
    major_counts = labels.loc[test_atac, "cell_type"].value_counts()
    major_types = major_counts[major_counts > MIN_CELLS].index.tolist()
    test_major = labels.index[test_mask & labels["cell_type"].isin(major_types)]
    return labels, test_rna, test_atac, test_major, major_types


def read_latent(path):
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip('"')
    return df


def cluster_reps(latents, test_major, nclust):
    clusters = {}
    fixed_res = None
    if FIX_RES:
        first = next(iter(latents.values()))
        lt = first.loc[first.index.isin(test_major)]
        adata = sc.AnnData(lt)
        adata.obsm["X_emb"] = lt.values
        fixed_res = benchmark_fun.find_louvain_res(adata, nclust)

    for rep, df in latents.items():
        lt = df.loc[df.index.isin(test_major)]
        adata = sc.AnnData(lt)
        adata.obsm["X_emb"] = lt.values
        res = fixed_res if FIX_RES else benchmark_fun.find_louvain_res(adata, nclust)
        app = adata.copy()
        sc.pp.neighbors(app)
        sc.tl.louvain(app, resolution=res, random_state=SEED)
        clusters[rep] = pd.Series(app.obs["louvain"].astype(str).values, index=lt.index)
        print(f"  rep{rep}: res={res:.4f} -> {app.obs['louvain'].nunique()} clusters (target nclust={nclust})", flush=True)
    return clusters


def predict_atac_labels(latents, labels, test_rna, test_atac):
    preds = {}
    for rep, df in latents.items():
        rna = [b for b in test_rna if b in df.index]
        atac = [b for b in test_atac if b in df.index]
        if not rna or not atac:
            continue
        clf = KNeighborsClassifier(n_neighbors=K, metric="cosine", weights="distance")
        clf.fit(df.loc[rna], labels.loc[rna, "cell_type"])
        preds[rep] = pd.Series(clf.predict(df.loc[atac]), index=atac)
    return preds


def subset_nmi(clusters, combo):
    cc = pd.concat([clusters[r] for r in combo], axis=1, join="inner")
    cc.columns = list(combo)
    return float(
        np.mean(
            [
                normalized_mutual_info_score(cc[a], cc[b])
                for a, b in itertools.combinations(combo, 2)
            ]
        )
    )


def subset_repro(preds, labels, combo):
    if any(rep not in preds for rep in combo):
        return np.nan
    knn = pd.concat([preds[rep] for rep in combo], axis=1, join="inner")
    knn.columns = list(combo)
    actual = labels.loc[knn.index, "cell_type"].values
    counts = pd.Series(actual).value_counts()
    categories = counts.index.tolist()
    major = [c for c in categories if counts[c] > MIN_CELLS]
    cms = []
    for rep in combo:
        cm = confusion_matrix(actual, knn[rep].values, labels=categories).astype(float)
        cms.append(cm / cm.sum(axis=1, keepdims=True))
    cm_std = np.std(np.stack(cms, 0), axis=0)
    return 1 - float(np.mean([np.mean(cm_std[categories.index(ct), :]) for ct in major]))


def write_method(method, paths, labels, test_rna, test_atac, test_major, major_types):
    print(f"\n{method}: {len(paths)} staged reps", flush=True)
    latents = {rep_number(p): read_latent(p) for p in paths}
    clusters = cluster_reps(latents, test_major, len(major_types))
    preds = predict_atac_labels(latents, labels, test_rna, test_atac)

    pair_rows = []
    for a, b in itertools.combinations(sorted(latents), 2):
        cc = pd.concat([clusters[a], clusters[b]], axis=1, join="inner")
        cc.columns = [a, b]
        pair_rows.append({
            "Method": method,
            "Pair": f"{method}-{a} vs {method}-{b}",
            "NMI": float(normalized_mutual_info_score(cc[a], cc[b])),
        })
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(OUT / f"{method}_pairwise_nmi.csv", index=False)

    rows = []
    for combo in itertools.combinations(sorted(latents), 3):
        rows.append({
            "Method": method,
            "subset": "rep " + "+".join(map(str, combo)),
            "Fig4A_repro": round(subset_repro(preds, labels, combo), 4),
            "Fig4B_NMI": round(subset_nmi(clusters, combo), 4),
        })
    out = pd.DataFrame(rows).sort_values("Fig4B_NMI", ascending=True)
    out.to_csv(OUT / f"{method}_pick3.csv", index=False)
    print(out.to_string(index=False), flush=True)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    labels, test_rna, test_atac, test_major, major_types = load_labels()
    print(f"label={LABEL}")
    print(f"stage={STAGE}")
    print(f"out={OUT}")
    print(f"FIX_RES={int(FIX_RES)} | major cell types (> {MIN_CELLS} ATAC cells): {major_types}")

    summaries = []
    for method in METHODS:
        paths = sorted(
            glob.glob(str(STAGE / method / "rep*.csv")),
            key=rep_number,
        )
        if len(paths) < 3:
            print(f"\n{method}: only {len(paths)} staged reps -> skipped")
            continue
        summaries.append(
            write_method(method, paths, labels, test_rna, test_atac, test_major, major_types)
        )

    if summaries:
        summary = pd.concat(summaries, ignore_index=True)
        summary.to_csv(OUT / "pick3_summary.csv", index=False)
        best_low = summary.sort_values(["Method", "Fig4B_NMI"]).groupby("Method").head(1)
        print("\nLowest-NMI subset per method:")
        print(best_low.to_string(index=False))
        print(f"\nwrote {OUT / 'pick3_summary.csv'}")


if __name__ == "__main__":
    main()
