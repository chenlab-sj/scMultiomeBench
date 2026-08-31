#!/usr/bin/env python3
"""
BMMC_d1 stage 2: score the staged latents with the ORIGINAL batch-aware benchmark_metrics() (3 test
batches s2d1/s4d1/s1d1; computes sample_asw + ks_sample for the batch-correction term, via
get_parallel_dist3 / cluster_cosistent3). The metric function is imported UNCHANGED from benchmark_metrics_lib.py;
only the input wiring (staged latents + the current 13 base + 6 (batch) method list) is new.

Output (consumed by the plot): sum_metrics.csv (incl. sample_asw) + celltype_metrics.csv (incl. ks_sample).
Run in benchmark_env (scib + benchmark_fun) on the cluster.  python 02_run_bmmc_metrics.py
"""
import os
import re
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_bmmc import benchmark_metrics          # original batch-aware fn (pulls in scib + benchmark_fun)

HERE  = os.path.dirname(os.path.abspath(__file__))
STAGE = os.path.join(HERE, "staged")
LABEL = os.path.join(HERE, "old", "BMMC_d1", "label.csv")
OUT   = HERE + "/"                                   # -> sum_metrics.csv + celltype_metrics.csv here

PIPELINES = ["BindSC", "Seurat(CCA)", "MaxFuse", "MIDAS", "scBridge", "scVI", "scButterfly", "scDART",
             "scglue", "scglue(multiome)", "scJoint", "simba", "Portal",
             "BindSC(batch)", "scVI(batch)", "scglue(batch)", "scglue(multiome,batch)",
             "scJoint(batch)", "MIDAS(batch)"]      # Portal re-added for Fig6; Cobolt to be added when its latent lands

# Portal-only mode (PORTAL_ONLY=1): score just Portal into a portal_only/ subdir (benchmark_metrics
# OVERWRITES sum_metrics.csv/celltype_metrics.csv, so we must NOT point it at HERE), then APPEND Portal's
# rows to the existing HERE/sum_metrics.csv + celltype_metrics.csv -- the 18 done methods aren't re-run.
PORTAL_ONLY = os.environ.get("PORTAL_ONLY") == "1"
if PORTAL_ONLY:
    PIPELINES = ["Portal"]
    OUT = os.path.join(HERE, "portal_only") + "/"
    os.makedirs(OUT, exist_ok=True)

# ---- label-derived inputs (mirror the original main() setup) ----
labels_annot = pd.read_csv(LABEL, index_col=0)
bc_test      = labels_annot[labels_annot["set"] != "train multiomics"].index.values
labels_test  = labels_annot.loc[bc_test]
bc_test_rna  = [b for b in bc_test if re.search(r"_rna.*$", b)]
bc_test_atac = [b for b in bc_test if re.search(r"_atac.*$", b)]
label_dic    = dict(zip(labels_test.index, labels_test["cell_type"]))
celltype     = labels_test["cell_type"].unique().tolist()

# ---- staged latents -> lat_dfs / pipelines ----
lat_dfs, pipelines = [], []
for m in PIPELINES:
    f = os.path.join(STAGE, m, "latent.csv")
    if not os.path.exists(f):
        print("SKIP (missing staged):", m); continue
    lat_dfs.append(pd.read_csv(f, index_col=0)); pipelines.append(m)
print(f"scoring {len(pipelines)} methods: {pipelines}", flush=True)

# ---- celltest{1,2,3}: RNA-cell stems per test batch (from bc_test_rna), exactly as the original main() ----
celltest1, celltest2, celltest3 = [], [], []
for item in bc_test_rna:
    mo = re.match(r"^(.*)_rna(\d+)$", item)
    if mo:
        pre, suf = mo.group(1), int(mo.group(2))
        (celltest1 if suf == 1 else celltest2 if suf == 2 else celltest3).append(pre)

benchmark_metrics(lat_dfs, pipelines, labels_annot, celltype, bc_test, bc_test_rna, bc_test_atac,
                  celltest1, celltest2, celltest3, label_dic, OUT)

if PORTAL_ONLY:
    # append Portal's rows to the main CSVs (drop any prior Portal row first -> idempotent re-runs)
    for fn, mcol in [("sum_metrics.csv", "method"), ("celltype_metrics.csv", "method")]:
        main_f = os.path.join(HERE, fn)
        new = pd.read_csv(os.path.join(OUT, fn))
        if os.path.exists(main_f):
            base = pd.read_csv(main_f)
            base = base[base[mcol] != "Portal"]
            pd.concat([base, new], ignore_index=True).to_csv(main_f, index=False)
        else:
            new.to_csv(main_f, index=False)
        print(f"appended Portal -> {main_f} (+{len(new)} rows)", flush=True)
print("DONE: sum_metrics.csv + celltype_metrics.csv", flush=True)
