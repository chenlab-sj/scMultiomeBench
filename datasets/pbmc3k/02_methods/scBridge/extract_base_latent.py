#!/usr/bin/env python3
"""
Recover the scBridge base (rep1) latent into the repo as scBridge/latent.csv, so the repo carries
rep1 alongside rep2/rep3 and fig4/reproduce_metrics.py can read all three from {REPO}/scBridge/...
instead of the cluster {LSA}/scBridge/pbmc3k/rep1/latent.csv.

No re-run: the base run already computed the harmony embedding and stored it in .obsm["Embedding"]
of the two -integrated.h5ad files (2_save_scBridge_result.py:35-36,51-52). This just reads it back
and re-does the source(rna)+target(atac) concat exactly as 2_save_scBridge_result.py:56-61.
"""
import os
import scanpy as sc
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
src = sc.read_h5ad(os.path.join(HERE, "scBridge_rna-integrated.h5ad"))    # source = RNA (annotated)
tgt = sc.read_h5ad(os.path.join(HERE, "scBridge_atac-integrated.h5ad"))   # target = ATAC (query)

src_lat = pd.DataFrame(src.obsm["Embedding"], index=src.obs.index)
tgt_lat = pd.DataFrame(tgt.obsm["Embedding"], index=tgt.obs.index)
combined = pd.concat([src_lat, tgt_lat], axis=0)

out = os.path.join(HERE, "latent.csv")
combined.to_csv(out)
print(f"wrote {out}: {combined.shape[0]} cells x {combined.shape[1]} dims "
      f"(rna {src_lat.shape[0]} + atac {tgt_lat.shape[0]})")
print("first barcodes:", list(combined.index[:2]), "...", list(combined.index[-2:]))
