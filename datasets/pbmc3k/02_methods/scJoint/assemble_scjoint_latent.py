#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""Assemble scJoint_latent.csv from the raw embeddings scJoint writes to <SCJOINT_ROOT>/output/.
Reproduces the rep1-3 latent format: rna embeddings indexed <barcode>_rna stacked on atac embeddings
indexed <barcode>_atac (barcodes from scJoint_need/{rna,atac}_bc.txt, order matches the .npz inputs).
Run after main.py. Env: SCJOINT_ROOT (default /path/to/tools/scJoint)."""
import os
import numpy as np
import pandas as pd

ROOT = os.environ.get("SCJOINT_ROOT", "/path/to/tools/scJoint")
OUT  = os.path.join(ROOT, "output")
NEED = os.path.join(ROOT, "pbmc3k", "scJoint_need")

rna  = np.loadtxt(os.path.join(OUT, "pbmc3k_rna_scjoint_embeddings.txt"))
atac = np.loadtxt(os.path.join(OUT, "pbmc3k_atac_gene_scjoint_embeddings.txt"))
rna_bc  = [l.strip() for l in open(os.path.join(NEED, "rna_bc.txt"))  if l.strip()]
atac_bc = [l.strip() for l in open(os.path.join(NEED, "atac_bc.txt")) if l.strip()]
assert len(rna_bc)  == rna.shape[0],  f"rna bc {len(rna_bc)} != emb {rna.shape[0]}"
assert len(atac_bc) == atac.shape[0], f"atac bc {len(atac_bc)} != emb {atac.shape[0]}"

rna_df  = pd.DataFrame(rna,  index=[b + "_rna"  for b in rna_bc])
atac_df = pd.DataFrame(atac, index=[b + "_atac" for b in atac_bc])
lat = pd.concat([rna_df, atac_df])
lat.to_csv(os.path.join(OUT, "scJoint_latent.csv"))
print(f"wrote {OUT}/scJoint_latent.csv  shape={lat.shape}")
