#!/usr/bin/env python3
# scJoint assembler (pbmc3k ArchR gene-activity variant): build scJoint_latent.csv from main.py's
# embeddings + the stage-0 barcodes.
# - main.py writes embeddings to the SHARED package output dir, OVERWRITTEN by every scJoint run; the
#   main sub snapshots that into ./output before this runs, so we read the snapshot (race-safe).
# - Embedding rows are position-aligned to rna_bc.txt / atac_bc.txt (scJoint preserves input order).
# - pbmc3k is a paired multiome SPLIT, so we append _rna / _atac to the (plain) barcodes to match
#   metrics/label.csv -- compute_metrics keys cells off the <stem>_rna / <stem>_atac suffix.
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "output")          # snapshot of the package output (see pbmc3k_scjoint_main_sub.sh)
NEED = os.path.join(HERE, "scJoint_need")

# Read the EXACT pbmc3k_* embedding names (deterministic from the config npz basenames). The shared
# package output also contains STALE non-prefixed embeddings from other scJoint runs (e.g. cross-donor
# atac = 5783 cells) -- a loose glob would grab those and mismatch the 1642 pbmc3k barcodes.
rna_e  = pd.read_csv(os.path.join(OUT, "pbmc3k_rna_scjoint_embeddings.txt"),       sep=r"\s+", header=None)
atac_e = pd.read_csv(os.path.join(OUT, "pbmc3k_atac_gene_scjoint_embeddings.txt"), sep=r"\s+", header=None)
rna_bc  = pd.read_csv(os.path.join(NEED, "rna_bc.txt"),  header=None)[0].astype(str) + "_rna"
atac_bc = pd.read_csv(os.path.join(NEED, "atac_bc.txt"), header=None)[0].astype(str) + "_atac"
assert len(rna_e)  == len(rna_bc),  (len(rna_e),  len(rna_bc))
assert len(atac_e) == len(atac_bc), (len(atac_e), len(atac_bc))

rna_e.index = rna_bc.values
atac_e.index = atac_bc.values
lat = pd.concat([rna_e, atac_e]); lat.columns = range(lat.shape[1])
lat.to_csv(os.path.join(OUT, "scJoint_latent.csv"))
print(f"wrote {OUT}/scJoint_latent.csv  shape={lat.shape}  ({len(rna_e)} rna + {len(atac_e)} atac)")
