#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# scJoint stage 4: assemble latent.csv from main.py's embeddings + the s1 barcodes.
# Reads THIS experiment's own copy of the embeddings (s3 does `cp -r output/ $outdir`), not the
# shared /path/to/tools/scJoint/output. That shared dir is reused by every scJoint run in the
# repo and still held pbmc3k_* leftovers, so reading it directly meant a stale-embedding read
# whenever main.py failed or another experiment ran in between (the 8291-vs-8442 assert).
# The RNA/ATAC embedding rows are position-aligned to rna_bc.txt / atac_bc.txt (scJoint preserves
# input order).
import os
import pandas as pd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
XD  = "/path/to/multiomeBench/BMMC_d1/crosssite/scjoint"

rna_e  = pd.read_csv(f"{OUT}/rna_scjoint_embeddings.txt",       sep=r"\s+", header=None)
atac_e = pd.read_csv(f"{OUT}/atac_gene_scjoint_embeddings.txt", sep=r"\s+", header=None)
rna_bc  = pd.read_csv(f"{XD}/rna_bc.txt",  header=None)[0].astype(str)
atac_bc = pd.read_csv(f"{XD}/atac_bc.txt", header=None)[0].astype(str)
assert len(rna_e) == len(rna_bc),  (len(rna_e), len(rna_bc))
assert len(atac_e) == len(atac_bc), (len(atac_e), len(atac_bc))

rna_e.index = rna_bc.values
atac_e.index = atac_bc.values
lat = pd.concat([rna_e, atac_e]); lat.columns = range(lat.shape[1])
lat.to_csv(f"{XD}/scJoint_latent.csv")
print(f"wrote {XD}/scJoint_latent.csv  shape={lat.shape}")
