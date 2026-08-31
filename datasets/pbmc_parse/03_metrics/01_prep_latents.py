#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Pre-normalize the two latents whose barcodes don't already match label.csv (<bc>_rna / <bc>_atac).
# 00_compute_metrics.py does NOT normalize barcodes (it drops anything not in label.csv), so scVI
# (_paired/_expression/_accessibility) and cobolt (train~/testrna~/testatac~ prefixes) need fixing
# first -- exactly the pbmc3k convention of pointing the methods array at scvi_latent.csv /
# cobolt_latent.csv. The other 12 methods are already <bc>_rna/_atac and need no prep.
import pandas as pd

ROOT = "/path/to/multiomeBench/pbmc_parse"

scvi = pd.read_csv(f"{ROOT}/scVI/latent.csv", index_col=0)
scvi.index = (scvi.index.str.replace("_paired", "-1")
                        .str.replace("_expression", "_rna")
                        .str.replace("_accessibility", "-1_atac"))
scvi.to_csv(f"{ROOT}/scVI/scvi_latent.csv")

cobolt = pd.read_csv(f"{ROOT}/cobolt/latent.csv", index_col=0)
cobolt.index = cobolt.index.str.replace(r"^.*~", "", regex=True)   # strip train~/testrna~/testatac~
cobolt.to_csv(f"{ROOT}/cobolt/cobolt_latent.csv")

print(f"scvi_latent.csv {scvi.shape}  cobolt_latent.csv {cobolt.shape}")
