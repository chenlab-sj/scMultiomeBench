#!/usr/bin/env python3
# Normalize the scVI (MultiVI) latents whose barcodes don't match label.csv.
# cobolt/cobolt_latent.csv and scglue_paired/latent.csv are ALREADY label-ready (their test cells
# carry _rna/_atac; only train cells are raw) -> left untouched. Only scVI needs normalization:
#   MultiVI barcodes are <bc>_paired (train) / <bc>_expression (test RNA) / <bc>_accessibility (test
#   ATAC), with the 10x "-1" dropped. Re-add it so they match label.csv (<bc>-1 / <bc>-1_rna / <bc>-1_atac):
#       _paired -> -1 ,  _expression -> -1_rna ,  _accessibility -> -1_atac
# Runs for the main seed-420 run + any seed reps (rep2 seed0 / rep3 seed40) that have finished.
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HT   = os.path.normpath(os.path.join(HERE, "..", "..", "HT137B1-S1H7"))


def norm_scvi(d):
    src = os.path.join(d, "latent.csv")
    if not os.path.exists(src):
        print(f"  skip (no latent.csv yet): {os.path.relpath(d, HT)}")
        return
    lat = pd.read_csv(src, index_col=0)
    lat.index = (lat.index.str.replace("_paired", "-1")
                          .str.replace("_expression", "-1_rna")
                          .str.replace("_accessibility", "-1_atac"))
    out = os.path.join(d, "scvi_latent.csv")
    lat.to_csv(out)
    print(f"  wrote {os.path.relpath(out, HT)}  {lat.shape}")


for sub in ["scVI", "scVI/rep2", "scVI/rep3"]:
    norm_scvi(os.path.join(HT, sub))
