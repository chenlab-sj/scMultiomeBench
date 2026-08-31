#!/usr/bin/env python3
"""
Fix the RMS scVI latent barcodes WITHOUT re-running scVI.

01_run_scvi.py used set_axis(split("_",1)[0]) on adata_mvi.obs, which STRIPPED the
organize_multiome_anndatas modality tag and collapsed test RNA + ATAC to the same bare barcode
(<bc>) -- so the saved latent can't be joined to label.csv's <bc>-1_rna / <bc>-1_atac.

organize_multiome_anndatas orders rows deterministically as [train, test_rna, test_atac], and the
two test blocks carry identical barcodes in the same order, so we recover modality by ROW POSITION
and append the label suffix. Verified: 8055/8055 label rna AND atac cells match.

(BRCA did this differently but equivalently: its 01_run_scvi.py has the strip COMMENTED OUT, keeping the
<bc>_expression / <bc>_accessibility tag, then maps _expression->_rna, _accessibility->_atac.)

Writes scvi_latent.csv (test cells only, <bc>-1_rna / <bc>-1_atac) next to each latent.csv.
"""
import os
import pandas as pd

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "res_Mast607A")
TEST = 8291   # test-h5 cells per modality; train = n - 2*TEST (32115 = 15533 + 2*8291)

# rep1 (base) was already prepared as scVI/scvi_latent.csv by the original benchmark (train "_Mast39"
# + test "_rna"/"_atac"); the reproduce reps live under res_Mast607A/ and need this fix. rep4/rep5 are
# the seed-fragility spot-check (seeds 13/99); a missing rep is skipped, so it's safe to run anytime.
for rep in ["rep2", "rep3", "rep4", "rep5"]:
    src = os.path.join(BASE, rep, "latent.csv")
    if not os.path.exists(src):
        print(f"skip (missing): {src}")
        continue
    df = pd.read_csv(src, index_col=0)
    n = len(df)
    train = n - 2 * TEST
    rna_blk = df.index[train:train + TEST]
    atac_blk = df.index[train + TEST:n]
    assert (rna_blk.values == atac_blk.values).all(), \
        f"{src}: test blocks misaligned -- not [train, test_rna, test_atac] order"
    test = df.iloc[train:].copy()
    test.index = [f"{b}-1_rna" for b in rna_blk] + [f"{b}-1_atac" for b in atac_blk]
    dst = os.path.join(BASE, rep, "scvi_latent.csv")
    test.to_csv(dst)
    print(f"{src} -> {dst}  ({test.shape[0]} cells = {TEST} rna + {TEST} atac)")

print("DONE")
