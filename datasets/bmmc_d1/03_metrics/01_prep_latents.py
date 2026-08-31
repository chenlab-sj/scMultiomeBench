#!/usr/bin/env python3
"""
BMMC_d1 stage 1: normalize each method's latent barcodes to label.csv's convention (<bc>-1_rna{1,2,3} /
<bc>-1_atac{1,2,3} for the 3 test batches s2d1/s4d1/s1d1; train cells carry _s* suffixes and aren't scored)
and stage to staged/<method>/latent.csv. Most methods already match; only scVI/scVI2 need remapping
(<bc>_rna1_expression -> <bc>-1_rna1). Prints, per method, how many cells land on the label test set.

Row list = 13 base + 6 (batch) variants (Cobolt to be added later; conos/MinNet2 have no latent).
  python 01_prep_latents.py
"""
import os
import re
import pandas as pd

HERE  = os.path.dirname(os.path.abspath(__file__))
SCR   = os.path.normpath(os.path.join(HERE, "..", "scripts"))
STAGE = os.path.join(HERE, "staged")
LABEL = os.path.join(HERE, "old", "BMMC_d1", "label.csv")

# display name -> (kind, path under scripts/)
METHODS = {
    "BindSC":                 ("standard", "bindsc/coembed_coor.csv"),
    "Seurat(CCA)":            ("standard", "seurat3/coembed_coor.csv"),
    "LIGER":                  ("standard", "liger/coembed_coor.csv"),
    "MaxFuse":                ("standard", "maxfuse/latent.csv"),
    "MIDAS":                  ("standard", "midas/latent.csv"),
    "scBridge":               ("standard", "scBridge/latent.csv"),
    "scVI":                   ("scvi",     "scVI/latent.csv"),
    "scButterfly":            ("standard", "scbutterfly/latent.csv"),
    "scDART":                 ("standard", "scdart/latent.csv"),
    "scglue":                 ("standard", "scglue/latent.csv"),
    "scglue(multiome)":       ("standard", "scglue_paired/latent.csv"),
    "scJoint":                ("standard", "scjoint/scJoint_latent.csv"),
    "simba":                  ("standard", "simba/latent.csv"),
    "Portal":                 ("standard", "portal/lat_df.csv"),   # re-added for Fig6 (published method)
    # ---- (batch) variants ----
    "BindSC(batch)":          ("standard", "bindsc2/coembed_coor.csv"),
    "scVI(batch)":            ("scvi",     "scVI2/latent.csv"),
    "scglue(batch)":          ("standard", "scglue2/latent.csv"),
    "scglue(multiome,batch)": ("standard", "scglue_paired2/latent.csv"),
    "scJoint(batch)":         ("standard", "scjoint2/scJoint_latent.csv"),
    "MIDAS(batch)":           ("standard", "midas2/latent.csv"),
}

_BATCH = re.compile(r"_(rna|atac)[123]$")

def normalize_bc(bc):
    # -> label convention <bc>-1_rna1 : strip scVI's expr/access tags, then insert -1 if it's missing
    # (idempotent for methods already in <bc>-1_rna1 form; fixes scVI <bc>_rna1_expression and MIDAS
    # <bc>_rna1, both of which drop the -1). Train cells (no _rna/_atac[123] tag) are left unchanged.
    bc = bc.replace("_expression", "").replace("_accessibility", "")
    m = _BATCH.search(bc)
    if m and not bc[:m.start()].endswith("-1"):
        return bc[:m.start()] + "-1" + bc[m.start():]
    return bc

def main():
    lab = set(pd.read_csv(LABEL, index_col=0).index.astype(str))
    test = {b for b in lab if _BATCH.search(b)}                         # labelled test cells (<bc>-1_rna1 ...)

    for name, (kind, rel) in METHODS.items():
        path = os.path.join(SCR, rel)
        if not os.path.exists(path):
            print(f"  MISSING {name}: {rel}"); continue
        df = pd.read_csv(path, index_col=0)
        df.index = [normalize_bc(b) for b in df.index.astype(str).str.strip('"')]
        df.columns = [f"latent_{i}" for i in range(df.shape[1])]
        outdir = os.path.join(STAGE, name); os.makedirs(outdir, exist_ok=True)
        df.to_csv(os.path.join(outdir, "latent.csv"))
        onlbl = sum(b in test for b in df.index)
        flag = "" if onlbl > 0.9 * len(test) else "  <-- LOW overlap, check barcodes"
        print(f"  {name:24s} {df.shape[0]:6d} cells -> {onlbl}/{len(test)} on labelled test{flag}")

if __name__ == "__main__":
    main()
