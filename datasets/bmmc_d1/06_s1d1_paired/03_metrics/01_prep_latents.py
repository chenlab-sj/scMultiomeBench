#!/usr/bin/env python3
# Normalize every BMMC method latent to BARE barcodes  <bc>-1_rna / <bc>-1_atac  -> ./latents/<key>.csv,
# so 00_compute_metrics.py / 01_adjust_accuracy.py / 02_peak_similarity.R / plot_metrics_matrix.R all run VERBATIM
# (they key on the bare _rna/_atac suffix). The raw latents come in THREE barcode families (verified by
# inspecting the actual latents -- do NOT assume):
#   * "strip": <bc>-1_rna3           (keeps 10x "-1")            -> drop the site digit:  _rna3 -> _rna
#              seurat3,bindsc,portal,scjoint,scglue,scglue_paired,simba,maxfuse,scBridge
#   * "noc1" : <bc>_rna3             (MultiVI/scButterfly dropped the "-1") -> ADD it: _rna3 -> -1_rna
#              midas, scbutterfly
#   * "scvi" : <bc>_rna3_expression / <bc>_atac2_accessibility / <bc>_paired  -> -1_rna / -1_atac / -1
#              scVI (MultiVI)
# The latent dir (.../BMMC_d1/<exp>) is derived from THIS file's path, so the script is identical for
# both the crosssite and s1d1_paired experiments.
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))                          # .../BMMC_d1/benchmark/<exp>
EXP  = os.path.basename(HERE)                                              # crosssite | s1d1_paired
XD   = os.path.normpath(os.path.join(HERE, "..", "..", EXP))               # .../BMMC_d1/<exp> (latents)
OUT  = os.path.join(HERE, "latents")
os.makedirs(OUT, exist_ok=True)

# key -> (latent path relative to XD, barcode family)
METHODS = {
    "seurat3":       ("seurat3/coembed_coor.csv",     "strip"),
    "bindsc":        ("bindsc/coembed_coor.csv",      "strip"),
    "portal":        ("portal/lat_df.csv",            "strip"),
    "scjoint":       ("scjoint/scJoint_latent.csv",   "strip"),
    "scglue":        ("scglue/latent.csv",            "strip"),
    "scglue_paired": ("scglue_paired/latent.csv",     "strip"),
    "simba":         ("simba/latent.csv",             "strip"),
    "maxfuse":       ("maxfuse/latent.csv",           "strip"),
    "scBridge":      ("scBridge/latent.csv",          "strip"),
    "midas":         ("midas/latent.csv",             "noc1"),
    "scbutterfly":   ("scbutterfly/latent.csv",       "noc1"),
    "scVI":          ("scVI/latent.csv",              "scvi"),
}


def normalize(idx, rule):
    if rule == "scvi":
        return (idx.str.replace(r"_rna\d+_expression$", "-1_rna", regex=True)
                   .str.replace(r"_atac\d+_accessibility$", "-1_atac", regex=True)
                   .str.replace(r"_paired$", "-1", regex=True))
    if rule == "noc1":                                     # no "-1" present -> add it
        return idx.str.replace(r"_(rna|atac)\d+$", r"-1_\1", regex=True)
    return idx.str.replace(r"_(rna|atac)\d+$", r"_\1", regex=True)   # "strip": keep "-1", drop digit


for key, (rel, rule) in METHODS.items():
    src = os.path.join(XD, rel)
    if not os.path.exists(src):
        print(f"  SKIP (missing): {key}  {src}")
        continue
    lat = pd.read_csv(src, index_col=0)
    lat.index = normalize(lat.index.astype(str), rule)
    dst = os.path.join(OUT, f"{key}.csv")
    lat.to_csv(dst)
    n_rna = int(lat.index.str.endswith("_rna").sum())
    n_atac = int(lat.index.str.endswith("_atac").sum())
    print(f"  {key:<14} {lat.shape[0]:>6} cells  (_rna={n_rna}, _atac={n_atac})  -> latents/{key}.csv")
