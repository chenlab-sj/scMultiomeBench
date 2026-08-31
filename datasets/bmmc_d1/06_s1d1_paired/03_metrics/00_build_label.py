#!/usr/bin/env python3
# Build label.csv (+ label_major.csv) for the BMMC s1d1 PAIRED-MULTIOME baseline:
#   test scRNA  = s1d1 multiome RNA   (index <bc>-1_rna)
#   test scATAC = s1d1 multiome ATAC  (index <bc>-1_atac; SAME cells as the RNA -> paired)
# This is the multiome (paired) reference for the cross-site scatter: RNA and ATAC come from the SAME
# s1d1 cells, so <bc>-1_rna and <bc>-1_atac share a stem and 00_compute_metrics.py sees them as paired
# (cell_test non-empty). The final score still uses only the 6 metrics shared with the unpaired run
# (plot_metrics_matrix.R with SPLICE_PUBLISHED=0 drops ks.statistic/ari/ami), so the baseline is scored fairly.
#
# Barcodes normalized to BARE <bc>-1_rna / <bc>-1_atac to match the normalized latents (01_prep_latents.py:
# _rna3->_rna, _atac3->_atac; scVI _expression/_accessibility). No train rows needed (compute filters them).
import os
import pandas as pd
import numpy as np

HERE      = os.path.dirname(os.path.abspath(__file__))
CELLTYPE  = os.path.normpath(os.path.join(HERE, "..", "..", "scripts", "celltype.csv"))
OUT       = os.path.join(HERE, "label.csv")
OUT_MAJOR = os.path.join(HERE, "label_major.csv")

RNA_BATCH      = "s1d1"     # test scRNA site
ATAC_BATCH     = "s1d1"     # test scATAC site (SAME site -> paired multiome baseline)
MIN_CELLS      = 10
MAJOR_ATAC_MIN = 100        # "major" cell type = > this many TEST-ATAC cells (rare-celltype composite fix)

np.random.seed(42)

ct = pd.read_csv(CELLTYPE)
ct.columns = ["barcode", "batch", "cell_type"]

# ---- RNA side (s1d1): barcode <bc>-1_rna3 -> <bc>-1_rna ----
rna = ct[ct["batch"] == RNA_BATCH].copy()
rna.index = rna["barcode"].str.replace(r"_rna\d+$", "_rna", regex=True)
rna = rna[["cell_type"]]
rna["modality"] = "test_scRNA"
rna = rna.groupby("cell_type").filter(lambda x: len(x) >= MIN_CELLS)

# ---- ATAC side (s1d1, SAME cells): the cell's own RNA label; barcode <bc>-1_rna3 -> <bc>-1_atac ----
atac = ct[ct["batch"] == ATAC_BATCH].copy()
atac.index = atac["barcode"].str.replace(r"_rna\d+$", "_atac", regex=True)
atac = atac[["cell_type"]]
atac["modality"] = "test_scATAC"

# ---- score the shared cell types (== all s1d1 types kept, since RNA and ATAC are the same cells) ----
shared = set(rna["cell_type"].dropna().unique()) & set(atac["cell_type"].dropna().unique())
rna_drop  = sorted(set(rna["cell_type"].unique())  - shared)
atac_drop = sorted(set(atac["cell_type"].unique()) - shared)
rna  = rna[rna["cell_type"].isin(shared)]
atac = atac[atac["cell_type"].isin(shared)]

# ---- random baselines ----
p = rna["cell_type"].value_counts(normalize=True)
atac["random_atac"] = np.random.choice(p.index, size=atac.shape[0], p=p.values)
labels = pd.concat([rna, atac], axis=0)
labels["random_celltype"] = labels["cell_type"].sample(frac=1, random_state=42).values
labels["random_modality"] = labels["modality"].sample(frac=1, random_state=42).values
labels = labels.reindex(columns=["cell_type", "modality", "random_atac",
                                 "random_celltype", "random_modality"])
labels.index.name = None
labels.to_csv(OUT)

# ---- major basis: cell types with > MAJOR_ATAC_MIN TEST-ATAC cells ----
atac_counts = atac["cell_type"].value_counts()
major_types = set(atac_counts[atac_counts > MAJOR_ATAC_MIN].index)
labels_major = labels[labels["cell_type"].isin(major_types)]
labels_major.to_csv(OUT_MAJOR)

print(f"wrote {OUT}")
print(f"  test scRNA (s1d1) : {rna.shape[0]:>6}")
print(f"  test scATAC (s1d1): {atac.shape[0]:>6}")
print(f"  SCORED shared types ({len(shared)}): {sorted(shared)}")
print(f"  dropped (RNA-only) : {rna_drop}")
print(f"  dropped (ATAC-only): {atac_drop}")
print(f"wrote {OUT_MAJOR}: MAJOR types (ATAC>{MAJOR_ATAC_MIN}) ({len(major_types)}): {sorted(major_types)}")
ctab = pd.DataFrame({"RNA(s1d1)": rna["cell_type"].value_counts(),
                     "ATAC(s1d1)": atac["cell_type"].value_counts()}).fillna(0).astype(int)
print(ctab.sort_values("ATAC(s1d1)", ascending=False).to_string())
