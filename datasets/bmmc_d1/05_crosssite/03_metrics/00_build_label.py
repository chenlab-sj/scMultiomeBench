#!/usr/bin/env python3
# Build label.csv (+ label_major.csv) for the BMMC same-donor CROSS-SITE run:
#   test scRNA  = s1d1 multiome RNA   (index <bc>-1_rna)
#   test scATAC = s4d1 multiome ATAC  (index <bc>-1_atac; label = that s4d1 cell's own RNA label)
# Same donor (d1), different processing sites (1 vs 4) -> isolates the technical/site batch effect with
# NO donor-biology confound. GENUINELY UNPAIRED (s1d1 RNA cells != s4d1 ATAC cells), so it mirrors the
# pbmc_parse/crossdonor pipeline: the same-cell metrics (ks.statistic/ari/ami) are dropped downstream.
#
# Barcodes are normalized to BARE <bc>-1_rna / <bc>-1_atac to match the normalized latents
# (01_prep_latents.py strips the site digit: _rna3->_rna, _atac2->_atac; scVI _expression/_accessibility).
# celltype.csv carries only RNA-suffixed barcodes (one per multiome cell); the ATAC label reuses the
# RNA label of the same cell (rewrite the suffix). No train rows are needed: 00_compute_metrics.py filters
# "train multiomics" out, and only test cells are scored.
import os
import pandas as pd
import numpy as np

HERE      = os.path.dirname(os.path.abspath(__file__))
CELLTYPE  = os.path.normpath(os.path.join(HERE, "..", "..", "scripts", "celltype.csv"))
OUT       = os.path.join(HERE, "label.csv")
OUT_MAJOR = os.path.join(HERE, "label_major.csv")

RNA_BATCH      = "s1d1"     # test scRNA site
ATAC_BATCH     = "s4d1"     # test scATAC site (different site, same donor d1)
MIN_CELLS      = 10         # drop a cell type with < this many RNA cells
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

# ---- ATAC side (s4d1): the same cell's RNA label; barcode <bc>-1_rna2 -> <bc>-1_atac ----
atac = ct[ct["batch"] == ATAC_BATCH].copy()
atac.index = atac["barcode"].str.replace(r"_rna\d+$", "_atac", regex=True)
atac = atac[["cell_type"]]
atac["modality"] = "test_scATAC"

# ---- score only the cell types SHARED by the s1d1 RNA and s4d1 ATAC sides ----
shared = set(rna["cell_type"].dropna().unique()) & set(atac["cell_type"].dropna().unique())
rna_drop  = sorted(set(rna["cell_type"].unique())  - shared)
atac_drop = sorted(set(atac["cell_type"].unique()) - shared)
rna  = rna[rna["cell_type"].isin(shared)]
atac = atac[atac["cell_type"].isin(shared)]

# ---- random baselines (random_atac ~ RNA cell-type distribution; shuffles for celltype/modality) ----
p = rna["cell_type"].value_counts(normalize=True)
atac["random_atac"] = np.random.choice(p.index, size=atac.shape[0], p=p.values)
labels = pd.concat([rna, atac], axis=0)
labels["random_celltype"] = labels["cell_type"].sample(frac=1, random_state=42).values
labels["random_modality"] = labels["modality"].sample(frac=1, random_state=42).values
labels = labels.reindex(columns=["cell_type", "modality", "random_atac",
                                 "random_celltype", "random_modality"])
labels.index.name = None
labels.to_csv(OUT)

# ---- major basis: cell types with > MAJOR_ATAC_MIN TEST-ATAC cells (the composite-score fix) ----
atac_counts = atac["cell_type"].value_counts()
major_types = set(atac_counts[atac_counts > MAJOR_ATAC_MIN].index)
labels_major = labels[labels["cell_type"].isin(major_types)]
labels_major.to_csv(OUT_MAJOR)

print(f"wrote {OUT}")
print(f"  test scRNA (s1d1) : {rna.shape[0]:>6}")
print(f"  test scATAC (s4d1): {atac.shape[0]:>6}")
print(f"  SCORED shared types ({len(shared)}): {sorted(shared)}")
print(f"  dropped (RNA-only) : {rna_drop}")
print(f"  dropped (ATAC-only): {atac_drop}")
print(f"wrote {OUT_MAJOR}: MAJOR types (ATAC>{MAJOR_ATAC_MIN}) ({len(major_types)}): {sorted(major_types)}")
ctab = pd.DataFrame({"RNA(s1d1)": rna["cell_type"].value_counts(),
                     "ATAC(s4d1)": atac["cell_type"].value_counts()}).fillna(0).astype(int)
print(ctab.sort_values("ATAC(s4d1)", ascending=False).to_string())
