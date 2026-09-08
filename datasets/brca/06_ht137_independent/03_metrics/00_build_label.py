#!/usr/bin/env python3
# Build label.csv for the BRCA HT137B1-S1H7 cross-modality benchmark (Fig S4B).
#
# HT137 is same-donor INDEPENDENT sequencing: test scRNA and test scATAC are DIFFERENT cells
# (only ~5 barcode strings coincide by chance) -> genuinely UNPAIRED, like the Parse experiment.
# The RNA and ATAC annotations use DIFFERENT vocabularies, so the ATAC labels are COLLAPSED to the
# RNA vocabulary and the benchmark scores only the cell types shared by both modalities (the user's
# "collapse to RNA vocab" choice). ATAC-only types with no RNA counterpart (e.g. Luminal mature) are
# dropped -- a method can't place a cell type it never saw in the RNA reference, so scoring it is moot.
#
#   test scRNA  (rna_celltype.csv)  -> cell_type = RNA labels,  index <bc>_rna,  modality "test_scRNA"
#   test scATAC (atac_celltype.csv) -> cell_type = ATAC labels HARMONIZED to RNA vocab, index <bc>_atac
#   + random baselines (random_atac / random_celltype / random_modality) for the ADJUSTED metrics
#   + the HT235 train-multiome barcodes as "train multiomics" rows (excluded from scoring; present so
#     category-3 methods' train cells are recognised-and-dropped rather than unmatched).
import os
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))                 # .../BRCA/benchmark/figS4b
HT   = os.path.normpath(os.path.join(HERE, "..", "..", "HT137B1-S1H7"))
RNA_CT   = os.path.join(HT, "rna_celltype.csv")                   # bc(<bc>_rna), cell_type, modality
ATAC_CT  = os.path.join(HT, "atac_celltype.csv")                  # bc(<bc>_atac), cell_type, modality
TRAIN_BC = os.path.join(HT, "train_forHT137B1-S1H7", "barcodes.tsv")
OUT      = os.path.join(HERE, "label.csv")
MIN_CELLS = 10            # keep RNA types with >=10 cells (keeps DC=13; the shared set is ~9 types)

# ATAC finer vocab -> RNA coarse vocab (identity for Tumor/Basal progenitor/Endothelial/Plasma/DC).
ATAC_TO_RNA = {
    "Lum progenitor": "Luminal progenitor",
    "Lum mature":     "Luminal mature",        # no RNA counterpart -> dropped at the intersect below
    "mCAF":           "Fibroblasts",
    "vCAF":           "Fibroblasts",
    "T_NK":           "T-cells",
    "Macrophage":     "Macrophages",
}

np.random.seed(42)

# ---- test scRNA (RNA vocab; barcodes already carry _rna) ----
rna = pd.read_csv(RNA_CT).set_index("bc")[["cell_type", "modality"]]
rna = rna.groupby("cell_type").filter(lambda x: len(x) >= MIN_CELLS)

# ---- test scATAC (harmonise ATAC vocab -> RNA) ----
atac = pd.read_csv(ATAC_CT).set_index("bc")[["cell_type", "modality"]]
atac["cell_type"] = atac["cell_type"].replace(ATAC_TO_RNA)

# ---- keep ONLY cell types present in BOTH modalities (the scored set) ----
# A type in just one modality (RNA-only B-cells/Mast, or ATAC-only Luminal mature) cannot be scored
# cross-modally AND breaks the inter-omics distance: compute_metrics' celltype loop would hand
# get_celltype_dist a 0-column matrix for the missing side -> results.loc[0] KeyError. So restrict
# SYMMETRICALLY to the intersection (mirrors the crossdonor build_label).
shared    = set(rna["cell_type"].dropna().unique()) & set(atac["cell_type"].dropna().unique())
rna_drop  = sorted(set(rna["cell_type"].unique())  - shared)
atac_drop = sorted(set(atac["cell_type"].unique()) - shared)
rna  = rna[rna["cell_type"].isin(shared)]
atac = atac[atac["cell_type"].isin(shared)]

# ---- random baselines: random_atac ~ RNA cell-type distribution; the others are shuffles ----
p = rna["cell_type"].value_counts(normalize=True)
atac["random_atac"] = np.random.choice(p.index, size=atac.shape[0], p=p.values)
labels_test = pd.concat([rna, atac], axis=0)
labels_test["random_celltype"] = labels_test["cell_type"].sample(frac=1, random_state=42).values
labels_test["random_modality"] = labels_test["modality"].sample(frac=1, random_state=42).values

# ---- train-multiome rows (HT235): barcodes only; cell_type left blank (train is excluded from scoring) ----
train_bc = pd.read_csv(TRAIN_BC, header=None)[0].astype(str)
train = pd.DataFrame(index=pd.Index(train_bc, name="bc"),
                     data={"cell_type": np.nan, "modality": "train multiomics"})

cols = ["cell_type", "modality", "random_atac", "random_celltype", "random_modality"]
labels = pd.concat([train, labels_test], axis=0).reindex(columns=cols)
labels.index.name = None
labels.to_csv(OUT)

print(f"wrote {OUT}")
print(f"  test scRNA  : {rna.shape[0]:>6}  ({rna['cell_type'].nunique()} types)")
print(f"  test scATAC : {atac.shape[0]:>6}  ({atac['cell_type'].nunique()} types, harmonised+intersected)")
print(f"  train rows  : {train.shape[0]:>6}")
print(f"  scored cell types (shared, {len(shared)}): {sorted(shared)}")
print(f"  dropped RNA-only types : {rna_drop}")
print(f"  dropped ATAC-only types: {atac_drop}")
print("  per-type counts (RNA / ATAC):")
ct = pd.DataFrame({"RNA": rna['cell_type'].value_counts(),
                   "ATAC": atac['cell_type'].value_counts()}).fillna(0).astype(int)
print(ct.to_string())
