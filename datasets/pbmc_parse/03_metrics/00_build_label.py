#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Build label.csv for the cross-platform benchmark metrics (replicates fileprep.ipynb).
#   test scRNA (Parse)  -> cell_type = Azimuth labels  (rna_celltype), index <bc>_rna
#   test scATAC (pbmc3k) -> cell_type = pbmc3k labels   (atac_celltype), index <bc>_atac
#   + random baselines (random_atac / random_celltype / random_modality) for the ADJUSTED metrics
#   + prepend the pbmc3k train-multiome label.csv
# NB: cell types with < MIN_CELLS cells in the Parse RNA are dropped. MIN_CELLS=30 keeps CD16
#     Monocytes (35 -- a distinct, hard type) while excluding Myeloid DC (27), pDC (8), naive CD4 T (7).
import pandas as pd
import numpy as np

ROOT = "/path/to/multiomeBench"
TEST_ANNOT  = f"{ROOT}/pbmc_parse/input_azimuth/muti_celltype.csv"
TRAIN_ANNOT = "/path/to/multiomeBench/common/pbmc3k/label.csv"
OUT         = f"{ROOT}/pbmc_parse/benchmark/metrics/label.csv"
MIN_CELLS   = 30

np.random.seed(42)
test = pd.read_csv(TEST_ANNOT, index_col=0)

# RNA (Parse, Azimuth labels)
rna = test[test["modality"] == "parse scRNA"].copy()
rna["cell_type"] = rna["rna_celltype"]; rna["modality"] = "test scRNA"
rna = rna[["cell_type", "modality"]]
rna.index = [b + "_rna" for b in rna.index]
rna = rna.groupby("cell_type").filter(lambda x: len(x) >= MIN_CELLS)

# ATAC (pbmc3k labels), restricted to types present in the RNA-kept set
atac = test[test["modality"] == "pbmc3k_test scATAC"].copy()
atac["cell_type"] = atac["atac_celltype"]; atac["modality"] = "test scATAC"
atac = atac[["cell_type", "modality"]]
atac.index = [b + "_atac" for b in atac.index]
atac = atac[atac["cell_type"].isin(set(rna["cell_type"].dropna().unique()))]

# random_atac ~ RNA cell-type distribution; random_celltype/modality = shuffles
p = rna["cell_type"].value_counts(normalize=True)
atac["random_atac"] = np.random.choice(p.index, size=atac.shape[0], p=p.values)
labels_test = pd.concat([rna, atac], axis=0)
labels_test["random_celltype"] = labels_test["cell_type"].sample(frac=1, random_state=42).values
labels_test["random_modality"] = labels_test["modality"].sample(frac=1, random_state=42).values

# keep ONLY the train-multiome rows -- the pbmc3k label file also contains pbmc3k's own
# test_scATAC cells whose barcodes collide with our cross-platform ATAC (same pbmc3k ATAC).
train = pd.read_csv(TRAIN_ANNOT, index_col=0)
train = train[train["modality"] == "train multiomics"]
labels = pd.concat([train, labels_test], axis=0)
labels.to_csv(OUT)
print(f"wrote {OUT}: rna={rna.shape[0]} atac={atac.shape[0]} test={labels_test.shape[0]} total={labels.shape[0]}")
