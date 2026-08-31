# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIDAS (scmidas) on BMMC_d1 s1d1 PAIRED multiome (RNA + ATAC, SAME s1d1 cells) -> latent.csv for the benchmark.

WITHOUT batch variant (the donor/site batch is IGNORED).
Mirrors BMMC scVI (01_run_scvi.py, batch_key="modality" only) vs scVI2 (adds
categorical_covariate_keys=["data"]). Here the MIDAS batch label marks ONLY the
modality-source (train vs test_rna vs test_atac), exactly like the RMS MIDAS script.
The donor/site (obs["data"]=test3) is NOT encoded in the batch_key, so
MIDAS does not correct for the donor batch.

Input contract mirrors BMMC scVI (read_10x_multiome):
  TRAIN = s1d2_train                 (paired multiome anchor; obs["data"]="train")
  TEST3 = train_test3_s1d1 -> _rna3/_atac3, obs["data"]="test3"  (s1d1 PAIRED baseline)
The test is split by var.modality (Gene Expression -> rna, Peaks -> atac):
te_rna = s1d1 RNA (_rna3), te_atac = s1d1 ATAC (_atac3) -- the SAME s1d1 cells.
(test1 s2d1 and test2 s4d1 are DROPPED for the paired baseline.)

Both train and test are RAW counts (MIDAS log1p's RNA + binarizes ATAC internally).
Train/test features intersected. Output latent.csv is indexed by the suffixed
barcodes (matching celltype.csv).
"""
import os
import re
import timeit
import numpy as np
import pandas as pd
import scanpy as sc
import mudata as mu
import scmidas
from scmidas.config import load_config

MIDAS = getattr(scmidas, "MIDAS", None) or scmidas.model.MIDAS


def read_10x_multiome(base_path):
    """env-safe VERBATIM copy of scvi.data.read_10x_multiome (uses only mmread + pandas + AnnData; scvi
    is absent in midas-env). The BMMC_d1 dirs are uncompressed v3 triplets (matrix.mtx/features.tsv/
    barcodes.tsv) that sc.read_10x_mtx can't parse but scvi (and this copy) can. Output matches the scVI
    baseline exactly: cells x features, var['modality'] from feature type, barcodes split on '-' (drops -1)."""
    import os
    from scipy.io import mmread
    from anndata import AnnData
    data = mmread(os.path.join(base_path, "matrix.mtx")).transpose()
    features = pd.read_csv(os.path.join(base_path, "features.tsv"), sep="\t", header=None, index_col=1)
    features.rename({0: "ID", 2: "modality", 3: "chr", 4: "start", 5: "end"}, axis="columns", inplace=True)
    features.index.name = None
    cell_annot = pd.read_csv(os.path.join(base_path, "barcodes.tsv"), sep="-", header=None, index_col=None)
    cell_annot.rename({0: "barcode", 1: "batch_id"}, axis="columns", inplace=True)
    cell_annot.set_index("barcode", inplace=True)
    cell_annot.index = cell_annot.index.astype(str)
    return AnnData(data.tocsr(), var=features, obs=cell_annot)

######################################################
## input args (mirror BMMC scVI 01_run_scvi.py)
train_10x = "/path/to/data/BMMC_d1/s1d2_train"           ## train data with combined peak
test3_10x = "/path/to/data/BMMC_d1/train_test3_s1d1"     ## test3 s1d1 (PAIRED multiome baseline)

out_dir = "/path/to/multiomeBench/BMMC_d1/s1d1_paired/midas/"
########################################################

SEED = int(os.environ.get("SEED", "420"))
EPOCHS = int(os.environ.get("MIDAS_EPOCHS", "2000"))
os.makedirs(out_dir, exist_ok=True)
start = timeit.default_timer()

np.random.seed(SEED)
try:
    import torch; torch.manual_seed(SEED)
    import lightning.pytorch as pl; pl.seed_everything(SEED, workers=True)
except Exception:
    try:
        import pytorch_lightning as pl; pl.seed_everything(SEED)
    except Exception:
        pass


def chrom_of(var):
    if "interval" in var.columns:
        src = var["interval"].astype(str)
    elif "chrom" in var.columns:
        return var["chrom"].astype(str).values
    else:
        src = pd.Series(var.index.astype(str), index=var.index)
    return src.str.extract(r"^([^\s:_-]+)")[0].fillna("NA").values


def order_atac_by_chrom(atac):
    chrom = chrom_of(atac.var)

    def chrom_key(c):
        m = re.match(r"^chr?(\w+)$", str(c), re.I); tok = m.group(1) if m else str(c)
        return (0, int(tok)) if tok.isdigit() else (1, str(tok))

    order = sorted(range(atac.n_vars), key=lambda i: (chrom_key(chrom[i]), i))
    atac = atac[:, order].copy(); chrom = chrom[order]
    sizes, run, prev = [], 0, None
    for c in chrom:
        if prev is not None and c != prev:
            sizes.append(run); run = 0
        run += 1; prev = c
    if run:
        sizes.append(run)
    return atac, [int(s) for s in sizes]


## ---- load (read_10x_multiome; var.modality Gene Expression / Peaks) ----
adata_train = read_10x_multiome(train_10x)
adata_train.var_names_make_unique()
adata_test3 = read_10x_multiome(test3_10x)
adata_test3.var_names_make_unique()

## split train into paired rna/atac (no suffix; paired anchor)
tr_rna  = adata_train[:, adata_train.var.modality == "Gene Expression"].copy()
tr_atac = adata_train[:, adata_train.var.modality == "Peaks"].copy()
tr_rna.obs["data"] = "train"
tr_atac.obs["data"] = "train"

## split the s1d1 test by modality, suffix barcodes, tag donor batch via obs["data"]
adata_test3_rna  = adata_test3[:, adata_test3.var.modality == "Gene Expression"].copy()
adata_test3_atac = adata_test3[:, adata_test3.var.modality == "Peaks"].copy()

adata_test3_rna.obs.index  = adata_test3_rna.obs.index  + "_rna3";  adata_test3_rna.obs["data"]  = "test3"
adata_test3_atac.obs.index = adata_test3_atac.obs.index + "_atac3"; adata_test3_atac.obs["data"] = "test3"

## s1d1 PAIRED multiome: te_rna = s1d1 test3 RNA (_rna3), te_atac = s1d1 test3 ATAC (_atac3) -- SAME cells
te_rna  = adata_test3_rna
te_atac = adata_test3_atac

## intersect train/test features per modality
g = tr_rna.var_names.intersection(te_rna.var_names)
p = tr_atac.var_names.intersection(te_atac.var_names)
tr_rna, te_rna   = tr_rna[:, g].copy(),  te_rna[:, g].copy()
tr_atac, te_atac = tr_atac[:, p].copy(), te_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  test {te_rna.n_obs}+{te_atac.n_obs}")

## build the mosaic: rna = train + all test_rna ; atac = train + all test_atac
rna  = sc.concat([tr_rna,  te_rna])
atac = sc.concat([tr_atac, te_atac])
atac, atac_chunks = order_atac_by_chrom(atac)

## ----------------------------------------------------------------------------
## WITHOUT batch: the MIDAS batch label marks ONLY modality-source
## (train vs test_rna vs test_atac). The donor/site (test3) is NOT encoded,
## so MIDAS does not correct for the donor batch. (mirrors scVI batch_key="modality")
## ----------------------------------------------------------------------------
rna.obs["batch"]  = ["train"] * tr_rna.n_obs  + ["test_rna"]  * te_rna.n_obs
atac.obs["batch"] = ["train"] * tr_atac.n_obs + ["test_atac"] * te_atac.n_obs
for ad in (rna, atac):
    ad.X = ad.X.astype(np.float32)

mdata = mu.MuData({"rna": rna, "atac": atac})
dims_x = {"rna": [int(rna.n_vars)], "atac": atac_chunks}

configs = load_config()
MIDAS.setup_mudata(mdata, batch_key="batch", dims_x=dims_x)
model = MIDAS(mdata, configs=configs, save_model_path=os.path.join(out_dir, "saved_models"))
model.train(max_epochs=EPOCHS, accelerator="gpu", devices=1)

lat = model.get_latent_representation(kind="c")
lat = pd.DataFrame(np.asarray(lat), index=list(mdata.obs_names))
test_bc = list(te_rna.obs_names) + list(te_atac.obs_names)
lat = lat.loc[lat.index.intersection(test_bc)].dropna(how="all")
lat.to_csv(os.path.join(out_dir, "latent.csv"))

stop = timeit.default_timer()
with open(os.path.join(out_dir, "midas_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={lat.shape}  atac_chunks={len(atac_chunks)}  time={stop-start:.1f}s")
