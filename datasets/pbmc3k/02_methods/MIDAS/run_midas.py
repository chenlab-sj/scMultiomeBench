# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIDAS (scmidas 0.3.0) on pbmc3k (UNPAIRED RNA + ATAC) -> latent.csv for benchmark_metrics.py.

MIDAS does mosaic integration: a MuData of per-modality AnnDatas + a batch_key, where
different batches measure different modality combos. We use 3 batches -- paired TRAIN
(anchor) + unpaired test RNA + unpaired test ATAC -- then read the joint latent for the
test cells (<bc>_rna / <bc>_atac). Mirrors the Multigrate mosaic setup.

Data contract (scmidas 0.3.0, verified against src/scmidas/model.py + demo3/data_layout):
  * One AnnData per modality, RAW COUNTS in .X (NOT a layer). Do NOT normalize / log1p:
    MIDAS log1p-transforms RNA internally and BINARIZES ATAC internally (BERNOULLI decoder).
  * A 'batch' column must live in EACH modality's .obs (read per-modality, not mdata.obs).
  * Paired cells are detected by SHARED obs_names across modalities (the train anchor).
  * ATAC is a PEAK matrix encoded BY CHROMOSOME CHUNK; chunk sizes are declared via
    setup_mudata(..., dims_x={'rna': [n_rna], 'atac': [chunk1, chunk2, ...]}). Peaks must be
    ordered so each chromosome's features are contiguous (we sort ATAC vars by chromosome).

API (scmidas 0.3.0, demo3 mosaic RNA+ATAC):
  configs = load_config()
  scmidas.MIDAS.setup_mudata(mdata, batch_key='batch', dims_x={'rna':[...], 'atac':[...]})
  model = scmidas.MIDAS(mdata, configs=configs, save_model_path=...)
  model.train(max_epochs=2000, accelerator='gpu', devices=1)
  lat = model.get_latent_representation(kind='c')   # joint biological z_c, aligned to obs_names
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

ROOT = "/path/to/multiomeBench"
DATA = f"{ROOT}/pbmc/pbmc3k/Data"
TRAIN_H5 = f"{DATA}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_train.h5"
TEST_H5  = f"{DATA}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
# reproducibility triplicate: REP=""/rep2/rep3 + SEED 420/0/42 -> midas/{,rep2,rep3}/latent.csv
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
OUT  = os.path.join(f"{ROOT}/pbmc/pbmc3k/scripts/midas", REP)
os.makedirs(OUT, exist_ok=True)
start = timeit.default_timer()

EPOCHS = int(os.environ.get("MIDAS_EPOCHS", "2000"))  # real run=2000; smoke: export MIDAS_EPOCHS=10 before bsub

np.random.seed(SEED)
try:
    import torch; torch.manual_seed(SEED)
    import lightning.pytorch as pl; pl.seed_everything(SEED, workers=True)
except Exception:
    try:
        import pytorch_lightning as pl; pl.seed_everything(SEED)
    except Exception:
        pass


def load_split(h5, suffix=None):
    a = sc.read_10x_h5(h5, gex_only=False); a.var_names_make_unique()
    rna  = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    atac = a[:, a.var["feature_types"] == "Peaks"].copy()
    if suffix:
        rna.obs_names  = [f"{b}_rna"  for b in rna.obs_names]
        atac.obs_names = [f"{b}_atac" for b in atac.obs_names]
    return rna, atac


def chrom_of(var):
    """Best-effort chromosome label per ATAC peak (var index 'chr-start-end' / 'chr:start-end',
    or an 'interval'/'chrom' var column). Falls back to a single chunk if unparseable."""
    if "interval" in var.columns:
        src = var["interval"].astype(str)
    elif "chrom" in var.columns:
        return var["chrom"].astype(str).values
    else:
        src = pd.Series(var.index.astype(str), index=var.index)
    out = src.str.extract(r"^([^\s:_-]+)")[0]
    return out.fillna("NA").values


def order_atac_by_chrom(atac):
    """Sort peaks so each chromosome's features are contiguous; return sorted AnnData + chunk
    sizes (n peaks per chromosome) for dims_x. Chromosomes ordered by natural chr number."""
    chrom = chrom_of(atac.var)

    def chrom_key(c):
        m = re.match(r"^chr?(\w+)$", str(c), re.I)
        tok = m.group(1) if m else str(c)
        return (0, int(tok)) if tok.isdigit() else (1, str(tok))

    order = sorted(range(atac.n_vars), key=lambda i: (chrom_key(chrom[i]), i))
    atac = atac[:, order].copy()
    chrom = chrom[order]
    # contiguous run lengths == per-chromosome chunk sizes
    sizes, run, prev = [], 0, None
    for c in chrom:
        if prev is not None and c != prev:
            sizes.append(run); run = 0
        run += 1; prev = c
    if run:
        sizes.append(run)
    return atac, [int(s) for s in sizes]


tr_rna, tr_atac = load_split(TRAIN_H5)              # paired anchor (shared obs_names)
te_rna, te_atac = load_split(TEST_H5, suffix=True)  # unpaired query

# RNA modality = train + test-rna ; ATAC modality = train + test-atac (RAW COUNTS, no norm/log).
rna  = sc.concat([tr_rna,  te_rna])
atac = sc.concat([tr_atac, te_atac])

# ATAC: order peaks by chromosome and record chunk sizes for the published chunked encoder.
atac, atac_chunks = order_atac_by_chrom(atac)

# batch must live in EACH modality's .obs (setup_mudata reads per-modality, NOT mdata.obs).
# 'train' = the paired OVERLAP (same cells in both modalities); test cells appear in one
# modality only -> mosaic. This is the "train as overlap" design.
rna.obs["batch"]  = ["train"] * tr_rna.n_obs  + ["test_rna"]  * te_rna.n_obs
atac.obs["batch"] = ["train"] * tr_atac.n_obs + ["test_atac"] * te_atac.n_obs

# Ensure raw integer-count dtype in .X (MIDAS expects counts; it transforms internally).
for ad in (rna, atac):
    if hasattr(ad.X, "toarray"):
        pass  # keep sparse; MIDAS handles sparse counts
    ad.X = ad.X.astype(np.float32)

mdata = mu.MuData({"rna": rna, "atac": atac})

# dims_x: RNA is one chunk; ATAC is per-chromosome chunks (published architecture).
dims_x = {"rna": [int(rna.n_vars)], "atac": atac_chunks}

configs = load_config()
MIDAS.setup_mudata(mdata, batch_key="batch", dims_x=dims_x)
model = MIDAS(
    mdata,
    configs=configs,
    save_model_path=os.path.join(OUT, "saved_models"),
)
# train() forwards to a Lightning Trainer (no built-in default epochs) -> pass max_epochs.
# GPU node (BSUB -gpu num=1): accelerator='gpu', devices=1. ATAC binarized internally.
model.train(max_epochs=EPOCHS, accelerator="gpu", devices=1)

# Joint biological latent z_c, returned aligned to mdata.obs_names (test cells get a joint
# embedding via the cross-modal posterior even though each is single-modality).
lat = model.get_latent_representation(kind="c")

# ---- latent.csv for the TEST cells only ----
lat = pd.DataFrame(np.asarray(lat), index=list(mdata.obs_names))
test_bc = list(te_rna.obs_names) + list(te_atac.obs_names)
lat = lat.loc[lat.index.intersection(test_bc)]
# drop any all-NaN rows (cells absent from training would be NaN; should not happen here)
lat = lat.dropna(how="all")
lat.to_csv(os.path.join(OUT, "latent.csv"))

stop = timeit.default_timer()
with open(os.path.join(OUT, "midas_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={lat.shape}  atac_chunks={len(atac_chunks)}  time={stop-start:.1f}s")
