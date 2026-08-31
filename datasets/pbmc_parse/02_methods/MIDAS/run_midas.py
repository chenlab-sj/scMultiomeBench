# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIDAS (scmidas 0.3.0) on the Parse CROSS-PLATFORM dataset -> latent.csv for benchmark_metrics.py.

This is the pbmc3k MIDAS run adapted to the cross-platform (Parse) setup. The MODELING
logic / params are mirrored verbatim from
  pbmc/pbmc3k/scripts/midas/run_midas.py
Only the INPUTS changed, because here RNA and ATAC are NOT the same split multiome:

  * TRAIN  = paired multiome reference (RNA+ATAC, SAME cells) -> the anchor batch.
             10x mtx dir: input_azimuth/pbmc3k_filtered_feature_bc_matrix_train/
  * test RNA  = independent Parse whole-cell RNA (DIFFERENT cells, RNA only).
             10x mtx dir: input_azimuth/pbmc_parse_d1/        -> obs '<bc>_rna'
  * test ATAC = pbmc3k test ATAC peaks (DIFFERENT cells, ATAC only).
             10x mtx dir: input_azimuth/pbmc3k_test_atac/     -> obs '<bc>_atac'

So instead of one combined test h5 (split into RNA/ATAC of the SAME cells, as in pbmc3k),
the two query batches are read from two separate 10x dirs with unrelated barcodes. RNA
genes are reconciled by sc.concat (inner join on gene symbols; Parse genes are a subset of
the train genes). ATAC peaks are identical between train and test (same peak set), so the
ATAC concat aligns peak-for-peak. Everything downstream (chromosome chunking, batch labels,
dims_x, MIDAS train, joint latent, test-cell filtering, latent.csv) is unchanged.

MIDAS does mosaic integration: a MuData of per-modality AnnDatas + a batch_key, where
different batches measure different modality combos. We use 3 batches -- paired TRAIN
(anchor) + unpaired test RNA + unpaired test ATAC -- then read the joint latent for the
test cells (<bc>_rna / <bc>_atac).

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
INPUT = f"{ROOT}/pbmc_parse/input_azimuth"
# NB: read the .h5 (not the mtx dirs) -- the mtx dirs hold uncompressed matrix.mtx, which
# scanpy's read_10x_mtx rejects (wants matrix.mtx.gz). The .h5 files are equivalent and match
# the original pbmc3k MIDAS reader.
TRAIN_DIR     = f"{INPUT}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_train.h5"  # paired multiome (RNA+ATAC)
TEST_RNA_DIR  = f"{INPUT}/pbmc_parse_d1.h5"                          # Parse RNA only
TEST_ATAC_DIR = f"{INPUT}/pbmc3k_test_atac.h5"                       # pbmc3k ATAC only
OUT  = f"{ROOT}/pbmc_parse/MIDAS"
os.makedirs(OUT, exist_ok=True)
start = timeit.default_timer()

EPOCHS = int(os.environ.get("MIDAS_EPOCHS", "2000"))  # real run=2000; smoke: export MIDAS_EPOCHS=10 before bsub


def load_split(d, suffix=None):
    """Read a 10x mtx dir and split into (rna, atac) by feature_types (paired multiome)."""
    a = sc.read_10x_h5(d, gex_only=False); a.var_names_make_unique()
    rna  = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    atac = a[:, a.var["feature_types"] == "Peaks"].copy()
    if suffix:
        rna.obs_names  = [f"{b}_rna"  for b in rna.obs_names]
        atac.obs_names = [f"{b}_atac" for b in atac.obs_names]
    return rna, atac


def load_modality(d, feature_type, suffix):
    """Read a single-modality 10x mtx dir (test RNA = Parse, or test ATAC = pbmc3k),
    keep only the requested feature_type, and suffix the barcodes (_rna / _atac).
    Used for the UNPAIRED cross-platform query batches (different cells, one modality)."""
    a = sc.read_10x_h5(d, gex_only=False); a.var_names_make_unique()
    a = a[:, a.var["feature_types"] == feature_type].copy()
    a.obs_names = [f"{b}_{suffix}" for b in a.obs_names]
    return a


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


tr_rna, tr_atac = load_split(TRAIN_DIR)                              # paired anchor (shared obs_names)
te_rna  = load_modality(TEST_RNA_DIR,  "Gene Expression", "rna")     # Parse RNA only  (<bc>_rna)
te_atac = load_modality(TEST_ATAC_DIR, "Peaks",           "atac")    # pbmc3k ATAC only (<bc>_atac)

# RNA modality = train + test-rna ; ATAC modality = train + test-atac (RAW COUNTS, no norm/log).
# sc.concat inner-joins on var (gene symbols): Parse RNA genes are a subset of the train genes,
# and the ATAC peak set is identical between train and test, so concat aligns features.
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
