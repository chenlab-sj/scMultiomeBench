# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIDAS (scmidas) on BRCA HT243B1-S1H4 (UNPAIRED RNA + ATAC) -> latent.csv for the benchmark.

Same mosaic design as the pbmc3k MIDAS run (3 batches: paired TRAIN anchor + unpaired test RNA
+ unpaired test ATAC; read joint z_c for the test cells). Only the BRCA data loading differs:
  * TRAIN  = train_forHT243B1-S1H4.h5  -> combined multiome 10x h5, split by feature_types (paired)
  * TEST   = HT243_S1H4_rna.h5ad / HT243_S1H4_atac.h5ad  (the same files scglue/scDART use, so the
             barcodes match the benchmark label.csv)
MIDAS needs RAW COUNTS (it log1p's RNA + binarizes ATAC internally) -> raw_counts() prefers a
counts layer / .raw. Train & test features are intersected so RNA/ATAC dims are consistent.

VERIFY on first run: (1) the train h5 splits into Gene Expression + Peaks; (2) the test h5ads are
RAW counts (printed integer-check); (3) shared-gene / shared-peak counts are sane (not near-zero);
(4) latent.csv ~ 5783 _rna + 5783 _atac (matches label.csv).
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
BRCA = f"{ROOT}/BRCA/HT243B1-S1H4"
TRAIN_H5  = f"{BRCA}/train_forHT243B1-S1H4.h5"   # paired anchor (combined multiome)
TEST_RNA  = f"{BRCA}/HT243_S1H4_rna.h5ad"        # unpaired test RNA  (matches label.csv)
# common peaks = SAME peak set as the train (chr:start-end); HT243_S1H4_atac.h5ad uses a DIFFERENT
# peak set (chr-start-end, different coords) -> 0 train overlap. Same 5783 test cells.
TEST_ATAC = f"{BRCA}/HT243B1-S1H4_commonpeaks.h5"  # unpaired test ATAC (10x h5, raw counts)
# reproducibility triplicate: REP=""/rep2/rep3 + SEED 420/0/42 -> midas/{,rep2,rep3}/latent.csv
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
OUT  = os.path.join(f"{BRCA}/midas", REP)
os.makedirs(OUT, exist_ok=True)
start = timeit.default_timer()
EPOCHS = int(os.environ.get("MIDAS_EPOCHS", "2000"))  # smoke: export MIDAS_EPOCHS=10 before bsub

np.random.seed(SEED)
try:
    import torch; torch.manual_seed(SEED)
    import lightning.pytorch as pl; pl.seed_everything(SEED, workers=True)
except Exception:
    try:
        import pytorch_lightning as pl; pl.seed_everything(SEED)
    except Exception:
        pass


def load_train(h5):
    a = sc.read_10x_h5(h5, gex_only=False); a.var_names_make_unique()
    return (a[:, a.var["feature_types"] == "Gene Expression"].copy(),
            a[:, a.var["feature_types"] == "Peaks"].copy())


def raw_counts(ad):
    """Return an AnnData whose .X is RAW counts (MIDAS transforms internally)."""
    if "counts" in ad.layers:
        ad.X = ad.layers["counts"]
    elif ad.raw is not None and ad.raw.shape[1] == ad.shape[1]:
        ad.X = ad.raw.X
    x = ad.X[:50].toarray() if hasattr(ad.X, "toarray") else np.asarray(ad.X[:50])
    if not np.allclose(x, np.round(x)):
        print(f"  WARNING: {ad.shape} .X is not integer counts -> MIDAS expects raw counts; check the h5ad")
    return ad


def load_test(path, suffix):
    ad = sc.read_h5ad(path); ad.var_names_make_unique(); ad = raw_counts(ad)
    ad.obs_names = [f"{b}_{suffix}" for b in ad.obs_names]
    return ad


def load_test_atac(path, suffix):
    # commonpeaks 10x h5: RAW counts, peaks chr:start-end matching the train.
    a = sc.read_10x_h5(path, gex_only=False); a.var_names_make_unique(); a = raw_counts(a)
    a.obs_names = [f"{b}_{suffix}" for b in a.obs_names]
    return a


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


tr_rna, tr_atac = load_train(TRAIN_H5)
te_rna  = load_test(TEST_RNA, "rna")
te_atac = load_test_atac(TEST_ATAC, "atac")

# align features train<->test so RNA/ATAC dims are consistent for the mosaic concat
g = tr_rna.var_names.intersection(te_rna.var_names)
p = tr_atac.var_names.intersection(te_atac.var_names)
tr_rna, te_rna   = tr_rna[:, g].copy(),  te_rna[:, g].copy()
tr_atac, te_atac = tr_atac[:, p].copy(), te_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  test {te_rna.n_obs}+{te_atac.n_obs}")

rna  = sc.concat([tr_rna,  te_rna])
atac = sc.concat([tr_atac, te_atac])
atac, atac_chunks = order_atac_by_chrom(atac)

# batch must live in EACH modality's .obs; 'train' = the paired overlap (same cells both modalities)
rna.obs["batch"]  = ["train"] * tr_rna.n_obs  + ["test_rna"]  * te_rna.n_obs
atac.obs["batch"] = ["train"] * tr_atac.n_obs + ["test_atac"] * te_atac.n_obs
for ad in (rna, atac):
    ad.X = ad.X.astype(np.float32)

mdata = mu.MuData({"rna": rna, "atac": atac})
dims_x = {"rna": [int(rna.n_vars)], "atac": atac_chunks}

configs = load_config()
MIDAS.setup_mudata(mdata, batch_key="batch", dims_x=dims_x)
model = MIDAS(mdata, configs=configs, save_model_path=os.path.join(OUT, "saved_models"))
model.train(max_epochs=EPOCHS, accelerator="gpu", devices=1)

lat = model.get_latent_representation(kind="c")
lat = pd.DataFrame(np.asarray(lat), index=list(mdata.obs_names))
test_bc = list(te_rna.obs_names) + list(te_atac.obs_names)
lat = lat.loc[lat.index.intersection(test_bc)].dropna(how="all")
lat.to_csv(os.path.join(OUT, "latent.csv"))

stop = timeit.default_timer()
with open(os.path.join(OUT, "midas_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={lat.shape}  atac_chunks={len(atac_chunks)}  time={stop-start:.1f}s")
