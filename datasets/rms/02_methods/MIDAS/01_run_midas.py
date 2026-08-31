# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIDAS (scmidas) on RMS Mast607A (UNPAIRED RNA + ATAC) -> latent.csv for the benchmark/reproducibility.

Same mosaic design as pbmc3k/BRCA. RMS data is CLEANER than BRCA: both train and test are combined
multiome 10x h5 (Gene Expression + common Peaks in one file), so we read_10x_h5 + split both.
  TRAIN = Mast39_Mast213F_..._Mast607A.h5   (paired anchor; what cobolt/scVI use)
  TEST  = ..._Mast607A_TB19_22652_commonpeaks.h5  (RNA + common peaks; unpaired query)
Both are RAW counts (MIDAS log1p's RNA + binarizes ATAC internally). Train/test features intersected.
Reproducibility triplicate via env SEED/REP (420/0/42 -> midas/{,rep2,rep3}/latent.csv).

VERIFY on first run: both h5 split into Gene Expression + Peaks; shared genes/peaks not near-zero;
latent ~ 8055 _rna + 8055 _atac (matches label.csv; the test h5 has ~8291, extra cells dropped downstream).
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
RMS  = f"{ROOT}/RMS/Mast607/data"
TRAIN_H5 = f"{RMS}/Mast39_Mast213F_train/Mast39_Mast213F_filtered_feature_bc_matrix_Mast607A.h5"
TEST_H5  = f"{RMS}/Mast607A_TB19_22652/filtered_feature_bc_matrix_Mast607A_TB19_22652_commonpeaks.h5"

# reproducibility triplicate: REP=""/rep2/rep3 + SEED 420/0/42 -> midas/{,rep2,rep3}/latent.csv
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
OUT  = os.path.join(f"{ROOT}/RMS/Mast607/script/midas", REP)
os.makedirs(OUT, exist_ok=True)
start = timeit.default_timer()
EPOCHS = int(os.environ.get("MIDAS_EPOCHS", "2000"))

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


tr_rna, tr_atac = load_split(TRAIN_H5)               # paired anchor (no suffix)
te_rna, te_atac = load_split(TEST_H5, suffix=True)   # unpaired test -> rna gets _rna, atac gets _atac
g = tr_rna.var_names.intersection(te_rna.var_names)
p = tr_atac.var_names.intersection(te_atac.var_names)
tr_rna, te_rna   = tr_rna[:, g].copy(),  te_rna[:, g].copy()
tr_atac, te_atac = tr_atac[:, p].copy(), te_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  test {te_rna.n_obs}+{te_atac.n_obs}")

rna  = sc.concat([tr_rna,  te_rna])
atac = sc.concat([tr_atac, te_atac])
atac, atac_chunks = order_atac_by_chrom(atac)

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
