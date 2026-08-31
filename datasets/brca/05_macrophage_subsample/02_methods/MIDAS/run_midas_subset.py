# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIDAS (scmidas) on the HT243 macrophage-subsample analysis -> per-(macrosub, sub) latent.csv.

Parameterized by env MACROSUB + SUB: HT243 TRAIN fixed, TEST = the macrophage subset (rna_sub{SUB} h5ad +
commonpeak_sub{SUB} 10x h5). Everything (saved_models, latent, runtime) writes under
{MACROSUB}/midas/sub{SUB}/, so the jobs are isolated and parallel-safe. Identical mosaic design to
HT243B1-S1H4/midas/run_midas.py; only the TEST inputs + OUT change.
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
SUBS = f"{ROOT}/BRCA/HT243-S1H4_subsample"
MACROSUB = os.environ["MACROSUB"]                 # e.g. HT243_S1H4_macrosub2
SUB      = os.environ["SUB"]                       # "1".."5"
TRAIN_H5  = f"{BRCA}/train_forHT243B1-S1H4.h5"      # fixed train (same as main HT243 run)
TEST_RNA  = f"{SUBS}/{MACROSUB}/HT243_S1H4_rna_sub{SUB}_rna.h5ad"
TEST_ATAC = f"{SUBS}/{MACROSUB}/HT243B1-S1H4_commonpeak_sub{SUB}.h5"
OUT = f"{SUBS}/{MACROSUB}/midas/sub{SUB}"
os.makedirs(OUT, exist_ok=True)
SEED = int(os.environ.get("SEED", "420"))
EPOCHS = int(os.environ.get("MIDAS_EPOCHS", "2000"))
start = timeit.default_timer()
print(f"[MIDAS subset] {MACROSUB} sub{SUB} -> {OUT}", flush=True)

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
    if "counts" in ad.layers:
        ad.X = ad.layers["counts"]
    elif ad.raw is not None and ad.raw.shape[1] == ad.shape[1]:
        ad.X = ad.raw.X
    x = ad.X[:50].toarray() if hasattr(ad.X, "toarray") else np.asarray(ad.X[:50])
    if not np.allclose(x, np.round(x)):
        print(f"  WARNING: {ad.shape} .X is not integer counts -> MIDAS expects raw counts")
    return ad


def load_test(path, suffix):
    ad = sc.read_h5ad(path); ad.var_names_make_unique(); ad = raw_counts(ad)
    ad.obs_names = [f"{b}_{suffix}" for b in ad.obs_names]
    return ad


def load_test_atac(path, suffix):
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

g = tr_rna.var_names.intersection(te_rna.var_names)
p = tr_atac.var_names.intersection(te_atac.var_names)
tr_rna, te_rna   = tr_rna[:, g].copy(),  te_rna[:, g].copy()
tr_atac, te_atac = tr_atac[:, p].copy(), te_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  test {te_rna.n_obs}+{te_atac.n_obs}", flush=True)

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
print(f"wrote {OUT}/latent.csv  shape={lat.shape}  atac_chunks={len(atac_chunks)}  time={stop-start:.1f}s")
