# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
scButterfly on the HT243 macrophage-subsample analysis -> per-(macrosub, sub) latent.csv.

Parameterized by env MACROSUB (HT243_S1H4_macrosub / macrosub2..5) + SUB (1..5): the HT243 TRAIN is fixed,
only the TEST changes to the macrophage subset (rna_sub{SUB} h5ad + commonpeak_sub{SUB} 10x h5). EVERY
output (model, latent, runtime) writes under {MACROSUB}/scbutterfly/sub{SUB}/, so the 25 jobs are fully
isolated and safe to run in parallel. Identical Path-A design to HT243B1-S1H4/scbutterfly/run_scbutterfly.py
(train the dual VAE on the paired train, pull the shared translator mu per modality) -- only inputs+OUT differ.
"""
import os, re, timeit
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scButterfly.butterfly import Butterfly

ROOT = "/path/to/multiomeBench"
BRCA = f"{ROOT}/BRCA/HT243B1-S1H4"
SUBS = f"{ROOT}/BRCA/HT243-S1H4_subsample"
MACROSUB = os.environ["MACROSUB"]                 # e.g. HT243_S1H4_macrosub2
SUB      = os.environ["SUB"]                       # "1".."5"
TRAIN_H5 = f"{BRCA}/train_forHT243B1-S1H4.h5"       # fixed train (same as main HT243 run)
TEST_RNA  = f"{SUBS}/{MACROSUB}/HT243_S1H4_rna_sub{SUB}_rna.h5ad"
TEST_ATAC = f"{SUBS}/{MACROSUB}/HT243B1-S1H4_commonpeak_sub{SUB}.h5"
OUT = f"{SUBS}/{MACROSUB}/scbutterfly/sub{SUB}"
os.makedirs(OUT, exist_ok=True)
# CRITICAL for parallel jobs: scButterfly's train_model() saves/loads its checkpoints to a RELATIVE
# ./model/ (the cwd). Without this chdir, all 25 jobs share one ./model at the submission dir and
# cross-load each other's (different chrom_list -> size-mismatch RuntimeError). chdir(OUT) makes ./model
# == OUT/model, isolated per (macrosub,sub). All inputs above are absolute, so this is safe.
os.chdir(OUT)
SEED = int(os.environ.get("SEED", "420"))
np.random.seed(SEED)
torch.manual_seed(SEED)
start = timeit.default_timer()
print(f"[scButterfly subset] {MACROSUB} sub{SUB} -> {OUT}", flush=True)


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
        print(f"  WARNING: {ad.shape} .X not integer counts -> scButterfly expects raw counts")
    return ad


def load_test(path):
    ad = sc.read_h5ad(path); ad.var_names_make_unique()
    return raw_counts(ad)


def load_test_atac(path):
    a = sc.read_10x_h5(path, gex_only=False); a.var_names_make_unique()
    return raw_counts(a)


def chrom_key(c):
    m = re.match(r"^chr?(\w+)$", str(c), re.I); t = m.group(1) if m else str(c)
    return (0, int(t)) if t.isdigit() else (1, str(t))


def chrom_of(var):
    iv = var["interval"].astype(str) if "interval" in var else pd.Series(var.index.astype(str), index=var.index)
    return iv.str.extract(r"^([^\s:_-]+)")[0].fillna("NA").values


def order_atac_by_chrom(atac):
    chrom = chrom_of(atac.var)
    order = sorted(range(atac.n_vars), key=lambda i: (chrom_key(chrom[i]), i))
    atac = atac[:, order].copy(); atac.var["chrom"] = chrom[order]
    return atac


def chrom_counts(var):
    chrom = var["chrom"].values if "chrom" in var else chrom_of(var)
    sizes, run, prev = [], 0, None
    for c in chrom:
        if prev is not None and c != prev: sizes.append(run); run = 0
        run += 1; prev = c
    if run: sizes.append(run)
    return [int(s) for s in sizes]


# ---- load train (paired) + test (the macrophage subset) ----
tr_rna, tr_atac = load_train(TRAIN_H5)
te_rna  = load_test(TEST_RNA)
te_atac = load_test_atac(TEST_ATAC)

common = te_rna.obs_names.intersection(te_atac.obs_names)
te_rna, te_atac = te_rna[common].copy(), te_atac[common].copy()
g = tr_rna.var_names.intersection(te_rna.var_names)
p = tr_atac.var_names.intersection(te_atac.var_names)
tr_rna, te_rna   = tr_rna[:, g].copy(),  te_rna[:, g].copy()
tr_atac, te_atac = tr_atac[:, p].copy(), te_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  test {te_rna.n_obs}", flush=True)

RNA_data  = sc.concat([tr_rna,  te_rna])
ATAC_data = sc.concat([tr_atac, te_atac])
ATAC_data = order_atac_by_chrom(ATAC_data)
n_train, n_test = tr_rna.n_obs, te_rna.n_obs
train_id = list(range(n_train))
test_id  = list(range(n_train, n_train + n_test))
test_bc  = list(RNA_data.obs_names[n_train:])

butterfly = Butterfly()
butterfly.load_data(RNA_data, ATAC_data, train_id, test_id, validation_id=None)
butterfly.data_preprocessing()
assert butterfly.RNA_data_p.n_obs == RNA_data.n_obs, "preprocessing dropped cells!"
chrom_list = chrom_counts(butterfly.ATAC_data_p.var)
assert sum(chrom_list) == butterfly.ATAC_data_p.n_vars, (sum(chrom_list), butterfly.ATAC_data_p.n_vars)
butterfly.construct_model(chrom_list)

model = butterfly.model
ckpt = os.path.join(OUT, "model")
if os.path.exists(os.path.join(ckpt, "RNA_encoder.pt")):
    print("found saved scButterfly model -> loading (skipping retrain)")
    for name, mod in (("RNA_encoder", model.RNA_encoder),
                      ("ATAC_encoder", model.ATAC_encoder),
                      ("translator",  model.translator)):
        mod.load_state_dict(torch.load(os.path.join(ckpt, f"{name}.pt"), map_location="cpu"))
else:
    butterfly.train_model()

for mod in (model.RNA_encoder, model.ATAC_encoder, model.translator):
    mod.eval()
dev = next(model.RNA_encoder.parameters()).device


def embed(adata_p, idx, encoder, modality):
    X = adata_p.X[idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    out = []
    with torch.no_grad():
        for i in range(0, X.shape[0], 512):
            xb = torch.tensor(X[i:i+512], dtype=torch.float32, device=dev)
            _, _, mu, _ = model.translator.test_model(encoder(xb), modality)
            out.append(mu.cpu().numpy())
    return np.concatenate(out, 0)


mu_r = embed(butterfly.RNA_data_p,  test_id, model.RNA_encoder,  "RNA")
mu_a = embed(butterfly.ATAC_data_p, test_id, model.ATAC_encoder, "ATAC")

lat = pd.concat([pd.DataFrame(mu_r, index=[f"{b}_rna"  for b in test_bc]),
                 pd.DataFrame(mu_a, index=[f"{b}_atac" for b in test_bc])])
lat.to_csv(os.path.join(OUT, "latent.csv"))
stop = timeit.default_timer()
with open(os.path.join(OUT, "scbutterfly_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote {OUT}/latent.csv  shape={lat.shape}  chroms={len(chrom_list)}  time={stop-start:.1f}s")
