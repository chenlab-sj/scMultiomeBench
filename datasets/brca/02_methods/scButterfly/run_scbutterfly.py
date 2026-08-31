# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
scButterfly on BRCA HT243B1-S1H4 (unpaired RNA + ATAC) -> latent.csv for the benchmark.

Same Path-A design as the pbmc3k run: train the dual-aligned VAE on the paired TRAIN, then pull
the SHARED translator latent mean (mu) for the test cells, encoding each modality SEPARATELY (the
test pairing is never used for the embedding -> fair unpaired co-embedding). Only BRCA data
loading differs:
  * TRAIN = train_forHT243B1-S1H4.h5  (combined multiome 10x h5 -> split; paired, row-aligned)
  * TEST  = HT243_S1H4_rna.h5ad / HT243_S1H4_atac.h5ad  (modality split of the same cells, so
            aligned by barcode; matches the benchmark label.csv)
scButterfly's data structure is paired/row-aligned, so test RNA & ATAC are aligned by barcode for
the data object, but mu_r / mu_a are read independently -> the pairing isn't used in the embedding.

VERIFY on first run: (1) train h5 splits into Gene Expression + Peaks; (2) test h5ads are raw
counts; (3) preprocessing keeps all cells; (4) sum(chrom_list)==ATAC_data_p.n_vars;
(5) latent.csv ~ 5783 _rna + 5783 _atac.
"""
import os, re, timeit
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scButterfly.butterfly import Butterfly

ROOT = "/path/to/multiomeBench"
BRCA = f"{ROOT}/BRCA/HT243B1-S1H4"
TRAIN_H5  = f"{BRCA}/train_forHT243B1-S1H4.h5"
TEST_RNA  = f"{BRCA}/HT243_S1H4_rna.h5ad"
# common peaks = the SAME peak set as the train (chr:start-end). HT243_S1H4_atac.h5ad uses a
# DIFFERENT peak set (chr-start-end, different coords) -> 0 train overlap. Same 5783 test cells.
TEST_ATAC = f"{BRCA}/HT243B1-S1H4_commonpeaks.h5"
# reproducibility triplicate: REP=""/rep2/rep3 + SEED 420/0/42 -> scbutterfly/{,rep2,rep3}/latent.csv
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
OUT  = os.path.join(f"{BRCA}/scbutterfly", REP)
os.makedirs(OUT, exist_ok=True)
np.random.seed(SEED)
torch.manual_seed(SEED)
FORCE_SEED = os.environ.get("FORCE_SEED", "0") == "1"  # override scButterfly internal seed 19193
start = timeit.default_timer()


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
    # commonpeaks 10x h5: RAW counts, peaks named chr:start-end matching the train.
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


# ---- load train (paired) + test (modality split of same cells) ----
tr_rna, tr_atac = load_train(TRAIN_H5)
te_rna  = load_test(TEST_RNA)
te_atac = load_test_atac(TEST_ATAC)

# test RNA & ATAC are the same cells -> align by barcode (same order) for the paired data object
common = te_rna.obs_names.intersection(te_atac.obs_names)
te_rna, te_atac = te_rna[common].copy(), te_atac[common].copy()
# feature-align train<->test (shared genes / peaks)
g = tr_rna.var_names.intersection(te_rna.var_names)
p = tr_atac.var_names.intersection(te_atac.var_names)
tr_rna, te_rna   = tr_rna[:, g].copy(),  te_rna[:, g].copy()
tr_atac, te_atac = tr_atac[:, p].copy(), te_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  test {te_rna.n_obs}")

RNA_data  = sc.concat([tr_rna,  te_rna])
ATAC_data = sc.concat([tr_atac, te_atac])
ATAC_data = order_atac_by_chrom(ATAC_data)
n_train, n_test = tr_rna.n_obs, te_rna.n_obs
train_id = list(range(n_train))
test_id  = list(range(n_train, n_train + n_test))
test_bc  = list(RNA_data.obs_names[n_train:])

# ---- scButterfly (aug_type=None -> no scvi) ----
butterfly = Butterfly()
butterfly.load_data(RNA_data, ATAC_data, train_id, test_id, validation_id=None)
butterfly.data_preprocessing()
assert butterfly.RNA_data_p.n_obs == RNA_data.n_obs, "preprocessing dropped cells!"
chrom_list = chrom_counts(butterfly.ATAC_data_p.var)
assert sum(chrom_list) == butterfly.ATAC_data_p.n_vars, (sum(chrom_list), butterfly.ATAC_data_p.n_vars)
if FORCE_SEED:
    from scButterfly.split_datasets import setup_seed
    setup_seed(SEED)                              # override the internal 19193 for weight init
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
    butterfly.train_model(seed=SEED) if FORCE_SEED else butterfly.train_model()

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
print(f"wrote latent.csv  shape={lat.shape}  chroms={len(chrom_list)}  time={stop-start:.1f}s")
