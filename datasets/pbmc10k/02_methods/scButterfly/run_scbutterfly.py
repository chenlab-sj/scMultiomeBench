# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
scButterfly on pbmc10k (unpaired RNA + ATAC) -> latent.csv for benchmark_metrics.py.

scButterfly is a dual-aligned VAE TRANSLATION method; its public API returns translated
profiles, not a latent. Path A: train on paired data, then pull the SHARED translator
latent mean (mu) for the test cells, embedding each modality SEPARATELY (so the test
pairing is never used -> a fair unpaired co-embedding):
    R2 = model.RNA_encoder(x);   _,_, mu_r,_ = model.translator.test_model(R2, 'RNA')
    A2 = model.ATAC_encoder(x);  _,_, mu_a,_ = model.translator.test_model(A2, 'ATAC')

Grounded against the source:
  * data_preprocessing() -> butterfly.RNA_data_p / butterfly.ATAC_data_p (AnnData): RNA = 3000
    HVG, normalized+log; ATAC = binarized/TF-IDF/[0,1]. These .X are the encoder inputs (NOT
    raw counts). Cell rows are preserved (preprocessing filters FEATURES, not cells).
  * We include the test cells in load_data so they get the SAME HVG/peak preprocessing as
    train; we train only on train_id; aug_type defaults to None -> the broken scvi is never
    touched (scButterfly imports scvi lazily only for MultiVI_augmentation).
  * construct_model(chrom_list) needs peak counts per chromosome of the PROCESSED ATAC
    (Split_Chrom_Encoder). We sort peaks by chromosome before load_data and compute
    chrom_list from ATAC_data_p afterward.

VERIFY on first run: (1) preprocessing keeps all cells (RNA_data_p.n_obs == RNA_data.n_obs);
(2) sum(chrom_list) == ATAC_data_p.n_vars; (3) RNA_encoder accepts RNA_data_p.X width.
"""
import os, re, timeit
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scButterfly.butterfly import Butterfly

ROOT = "/path/to/multiomeBench"
DATA = "/path/to/data/pbmc10k/pbmc10k_data_design0719"
TRAIN_H5 = f"{DATA}/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_train.h5"
TEST_H5  = f"{DATA}/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5"
# reproducibility triplicate: REP=""/rep2/rep3 + SEED 420/0/42 -> scbutterfly/{,rep2,rep3}/latent.csv
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
OUT  = os.path.join(f"{ROOT}/pbmc/pbmc10k/scripts/scbutterfly", REP)
os.makedirs(OUT, exist_ok=True)
np.random.seed(SEED)
torch.manual_seed(SEED)
# FORCE_SEED=1 overrides scButterfly's hard-coded internal seed (19193). Butterfly() and train_model()
# both call setup_seed(19193), discarding the seeds above -> triplicates are byte-identical. With
# FORCE_SEED=1 we re-seed with SEED before construct_model (weight init) AND pass seed=SEED to
# train_model, so the external SEED genuinely varies the run.
FORCE_SEED = os.environ.get("FORCE_SEED", "0") == "1"
start = timeit.default_timer()

def load_split(h5):
    a = sc.read_10x_h5(h5, gex_only=False); a.var_names_make_unique()
    rna  = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    atac = a[:, a.var["feature_types"] == "Peaks"].copy()
    return rna, atac

def chrom_key(c):
    m = re.match(r"^chr?(\w+)$", str(c), re.I); t = m.group(1) if m else str(c)
    return (0, int(t)) if t.isdigit() else (1, str(t))

def chrom_of(var):
    iv = var["interval"].astype(str) if "interval" in var else pd.Series(var.index.astype(str), index=var.index)
    return iv.str.extract(r"^([^\s:_-]+)")[0].fillna("NA").values

def order_atac_by_chrom(atac):
    chrom = chrom_of(atac.var)
    order = sorted(range(atac.n_vars), key=lambda i: (chrom_key(chrom[i]), i))
    atac = atac[:, order].copy()
    atac.var["chrom"] = chrom[order]
    return atac

def chrom_counts(var):
    """consecutive per-chromosome peak counts of an already chrom-ordered var."""
    chrom = var["chrom"].values if "chrom" in var else chrom_of(var)
    sizes, run, prev = [], 0, None
    for c in chrom:
        if prev is not None and c != prev: sizes.append(run); run = 0
        run += 1; prev = c
    if run: sizes.append(run)
    return [int(s) for s in sizes]

# ---- build paired data: all cells = train + test (raw counts) ----
tr_rna, tr_atac = load_split(TRAIN_H5)
te_rna, te_atac = load_split(TEST_H5)
RNA_data  = sc.concat([tr_rna,  te_rna])
ATAC_data = sc.concat([tr_atac, te_atac])
ATAC_data = order_atac_by_chrom(ATAC_data)        # peaks grouped by chromosome
n_train, n_test = tr_rna.n_obs, te_rna.n_obs
train_id = list(range(n_train))
test_id  = list(range(n_train, n_train + n_test))
test_bc  = list(RNA_data.obs_names[n_train:])     # raw test barcodes

# ---- scButterfly (aug_type=None -> no scvi) ----
butterfly = Butterfly()
butterfly.load_data(RNA_data, ATAC_data, train_id, test_id, validation_id=None)
butterfly.data_preprocessing()                    # -> RNA_data_p / ATAC_data_p
assert butterfly.RNA_data_p.n_obs == RNA_data.n_obs, "preprocessing dropped cells!"
chrom_list = chrom_counts(butterfly.ATAC_data_p.var)
assert sum(chrom_list) == butterfly.ATAC_data_p.n_vars, (sum(chrom_list), butterfly.ATAC_data_p.n_vars)
if FORCE_SEED:
    from scButterfly.split_datasets import setup_seed
    setup_seed(SEED)                              # override the internal 19193 for weight init
butterfly.construct_model(chrom_list)

# Reuse the trained model if it's already saved (./model/*.pt) -> skip the ~20min retrain.
# scButterfly saves each component's state_dict at the end of train_model().
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

# ---- extract shared latent mu for the TEST cells (RNA + ATAC separately) ----
# scButterfly's Model is a custom wrapper (NOT an nn.Module) -> put the sub-modules in
# eval mode and read the device off RNA_encoder's params.
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
