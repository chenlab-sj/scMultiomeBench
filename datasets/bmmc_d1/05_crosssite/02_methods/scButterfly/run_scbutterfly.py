# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
scButterfly on BMMC_d1 (multi-batch, unpaired RNA + ATAC) -> latent.csv. Path-A: train the
dual-aligned VAE on the paired TRAIN (s1d2), then pull the shared translator latent mean (mu) for the
test cells, encoding each modality SEPARATELY (test pairing never used for the embedding).

INPUT contract mirrors BMMC_d1 scVI (01_run_scvi.py) EXACTLY: read_10x_multiome on the train 10x-mtx dir +
3 test dirs, split each test by var.modality into RNA ("Gene Expression") and ATAC ("Peaks"), suffix the
barcodes (test1 -> _rna1/_atac1, test2 -> _rna2/_atac2, test3 -> _rna3/_atac3), set obs["data"], and concat
the 3 test RNAs into adata_rna and the 3 test ATACs into adata_atac (anndata.concat, axis=0, .var = test1's).
Each test dir is itself a paired multiome (same cells/order in its RNA & ATAC halves), so adata_rna and
adata_atac are ROW-ALIGNED -- same train/test index logic as the RMS 01_run_scbutterfly.py (RMS read one
paired multiome h5; here the "paired test" is the 3 concatenated test batches).

scButterfly treats the paired TRAIN as train_id and the (row-aligned) concatenated test RNA+ATAC as the
unpaired test set (test_id). NO dataset-batch parameter exists in scButterfly -> single run across all 3
test batches (do NOT add any batch_key). Reproducibility via env SEED/REP.

SKIP-RETRAIN TRAP (see scbutterfly-setup memory): reuses model/*.pt if present. Each rep has its OWN
model dir (OUT/model). If a run is killed mid-train, DELETE that rep's model/ before re-running, or it
loads a partial model -> degenerate latent. A real train is minutes; a tiny "time=" means it skipped.

OUTPUT: latent.csv indexed by the suffixed barcodes (<bc>_rna1 .. <bc>_atac3) so it matches
BMMC_d1/celltype.csv. out_dir = this scbutterfly folder.
"""
import os, re, timeit
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import anndata
from scButterfly.butterfly import Butterfly


def read_10x_multiome(base_path):
    """env-safe VERBATIM copy of scvi.data.read_10x_multiome (uses only mmread + pandas + AnnData; avoids
    the scvi -> chex/jax import cascade absent in scbutterfly-env). The BMMC_d1 dirs are uncompressed v3
    triplets (matrix.mtx/features.tsv/barcodes.tsv) that sc.read_10x_mtx can't parse but scvi (and this
    copy) can. Output matches the scVI baseline exactly: cells x features, var['modality'] from feature
    type, barcodes split on '-' (drops -1)."""
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

# ---- BMMC_d1 input dirs (10x-mtx triplets); mirror scVI 01_run_scvi.py ----
DATA = "/path/to/data/BMMC_d1"
train_10x = f"{DATA}/s1d2_train"          # train data with combined peaks (paired)
test1_10x = f"{DATA}/train_test1_s2d1"    # test1
test2_10x = f"{DATA}/train_test2_s4d1"    # test2
test3_10x = f"{DATA}/train_test3_s1d1"    # test3

SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
OUT  = os.path.join("/path/to/multiomeBench/BMMC_d1/crosssite/scbutterfly", REP)
os.makedirs(OUT, exist_ok=True)
np.random.seed(SEED)   # (scvi.settings.seed removed: scvi is not imported in scbutterfly-env; np+torch seeding suffices)
torch.manual_seed(SEED)
start = timeit.default_timer()


def chrom_key(c):
    m = re.match(r"^chr?(\w+)$", str(c), re.I); t = m.group(1) if m else str(c)
    return (0, int(t)) if t.isdigit() else (1, str(t))


def chrom_of(var):
    # peak var_names are "chr:start-end"; pull the leading chromosome token
    iv = pd.Series(var.index.astype(str), index=var.index)
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


# ---- load train (paired) + 3 tests (mirror scVI input loading EXACTLY) ----
adata_train = read_10x_multiome(train_10x)
adata_train.var_names_make_unique()

adata_test1 = read_10x_multiome(test1_10x)
adata_test1.var_names_make_unique()
adata_test2 = read_10x_multiome(test2_10x)
adata_test2.var_names_make_unique()
adata_test3 = read_10x_multiome(test3_10x)
adata_test3.var_names_make_unique()

adata_test1_rna  = adata_test1[:, adata_test1.var.modality == "Gene Expression"].copy()
adata_test1_atac = adata_test1[:, adata_test1.var.modality == "Peaks"].copy()
adata_test2_rna  = adata_test2[:, adata_test2.var.modality == "Gene Expression"].copy()
adata_test2_atac = adata_test2[:, adata_test2.var.modality == "Peaks"].copy()
adata_test3_rna  = adata_test3[:, adata_test3.var.modality == "Gene Expression"].copy()
adata_test3_atac = adata_test3[:, adata_test3.var.modality == "Peaks"].copy()

# split train the same way (paired -> RNA + ATAC halves) for scButterfly's paired train set
adata_train_rna  = adata_train[:, adata_train.var.modality == "Gene Expression"].copy()
adata_train_atac = adata_train[:, adata_train.var.modality == "Peaks"].copy()
adata_train.obs["data"] = "train"

# rename bc + add data batch (mirror scVI 01_run_scvi.py)
adata_test1_rna.obs.index = adata_test1_rna.obs.index + "_rna1"
adata_test1_rna.obs["data"] = "test1"
adata_test2_rna.obs.index = adata_test2_rna.obs.index + "_rna2"
adata_test2_rna.obs["data"] = "test2"
adata_test3_rna.obs.index = adata_test3_rna.obs.index + "_rna3"
adata_test3_rna.obs["data"] = "test3"

adata_test1_atac.obs.index = adata_test1_atac.obs.index + "_atac1"
adata_test1_atac.obs["data"] = "test1"
adata_test2_atac.obs.index = adata_test2_atac.obs.index + "_atac2"
adata_test2_atac.obs["data"] = "test2"
adata_test3_atac.obs.index = adata_test3_atac.obs.index + "_atac3"
adata_test3_atac.obs["data"] = "test3"

# CROSS-SITE UNPAIRED SPLIT: RNA = s1d1 ONLY (test3 / _rna3), ATAC = s4d1 ONLY (test2 / _atac2).
# s2d1 (test1 / _rna1 / _atac1) is DROPPED. The s1d2 paired TRAIN anchor above is kept unchanged.
adata_rna  = adata_test3_rna.copy()   # s1d1 RNA  (_rna3)
adata_atac = adata_test2_atac.copy()  # s4d1 ATAC (_atac2)
# adata_rna (s1d1) and adata_atac (s4d1) are DIFFERENT cells -> NOT row-aligned (unpaired test).

# ---- feature intersection (train <-> test), per modality (keeps chrom handling robust) ----
g = adata_train_rna.var_names.intersection(adata_rna.var_names)
p = adata_train_atac.var_names.intersection(adata_atac.var_names)
tr_rna,  te_rna  = adata_train_rna[:, g].copy(),  adata_rna[:, g].copy()
tr_atac, te_atac = adata_train_atac[:, p].copy(), adata_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  test_rna {te_rna.n_obs}  test_atac {te_atac.n_obs}")

# ---- build the combined RNA_data / ATAC_data (train rows then test rows), same as RMS ----
RNA_data  = anndata.concat([tr_rna,  te_rna])
ATAC_data = anndata.concat([tr_atac, te_atac])
RNA_data.var  = tr_rna.var
ATAC_data.var = tr_atac.var
ATAC_data = order_atac_by_chrom(ATAC_data)

n_train = tr_rna.n_obs
train_id = list(range(n_train))
# unpaired test: RNA (s1d1) and ATAC (s4d1) have different cell counts -> per-modality test indices
test_id_rna  = list(range(n_train, n_train + te_rna.n_obs))
test_id_atac = list(range(n_train, n_train + te_atac.n_obs))
test_id      = list(range(n_train, n_train + min(te_rna.n_obs, te_atac.n_obs)))  # butterfly.load_data only
rna_bc   = list(te_rna.obs_names)   # <bc>_rna3  (s1d1)
atac_bc  = list(te_atac.obs_names)  # <bc>_atac2 (s4d1)

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


# encode each test modality SEPARATELY: RNA encoder on the test rows -> _rnaN, ATAC encoder -> _atacN
mu_r = embed(butterfly.RNA_data_p,  test_id_rna,  model.RNA_encoder,  "RNA")
mu_a = embed(butterfly.ATAC_data_p, test_id_atac, model.ATAC_encoder, "ATAC")

lat = pd.concat([pd.DataFrame(mu_r, index=rna_bc),
                 pd.DataFrame(mu_a, index=atac_bc)])
lat.to_csv(os.path.join(OUT, "latent.csv"))
stop = timeit.default_timer()
with open(os.path.join(OUT, "scbutterfly_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={lat.shape}  chroms={len(chrom_list)}  time={stop-start:.1f}s")
