# NOTE: paths below are placeholders. See config/config.yaml and the README.
"""
scButterfly on BRCA HT137B1-S1H7 (GENUINELY UNPAIRED RNA + ATAC) -> latent.csv for the benchmark.

HT137 is same-donor INDEPENDENT sequencing: the test scRNA (15293 cells) and test scATAC (5127
cells) are DIFFERENT cells (only ~5 barcode strings coincide by chance). So -- unlike the HT243
run, where test RNA/ATAC were a modality split of the SAME paired cells -- we must NOT intersect
the test barcodes. This uses the cross-platform Parse design (pbmc_parse/scButterfly): concat
TRAIN+TEST PER MODALITY and embed each test modality SEPARATELY by its own row indices, so no test
pairing is ever assumed. Only the file loading is BRCA (train 10x h5 + RNA h5ad + commonpeak h5).

Path A (verbatim from the pbmc3k/Parse version): train the dual-aligned VAE on the paired TRAIN
multiome, then pull the SHARED translator latent mean (mu) for the test cells:
    R2 = model.RNA_encoder(x);   _,_, mu_r,_ = model.translator.test_model(R2, 'RNA')
    A2 = model.ATAC_encoder(x);  _,_, mu_a,_ = model.translator.test_model(A2, 'ATAC')

Inputs (BRCA, raw counts; barcodes are RAW -> we suffix _rna/_atac only at output, matching label.csv):
  * TRAIN     = train_forHT137B1-S1H7.h5   (combined multiome 10x h5 -> split; paired, row-aligned)
  * TEST RNA  = HT137B1-S1H7_rna.h5ad      (independent test RNA, the file scglue/scDART use)
  * TEST ATAC = HT137B1-S1H7_commonpeak.h5 (independent test ATAC, SAME peak set as the train)

VERIFY on first run: (1) train h5 splits into Gene Expression + Peaks; (2) test h5ad/h5 are raw
counts; (3) preprocessing keeps all cells (RNA_data_p.n_obs == RNA_data.n_obs, same for ATAC);
(4) sum(chrom_list)==ATAC_data_p.n_vars; (5) latent.csv ~ 15293 _rna + 5127 _atac (matches label.csv).
"""
import os, re, timeit
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scButterfly.butterfly import Butterfly

ROOT = "/path/to/multiomeBench"
BRCA = f"{ROOT}/BRCA/HT137B1-S1H7"
TRAIN_H5  = f"{BRCA}/train_forHT137B1-S1H7.h5"
TEST_RNA  = f"{BRCA}/HT137B1-S1H7_rna.h5ad"        # independent test RNA (raw counts, raw barcodes)
# common peaks = the SAME peak set as the train (chr:start-end). Read the MTX DIR, not the .h5:
# HT137B1-S1H7_commonpeak.h5 was written incompletely in 2024 (no /matrix/data,indices,indptr,shape)
# -> scanpy read_10x_h5 fails. The commonpeak/ MTX dir is the complete data (612783 peaks x 5127 cells),
# the same dir scVI reads. (HT137B1-S1H7_atac.h5ad is a DIFFERENT per-sample peak set -> not used here.)
TEST_ATAC = f"{BRCA}/HT137B1-S1H7_commonpeak"
# reproducibility hook (default = main, seed 420): REP=""/rep2/rep3 + SEED 420/0/40
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
FORCE_SEED = os.environ.get("FORCE_SEED", "0") == "1"  # override scButterfly internal seed 19193 (else rep2/3 byte-identical)
OUT  = os.path.join(f"{BRCA}/scbutterfly", REP)
os.makedirs(OUT, exist_ok=True)
np.random.seed(SEED)
torch.manual_seed(SEED)
start = timeit.default_timer()


def load_split(h5):
    a = sc.read_10x_h5(h5, gex_only=False); a.var_names_make_unique()
    return (a[:, a.var["feature_types"] == "Gene Expression"].copy(),
            a[:, a.var["feature_types"] == "Peaks"].copy())


def raw_counts(ad):
    """Return an AnnData whose .X is RAW counts (scButterfly transforms internally)."""
    if "counts" in ad.layers:
        ad.X = ad.layers["counts"]
    elif ad.raw is not None and ad.raw.shape[1] == ad.shape[1]:
        ad.X = ad.raw.X
    x = ad.X[:50].toarray() if hasattr(ad.X, "toarray") else np.asarray(ad.X[:50])
    if not np.allclose(x, np.round(x)):
        print(f"  WARNING: {ad.shape} .X not integer counts -> scButterfly expects raw counts")
    return ad


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
    """consecutive per-chromosome peak counts of an already chrom-ordered var."""
    chrom = var["chrom"].values if "chrom" in var else chrom_of(var)
    sizes, run, prev = [], 0, None
    for c in chrom:
        if prev is not None and c != prev: sizes.append(run); run = 0
        run += 1; prev = c
    if run: sizes.append(run)
    return [int(s) for s in sizes]


# ---- build per-modality data: TRAIN (paired) + TEST (independent cells) ----
# RNA = train_rna (+) HT137 test RNA ;  ATAC = train_atac (+) HT137 test ATAC.
# The two TEST blocks are DIFFERENT cells/counts -> RNA_data and ATAC_data have DIFFERENT n_obs.
tr_rna, tr_atac = load_split(TRAIN_H5)

te_rna  = sc.read_h5ad(TEST_RNA); te_rna.var_names_make_unique(); te_rna = raw_counts(te_rna)
# commonpeak MTX dir (prefer read_10x_mtx; manual scipy fallback for the 6-column features.tsv)
try:
    te_atac = sc.read_10x_mtx(TEST_ATAC, gex_only=False)
except Exception as e:
    print(f"  read_10x_mtx fallback ({type(e).__name__}: {e}) -> manual scipy load")
    import scipy.io, scipy.sparse as sp
    M = scipy.io.mmread(os.path.join(TEST_ATAC, "matrix.mtx")).T.tocsr()      # features x cells -> cells x features
    bc = pd.read_csv(os.path.join(TEST_ATAC, "barcodes.tsv"), header=None)[0].astype(str).values
    ft = pd.read_csv(os.path.join(TEST_ATAC, "features.tsv"), header=None, sep="\t")
    te_atac = sc.AnnData(X=sp.csr_matrix(M), obs=pd.DataFrame(index=bc),
                         var=pd.DataFrame({"feature_types": "Peaks"}, index=ft[0].astype(str).values))
te_atac.var_names_make_unique()
if "feature_types" in te_atac.var and (te_atac.var["feature_types"] == "Peaks").any():
    te_atac = te_atac[:, te_atac.var["feature_types"] == "Peaks"].copy()
te_atac = raw_counts(te_atac)

# feature-align train<->test per modality (shared genes / shared peaks)
g = tr_rna.var_names.intersection(te_rna.var_names)
p = tr_atac.var_names.intersection(te_atac.var_names)
tr_rna, te_rna   = tr_rna[:, g].copy(),  te_rna[:, g].copy()
tr_atac, te_atac = tr_atac[:, p].copy(), te_atac[:, p].copy()
print(f"shared genes {len(g)}  shared peaks {len(p)}  | train {tr_rna.n_obs}  "
      f"test {te_rna.n_obs}(rna)+{te_atac.n_obs}(atac)")

# concat per modality on the cell axis (inner join on shared features)
RNA_data  = sc.concat([tr_rna,  te_rna])
ATAC_data = sc.concat([tr_atac, te_atac])
ATAC_data = order_atac_by_chrom(ATAC_data)        # peaks grouped by chromosome

n_train      = tr_rna.n_obs                        # train cells are paired (same multiome)
n_test_rna   = te_rna.n_obs
n_test_atac  = te_atac.n_obs
train_id = list(range(n_train))
test_id_rna  = list(range(n_train, n_train + n_test_rna))
test_id_atac = list(range(n_train, n_train + n_test_atac))
test_bc_rna  = list(RNA_data.obs_names[n_train:])
test_bc_atac = list(ATAC_data.obs_names[n_train:])
# test_id passed to load_data is stored only (not used by Path-A extraction); use the shorter
# range so it indexes validly into BOTH modalities.
test_id = list(range(n_train, n_train + min(n_test_rna, n_test_atac)))

# ---- scButterfly (aug_type=None -> no scvi) ----
butterfly = Butterfly()
butterfly.load_data(RNA_data, ATAC_data, train_id, test_id, validation_id=None)
butterfly.data_preprocessing()                    # -> RNA_data_p / ATAC_data_p
assert butterfly.RNA_data_p.n_obs == RNA_data.n_obs, "preprocessing dropped RNA cells!"
assert butterfly.ATAC_data_p.n_obs == ATAC_data.n_obs, "preprocessing dropped ATAC cells!"
chrom_list = chrom_counts(butterfly.ATAC_data_p.var)
assert sum(chrom_list) == butterfly.ATAC_data_p.n_vars, (sum(chrom_list), butterfly.ATAC_data_p.n_vars)
if FORCE_SEED:
    from scButterfly.split_datasets import setup_seed
    setup_seed(SEED)                              # override the internal 19193 for weight init
butterfly.construct_model(chrom_list)

# Reuse the trained model if it's already saved (./model/*.pt) -> skip the retrain.
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


mu_r = embed(butterfly.RNA_data_p,  test_id_rna,  model.RNA_encoder,  "RNA")
mu_a = embed(butterfly.ATAC_data_p, test_id_atac, model.ATAC_encoder, "ATAC")

lat = pd.concat([pd.DataFrame(mu_r, index=[f"{b}_rna"  for b in test_bc_rna]),
                 pd.DataFrame(mu_a, index=[f"{b}_atac" for b in test_bc_atac])])
lat.to_csv(os.path.join(OUT, "latent.csv"))
stop = timeit.default_timer()
with open(os.path.join(OUT, "scbutterfly_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={lat.shape}  chroms={len(chrom_list)}  "
      f"n_rna={len(test_bc_rna)} n_atac={len(test_bc_atac)}  time={stop-start:.1f}s")
