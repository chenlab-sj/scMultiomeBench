# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
scButterfly on the Parse cross-platform dataset -> latent.csv for benchmark_metrics.py.

CROSS-PLATFORM SETUP (differs from the pbmc3k 10x version):
  TRAIN  = paired pbmc3k multiome (RNA+ATAC, SAME cells)  -> trains the translator (paired).
  TEST RNA  = Parse whole-cell RNA  (independent cells, gene-symbol features).
  TEST ATAC = pbmc3k ATAC peaks     (independent cells, DIFFERENT barcodes from TEST RNA).
The TEST RNA and TEST ATAC are UNPAIRED (different platforms, different cells, different counts),
exactly like the analogous category-3 methods (pbmc_parse/scVI, pbmc_parse/scglue_paired).

scButterfly is a dual-aligned VAE TRANSLATION method; its public API returns translated
profiles, not a latent. Path A (verbatim from the pbmc3k version): train on the paired TRAIN
multiome, then pull the SHARED translator latent mean (mu) for the test cells, embedding each
modality SEPARATELY (so no test pairing is ever used -> a fair unpaired co-embedding):
    R2 = model.RNA_encoder(x);   _,_, mu_r,_ = model.translator.test_model(R2, 'RNA')
    A2 = model.ATAC_encoder(x);  _,_, mu_a,_ = model.translator.test_model(A2, 'ATAC')

Grounded against the source:
  * data_preprocessing() -> butterfly.RNA_data_p / butterfly.ATAC_data_p (AnnData): RNA = 3000
    HVG, normalized+log; ATAC = binarized/TF-IDF/[0,1]. These .X are the encoder inputs (NOT
    raw counts). Cell rows are preserved (preprocessing filters FEATURES, not cells).
  * RNA_data / ATAC_data are preprocessed INDEPENDENTLY per modality, and load_data stores
    train/test ids per modality (train_id_r/_a, test_id_r/_a) -> the two modalities are NOT
    required to share cells. We therefore concat TRAIN+TEST PER MODALITY (RNA: train_rna +
    parse_test_rna ; ATAC: train_atac + pbmc3k_test_atac) even though the two test blocks have
    different cells/counts. Only the TRAIN block is paired (and only the TRAIN block is used for
    training, since train_id = range(n_train)).
  * construct_model(chrom_list) needs peak counts per chromosome of the PROCESSED ATAC
    (Split_Chrom_Encoder). We sort peaks by chromosome before load_data and compute chrom_list
    from ATAC_data_p afterward.
  * load_data() pairs RNA_data[train_id_r] with ATAC_data[train_id_a] by ROW POSITION during
    training; train_id is the same range for both modalities -> the TRAIN block MUST be paired
    and row-aligned. It is (it comes from the same multiome h5, split by feature_type, never
    reordered on the cell axis). test_id is stored but NOT used by our Path-A extraction (we
    embed each modality with its own explicit row indices), so the unequal test counts are fine.

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
INP  = f"{ROOT}/pbmc_parse/input_azimuth"
TRAIN_H5    = f"{INP}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_train.h5"  # paired multiome
TEST_RNA_H5 = f"{INP}/pbmc_parse_d1.h5"        # Parse whole-cell RNA (independent cells)
TEST_ATAC_H5 = f"{INP}/pbmc3k_test_atac.h5"    # pbmc3k ATAC peaks (independent cells)
# reproducibility reps: REP=""/rep2/rep3 + SEED -> pbmc_parse/scButterfly/{,rep2,rep3}/latent.csv.
# scButterfly hard-codes setup_seed(19193) internally, so plain reps are byte-identical; FORCE_SEED=1
# overrides it (re-seed before construct_model for weight init + pass seed to train_model) so the
# external SEED genuinely varies the run -> honest seed-robustness (mirrors the pbmc3k scButterfly setup).
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
FORCE_SEED = os.environ.get("FORCE_SEED", "0") == "1"
OUT  = os.path.join(f"{ROOT}/pbmc_parse/scButterfly", REP)
os.makedirs(OUT, exist_ok=True)
np.random.seed(SEED)
torch.manual_seed(SEED)
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

# ---- build per-modality data: TRAIN (paired) + TEST (independent cells) ----
# RNA = train_rna  (+) Parse test RNA ;  ATAC = train_atac (+) pbmc3k test ATAC.
# The two TEST blocks are different cells/counts -> RNA_data and ATAC_data have DIFFERENT n_obs.
tr_rna, tr_atac = load_split(TRAIN_H5)

te_rna  = sc.read_10x_h5(TEST_RNA_H5);                te_rna.var_names_make_unique()
te_atac = sc.read_10x_h5(TEST_ATAC_H5, gex_only=False)
# keep only Peaks if the test ATAC h5 carries a feature_types column
if "feature_types" in te_atac.var and (te_atac.var["feature_types"] == "Peaks").any():
    te_atac = te_atac[:, te_atac.var["feature_types"] == "Peaks"].copy()
te_atac.var_names_make_unique()

# concat per modality on the cell axis (inner join on shared features, mirrors pbmc3k version's sc.concat)
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
# test_id passed to load_data is stored only (not used by our Path-A extraction); use the
# shorter range so it indexes validly into BOTH modalities.
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
