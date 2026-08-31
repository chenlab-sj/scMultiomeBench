# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIRA UNPAIRED benchmark run on pbmc3k  (Option A from run_mira.py's NOTE).

Why this exists: make_joint_representation() needs PAIRED cells, so it cannot be
used on the unpaired test set. MIRA is natively a paired method; to place it in
the unpaired panel we bridge through the paired TRAIN multiome. This makes MIRA a
CATEGORY-3 method (it uses the multiome reference), NOT category-1 as a reviewer
suggested -- report it that way.

Pipeline:
  1. Fit RNA + ATAC topic models on the PAIRED TRAIN multiome (same cells).
  2. Use those SAME models to embed train RNA, train ATAC, test RNA, test ATAC
     -> obsm["X_umap_features"] (the features MIRA itself uses for neighbors).
  3. Fit CCA on the PAIRED TRAIN umap-features (RNA<->ATAC of the same cells)
     -> a shared canonical space.
  4. Project the UNPAIRED test RNA and test ATAC into that space SEPARATELY,
     then concatenate -> one latent.csv indexed by <bc>_rna / <bc>_atac.

NOTE on stability: MIRA's get_learning_rate_bounds() can pick a max LR that
overflows the gradient (-> ModelParamError). fit_stable() backs the max LR off
until training converges, instead of failing the whole job.
"""

import os
import timeit
import numpy as np
import pandas as pd
import scanpy as sc
import mira
from sklearn.cross_decomposition import CCA

try:
    from mira.topic_model.base import ModelParamError
except Exception:
    ModelParamError = Exception

# ------------------------------------------------------------------ I/O --
ROOT = "/path/to/multiomeBench"
DATA = f"{ROOT}/pbmc/pbmc3k/Data"
TRAIN_H5 = f"{DATA}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_train.h5"  # paired multiome reference
TEST_H5  = f"{DATA}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"   # unpaired test cells
OUT_DIR  = f"{ROOT}/pbmc/pbmc3k/scripts/MIRA"

SMOKE = False     # full run: no subsample (full train ~4k + full test ~1642+1642). True = fast smoke.
SEED  = 420
N_CCA = 15        # canonical dims for the cross-modal bridge (clipped to <= min topic dims)
mira.utils.pretty_print = lambda *a, **k: None

os.makedirs(OUT_DIR, exist_ok=True)
np.random.seed(SEED)
start = timeit.default_timer()

# --------------------------------------------------------------- load ----
def load_multiome(h5_path):
    a = sc.read_10x_h5(h5_path, gex_only=False)
    a.var_names_make_unique()
    rna_  = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    atac_ = a[:, a.var["feature_types"] == "Peaks"].copy()
    return rna_, atac_

rna_tr, atac_tr = load_multiome(TRAIN_H5)   # paired reference
rna_te, atac_te = load_multiome(TEST_H5)    # unpaired test

if SMOKE:
    itr = np.random.choice(rna_tr.n_obs, size=min(800, rna_tr.n_obs), replace=False)
    rna_tr, atac_tr = rna_tr[itr].copy(), atac_tr[itr].copy()
    ite = np.random.choice(rna_te.n_obs, size=min(400, rna_te.n_obs), replace=False)
    rna_te, atac_te = rna_te[ite].copy(), atac_te[ite].copy()

# unpaired test identity: <bc>_rna / <bc>_atac (the other methods' convention)
rna_te.obs_names  = [f"{bc}_rna"  for bc in rna_te.obs_names]
atac_te.obs_names = [f"{bc}_atac" for bc in atac_te.obs_names]

# test must share the train feature space (same model predicts both); keep train order
def align_features(ad_ref, ad_new):
    shared = ad_ref.var_names.intersection(ad_new.var_names)
    return ad_ref[:, shared].copy(), ad_new[:, shared].copy()

rna_tr,  rna_te  = align_features(rna_tr,  rna_te)
atac_tr, atac_te = align_features(atac_tr, atac_te)

# drop near-empty features (defined on TRAIN), then re-subset TEST to match.
# Degenerate all-zero / ultra-sparse features are a common cause of gradient blow-up.
sc.pp.filter_genes(rna_tr,  min_cells=3)
sc.pp.filter_genes(atac_tr, min_cells=3)
rna_te  = rna_te[:,  rna_tr.var_names].copy()
atac_te = atac_te[:, atac_tr.var_names].copy()
print(f"[smoke={SMOKE}] train rna {rna_tr.shape} atac {atac_tr.shape} | "
      f"test rna {rna_te.shape} atac {atac_te.shape}")

# ---------------------------------------------------- preprocessing ------
for ad_ in (rna_tr, rna_te, atac_tr, atac_te):
    ad_.layers["counts"] = ad_.X.copy()

# RNA HVGs are defined on TRAIN and carried to TEST (same model, same features)
sc.pp.normalize_total(rna_tr, target_sum=1e4); sc.pp.log1p(rna_tr)
sc.pp.highly_variable_genes(rna_tr, min_disp=0.5)
rna_tr.X = rna_tr.layers["counts"].copy()
rna_te.var["highly_variable"] = rna_tr.var["highly_variable"].values  # identical var order
rna_te.X = rna_te.layers["counts"].copy()

for ad_ in (atac_tr, atac_te):
    ad_.var["highly_variable"] = True

# ----------------------------------------------------- topic models ------
# DEFAULTS POLICY unchanged: MIRA defaults via make_model + recommend_parameters;
# only num_layers=3 pinned for the accessibility SkipEncoder (its documented default).
def make_rna_model():
    return mira.topics.make_model(
        rna_tr.n_obs, rna_tr.n_vars, feature_type="expression",
        highly_variable_key="highly_variable", counts_layer="counts", seed=SEED)

def make_atac_model():
    return mira.topics.make_model(
        atac_tr.n_obs, atac_tr.n_vars, feature_type="accessibility",
        counts_layer="counts", num_layers=3, seed=SEED)

def fit_stable(make_fn, ad, name):
    """Fit a MIRA topic model, backing off the max LR if the gradient overflows."""
    model = make_fn()
    bounds = model.get_learning_rate_bounds(ad)
    try:
        lo, hi = bounds
    except (TypeError, ValueError):
        lo = getattr(model, "min_learning_rate", 1e-5)
        hi = getattr(model, "max_learning_rate", 0.1)
    print(f"[{name}] auto LR bounds = ({lo:.2e}, {hi:.2e})")
    for backoff in (3.0, 6.0, 12.0, 30.0, 60.0):
        max_lr = hi / backoff
        try:
            if hasattr(model, "set_learning_rates"):
                model.set_learning_rates(lo, max_lr)
            else:
                model.min_learning_rate, model.max_learning_rate = lo, max_lr
            model.fit(ad)
            print(f"[{name}] converged at max_lr = {max_lr:.2e}")
            return model
        except ModelParamError:
            print(f"[{name}] gradient overflow at max_lr={max_lr:.2e}; backing off")
            model = make_fn()  # fresh params for the retry
    raise RuntimeError(f"{name} topic model did not converge; try a different SEED")

rna_model  = fit_stable(make_rna_model,  rna_tr, "rna")
atac_model = fit_stable(make_atac_model, atac_tr, "atac")

def umap_feats(model, ad_):
    model.predict(ad_)            # -> obsm["X_topic_compositions"]
    model.get_umap_features(ad_)  # -> obsm["X_umap_features"]
    return np.asarray(ad_.obsm["X_umap_features"])

R_tr = umap_feats(rna_model,  rna_tr)    # paired train RNA features
A_tr = umap_feats(atac_model, atac_tr)   # paired train ATAC features (same cells as R_tr)
R_te = umap_feats(rna_model,  rna_te)    # unpaired test RNA
A_te = umap_feats(atac_model, atac_te)   # unpaired test ATAC

# ----------------------------------------- CCA bridge (paired train) -----
# Standardize with TRAIN stats ourselves (scale=False) so the test projection uses
# the SAME transform and we depend only on the stable public rotations_ attributes.
def zfit(X):
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    return mu, sd

xmu, xsd = zfit(R_tr)
ymu, ysd = zfit(A_tr)
Rz_tr, Az_tr = (R_tr - xmu) / xsd, (A_tr - ymu) / ysd

k = int(min(N_CCA, R_tr.shape[1], A_tr.shape[1]))
cca = CCA(n_components=k, scale=False, max_iter=2000)
cca.fit(Rz_tr, Az_tr)

diff = np.max(np.abs(cca.transform(Rz_tr) - Rz_tr @ cca.x_rotations_))
print(f"[CCA] k={k}  manual-vs-transform max|diff|={diff:.2e} (should be ~0)")

rna_proj  = ((R_te - xmu) / xsd) @ cca.x_rotations_   # test RNA  -> shared space
atac_proj = ((A_te - ymu) / ysd) @ cca.y_rotations_   # test ATAC -> shared space

joint = pd.concat([
    pd.DataFrame(rna_proj,  index=rna_te.obs_names),
    pd.DataFrame(atac_proj, index=atac_te.obs_names),
])

out = os.path.join(OUT_DIR, "latent.csv")
joint.to_csv(out)
stop = timeit.default_timer()
with open(os.path.join(OUT_DIR, "mira_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote {out}  shape={joint.shape}  "
      f"(rna={rna_proj.shape[0]}, atac={atac_proj.shape[0]})  time={stop-start:.1f}s")
