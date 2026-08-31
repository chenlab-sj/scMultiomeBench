# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MIRA quick-run / smoke test on pbmc3k.

GOAL of THIS script: prove the env works end-to-end on a GPU node and that MIRA
trains topic models on your pbmc3k matrices + emits a barcode-indexed latent.csv
(the same I/O contract as 01_run_scvi.py / 01_run_scglue.py -> benchmark_metrics.py).

It runs the PAIRED MIRA workflow on the TRAIN multiome (the path MIRA officially
supports: mira.utils.make_joint_representation needs the SAME cells in both
modalities). That is a legitimate smoke test of "does MIRA run on my data".

It is NOT the unpaired benchmark run -- that is run_mira_unpaired.py (option A:
train-anchored CCA bridge, the script the benchmark actually uses). See the big
NOTE at the bottom for the design it implements.

API note: models are built with mira.topics.make_model(...) — the real constructor in
mira 2.1.1 (ExpressionTopicModel / AccessibilityTopicModel are doc-only faux classes).
"""

import os
import timeit
import numpy as np
import pandas as pd
import scanpy as sc
import mira

# ------------------------------------------------------------------ I/O --
# Cluster paths (LSF mounts these as /path/to/project/...).
ROOT = "/path/to/multiomeBench"
DATA = f"{ROOT}/pbmc/pbmc3k/Data"
TRAIN_H5 = f"{DATA}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_train.h5"  # paired multiome (GEX+Peaks)
TEST_H5  = f"{DATA}/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"   # test cells (GEX+Peaks); for the unpaired run
OUT_DIR  = f"{ROOT}/pbmc/pbmc3k/scripts/MIRA"

SMOKE = True      # True: subsample cells + few epochs so the first job returns fast
SEED  = 420
mira.utils.pretty_print = lambda *a, **k: None  # quiet; remove if you want MIRA's banners

os.makedirs(OUT_DIR, exist_ok=True)
np.random.seed(SEED)
start = timeit.default_timer()

# --------------------------------------------------------------- load ----
# Read the combined multiome .h5 (cleaner than the uncompressed mtx dirs, which
# sc.read_10x_mtx mis-detects as v3 and hunts for matrix.mtx.gz). Split by feature_types.
def load_multiome(h5_path):
    a = sc.read_10x_h5(h5_path, gex_only=False)
    a.var_names_make_unique()
    rna_  = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    atac_ = a[:, a.var["feature_types"] == "Peaks"].copy()
    return rna_, atac_

rna, atac = load_multiome(TRAIN_H5)   # paired train multiome for the smoke run

if SMOKE:
    idx = np.random.choice(rna.n_obs, size=min(800, rna.n_obs), replace=False)
    rna, atac = rna[idx].copy(), atac[idx].copy()
print(f"[smoke={SMOKE}] rna {rna.shape}  atac {atac.shape}")

# MIRA trains on raw counts kept in a layer; .X can be processed.
for ad_ in (rna, atac):
    ad_.layers["counts"] = ad_.X.copy()

# RNA: highly-variable genes drive the model's exogenous features
sc.pp.normalize_total(rna, target_sum=1e4); sc.pp.log1p(rna)
sc.pp.highly_variable_genes(rna, min_disp=0.5)
rna.X = rna.layers["counts"].copy()  # restore counts for the topic model

# ATAC: every peak is a feature (binary/count accessibility)
atac.var["highly_variable"] = True

# ----------------------------------------------------- topic models ------
# DEFAULTS POLICY: the benchmark runs every method at its OWN defaults. make_model
# auto-scales the architecture to dataset size via recommend_parameters() (it picked
# 14 topics / 24 epochs here), and we keep that. We pass only structural/data args +
# seed; we do NOT hand-tune num_topics, epochs, learning rates, or box_cox.
# EXCEPTION: MIRA's accessibility SkipEncoder requires num_layers >= 3, but on small
# (e.g. SMOKE) subsets recommend_parameters() drops it to 2 -> AssertionError. We pin
# the ATAC model to num_layers=3, which IS MIRA's documented default (so still "at
# defaults", just guarded against the small-data down-scaling). User-passed params win
# over recommend_parameters (parameter_recommendations.update(model_parameters)).

# NB: ExpressionTopicModel / AccessibilityTopicModel are doc-only "faux" classes in
# mira 2.1.1 — instantiating them raises NotImplementedError. The real constructor is
# mira.topics.make_model(n_samples, n_features, feature_type=..., ...).
rna_model = mira.topics.make_model(
    rna.n_obs, rna.n_vars,
    feature_type="expression",
    highly_variable_key="highly_variable",   # exogenous features = HVGs
    counts_layer="counts",
    seed=SEED,
)
atac_model = mira.topics.make_model(
    atac.n_obs, atac.n_vars,
    feature_type="accessibility",
    counts_layer="counts",
    num_layers=3,   # MIRA's default; the accessibility SkipEncoder requires >=3.
    seed=SEED,
)

for model, ad_ in [(rna_model, rna), (atac_model, atac)]:
    model.get_learning_rate_bounds(ad_)  # LR range test; auto-sets min/max LR for fit()
    model.fit(ad_)                        # trains at MIRA-default epochs/topics
    model.predict(ad_)                    # -> obsm["X_topic_compositions"]
    model.get_umap_features(ad_)          # default box_cox -> obsm["X_umap_features"]

# -------------------------------------- PAIRED joint rep (smoke only) ----
# Supported only because train RNA & ATAC are the SAME cells.
rna, atac = mira.utils.make_joint_representation(rna, atac)  # -> obsm["X_joint_umap_features"]
joint = pd.DataFrame(rna.obsm["X_joint_umap_features"], index=rna.obs_names)

latent_csv = os.path.join(OUT_DIR, "latent_smoke.csv")
joint.to_csv(latent_csv)
stop = timeit.default_timer()
with open(os.path.join(OUT_DIR, "mira_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote {latent_csv}  shape={joint.shape}  time={stop-start:.1f}s")

# =======================================================================
# NOTE — the UNPAIRED benchmark run is implemented in run_mira_unpaired.py:
#
# make_joint_representation() above needs PAIRED cells, so it cannot be used
# on the unpaired test set (separate `<bc>_rna` and `<bc>_atac` cells).
# MIRA out of the box is a paired method. To put it in your "category-1"
# unpaired panel you must add a cross-modal alignment step. Options:
#
#   A) Anchor on the TRAIN multiome: learn both topic models on paired train,
#      use it to define a shared space, project test RNA & test ATAC in, then
#      align topic spaces (e.g. CCA on paired-train umap features). Most defensible.
#      *** This is what run_mira_unpaired.py implements. ***
#   B) Translate ATAC->gene-activity, run a single expression topic model
#      (make_model feature_type='expression') over RNA + gene-activity together.
#      Loses MIRA's peak modeling (its point).
#   X) Paired "cheat" upper bound: run make_joint_representation on the test
#      cells using their known pairing -> NOT comparable to the other 18
#      methods (uses pairing they don't get). Useful only as a sanity ceiling.
# =======================================================================
