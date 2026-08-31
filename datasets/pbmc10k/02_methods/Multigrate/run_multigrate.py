# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
Multigrate on pbmc10k (unpaired RNA + ATAC) -> latent.csv for benchmark_metrics.py.

Multigrate is a scvi-tools VAE (like MultiVI in 01_run_scvi.py). It builds a joint latent
across modalities. We anchor on the paired TRAIN multiome and add the unpaired test
RNA / test ATAC, then read off the latent for the test cells (<bc>_rna / <bc>_atac).

API (multigrate 1.0.1, verified from the paired-integration tutorial):
  adata = mtg.data.organize_multimodal_anndatas(adatas=[[rna],[atac]], layers=[['counts'],['counts']])
  mtg.model.MultiVAE.setup_anndata(adata, rna_indices_end=<n_rna_genes>, categorical_covariate_keys=['group'])
  vae = mtg.model.MultiVAE(adata, losses=['nb','nb']); vae.train()
  vae.get_model_output()              # -> adata.obsm['X_multigrate']

NOTE (verify on first run): the exact nested structure of organize_multimodal_anndatas for
the PAIRED-anchor + UNPAIRED-query (mosaic) case is the one thing to confirm against your
installed version -- analogous to MultiVI's organize_multiome_anndatas(paired, rna, atac).
If pairing-by-obs_names isn't picked up, switch to the pure-diagonal form [[rna_test],[atac_test]].
"""
import os
import timeit
import numpy as np
import pandas as pd
import scanpy as sc
import multigrate as mtg

ROOT = "/path/to/multiomeBench"
DATA = "/path/to/data/pbmc10k/pbmc10k_data_design0719"
TRAIN_H5 = f"{DATA}/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_train.h5"  # paired multiome
TEST_H5  = f"{DATA}/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5"   # test cells
OUT  = f"{ROOT}/pbmc/pbmc10k/scripts/multigrate"
SEED = 420
os.makedirs(OUT, exist_ok=True)
start = timeit.default_timer()

def load_split(h5, suffix=None):
    a = sc.read_10x_h5(h5, gex_only=False); a.var_names_make_unique()
    rna  = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    atac = a[:, a.var["feature_types"] == "Peaks"].copy()
    if suffix:                                   # tag unpaired test cells
        rna.obs_names  = [f"{b}_rna"  for b in rna.obs_names]
        atac.obs_names = [f"{b}_atac" for b in atac.obs_names]
    return rna, atac

tr_rna, tr_atac = load_split(TRAIN_H5)              # paired anchor (tr_rna & tr_atac share obs_names)
te_rna, te_atac = load_split(TEST_H5, suffix=True)  # unpaired query (_rna / _atac)

for m in (tr_rna, tr_atac, te_rna, te_atac):        # raw counts for the nb loss
    m.layers["counts"] = m.X.copy()

# Mosaic layout = adatas[modality][group]; None -> zero-filled placeholder. organize
# requires non-None modalities WITHIN a group to share obs_names (i.e. be paired):
#   group 0 = paired TRAIN (tr_rna & tr_atac, same cells)
#   group 1 = test RNA only ; group 2 = test ATAC only
# organize creates the .obs['group'] column itself (don't set it manually).
adata = mtg.data.organize_multimodal_anndatas(
    adatas=[[tr_rna,  te_rna, None],      # RNA modality across the 3 groups
            [tr_atac, None,   te_atac]],  # ATAC modality
    layers=[["counts", "counts", None],
            ["counts", None,     "counts"]],
)
mtg.model.MultiVAE.setup_anndata(adata, rna_indices_end=tr_rna.n_vars,
                                 categorical_covariate_keys=["group"])
vae = mtg.model.MultiVAE(adata, losses=["nb", "nb"])   # ATAC loss may need tuning
vae.train()
vae.get_model_output()                                  # -> adata.obsm['X_multigrate']

# ---- latent.csv for the TEST cells only (the eval set) ----
lat = pd.DataFrame(adata.obsm["X_multigrate"], index=adata.obs_names)
test_bc = list(te_rna.obs_names) + list(te_atac.obs_names)
lat = lat.loc[lat.index.isin(test_bc)]
lat.to_csv(os.path.join(OUT, "latent.csv"))
stop = timeit.default_timer()
with open(os.path.join(OUT, "multigrate_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={lat.shape}  time={stop-start:.1f}s")
