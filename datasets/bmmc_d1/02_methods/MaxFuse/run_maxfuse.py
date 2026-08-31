# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MaxFuse on BMMC_d1 (unpaired RNA + ATAC, 3 test batches s2d1/s4d1/s1d1 concatenated) -> latent.csv.

Same Fusor workflow as the RMS Mast607A / pbmc3k MaxFuse runs (tutorial defaults), but the
Signac gene-activity prep is REUSED from the BMMC 00_file_prep.R outputs instead of being re-derived
in R. MaxFuse has NO dataset-batch correction parameter, so the 3 test batches are simply taken as
one concatenated RNA cell set vs one concatenated ATAC cell set (single run, no batch_key).

INPUTS (already exist in the BMMC_d1 folder, built by 00_file_prep.R; see 00_scjoint_s1.R for the same
recipe):
  test_rna.h5ad        cells x genes  RNA counts            obs = <bc>_rna1/_rna2/_rna3
  test_atac_gene.h5ad  cells x genes  ATAC GeneActivity     obs = <bc>_atac1/_atac2/_atac3
  atac_embed.h5ad      cells x 50     integrated ATAC LSI    obs = <bc>_atac1/_atac2/_atac3
The two count h5ads carry gene symbols as var_names, so the MaxFuse "shared" feature space is the
intersection of RNA genes and ATAC-activity genes (mirrors 00_file_prep.R / scjoint common.genes).

MaxFuse Fusor arrays (same roles as RMS 01_run_maxfuse.py):
  shared_arr1 = RNA   lognorm on shared genes      (rows = <bc>_rna*)
  shared_arr2 = ATAC  GeneActivity lognorm, same genes (rows = <bc>_atac*)
  active_arr1 = RNA   scaled HVG                    (rows = <bc>_rna*)
  active_arr2 = ATAC  integrated LSI (atac_embed)   (rows = <bc>_atac*)

OUTPUT: latent.csv (joint CCA embedding) indexed by the suffixed barcodes, matching celltype.csv.

Run in maxfuse-env (CPU, python3.8, numpy<1.24). See 01_run_maxfuse_sub.sh.
Reproducibility triplicate via env SEED/REP -> maxfuse/{,rep2,rep3}/latent.csv (prep is deterministic).
"""
import os
import timeit
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import maxfuse as mf

# Mirror BMMC scVI's data location (cluster paths). out_dir = the maxfuse folder.
DATA = "/path/to/data/BMMC_d1"
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
OUT  = os.path.join("/path/to/multiomeBench/BMMC_d1/scripts/maxfuse", REP)
os.makedirs(OUT, exist_ok=True)
np.random.seed(SEED)
start = timeit.default_timer()


def to_dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X)


# ---- load the pre-built BMMC gene-activity / RNA / LSI h5ads ----
rna       = sc.read_h5ad(f"{DATA}/test_rna.h5ad")        # RNA counts,  <bc>_rna*
atac_gene = sc.read_h5ad(f"{DATA}/test_atac_gene.h5ad")  # ATAC GeneActivity counts, <bc>_atac*
atac_lsi  = sc.read_h5ad(f"{DATA}/atac_embed.h5ad")      # integrated ATAC LSI (50d), <bc>_atac*

rna.var_names_make_unique()
atac_gene.var_names_make_unique()

# align the ATAC LSI embedding to the ATAC gene-activity cell order (both are <bc>_atac*)
atac_lsi = atac_lsi[atac_gene.obs_names].copy()

rna_bc  = rna.obs_names
atac_bc = atac_gene.obs_names

# ---- shared feature space = genes common to RNA and ATAC GeneActivity (mirrors scjoint) ----
shared_genes = rna.var_names.intersection(atac_gene.var_names)
print(f"shared genes: {len(shared_genes)} | RNA cells: {rna.n_obs} | ATAC cells: {atac_gene.n_obs}")

# RNA shared: lognorm on shared genes
rna_sh = rna[:, shared_genes].copy()
sc.pp.normalize_total(rna_sh, target_sum=1e4)
sc.pp.log1p(rna_sh)
s1 = to_dense(rna_sh.X).astype(np.float32)

# ATAC shared: GeneActivity lognorm on the same shared genes
atac_sh = atac_gene[:, shared_genes].copy()
sc.pp.normalize_total(atac_sh, target_sum=1e4)
sc.pp.log1p(atac_sh)
s2 = to_dense(atac_sh.X).astype(np.float32)

# RNA active: HVG -> scale (full-RNA lognorm first)
rna_act = rna.copy()
sc.pp.normalize_total(rna_act, target_sum=1e4)
sc.pp.log1p(rna_act)
sc.pp.highly_variable_genes(rna_act, n_top_genes=2000)
rna_act = rna_act[:, rna_act.var.highly_variable].copy()
sc.pp.scale(rna_act)
a1 = to_dense(rna_act.X).astype(np.float32)

# ATAC active: integrated LSI embedding (drop-in for the Signac LSI in the RMS run)
a2 = to_dense(atac_lsi.X).astype(np.float32)

assert s1.shape[1] == s2.shape[1], "shared arrays must have the same feature columns"
print(f"RNA {s1.shape}/{a1.shape}  ATAC {s2.shape}/{a2.shape}")

# ---- MaxFuse Fusor (identical workflow + params to RMS/pbmc3k) ----
fusor = mf.model.Fusor(shared_arr1=s1, shared_arr2=s2,
                       active_arr1=a1, active_arr2=a2, labels1=None, labels2=None)
fusor.split_into_batches(max_outward_size=5000, matching_ratio=3, metacell_size=2, verbose=True)
fusor.construct_graphs(n_neighbors1=15, n_neighbors2=15,
                       svd_components1=30, svd_components2=30,
                       resolution1=2, resolution2=2, resolution_tol=0.1, verbose=True)
fusor.find_initial_pivots(wt1=0.7, wt2=0.7, svd_components1=25, svd_components2=20)
fusor.refine_pivots(wt1=0.7, wt2=0.7, svd_components1=30, svd_components2=30,
                    cca_components=20, n_iters=3, randomized_svd=False, svd_runs=1, verbose=True)
fusor.filter_bad_matches(target="pivot", filter_prop=0.3)

rna_cca, atac_cca = fusor.get_embedding(active_arr1=fusor.active_arr1,
                                        active_arr2=fusor.active_arr2)

joint = pd.concat([pd.DataFrame(np.asarray(rna_cca),  index=rna_bc),
                   pd.DataFrame(np.asarray(atac_cca), index=atac_bc)])
joint.to_csv(os.path.join(OUT, "latent.csv"))
stop = timeit.default_timer()
with open(os.path.join(OUT, "maxfuse_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={joint.shape}  time={stop-start:.1f}s")
