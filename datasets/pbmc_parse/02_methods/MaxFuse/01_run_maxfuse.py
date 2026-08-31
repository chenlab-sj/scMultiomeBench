# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
MaxFuse on pbmc_parse cross-platform (unpaired RNA + ATAC) -> latent.csv.

MaxFuse IS a diagonal/unpaired integration method (category 1): it matches cells
across modalities and returns a shared CCA embedding for BOTH -> exactly the joint
space the benchmark scores. No pairing info is used.

CROSS-PLATFORM DIFFERENCE vs pbmc3k: there RNA and ATAC were a modality split of the
SAME multiome cells; HERE the RNA is independent Parse cells (1481) and the ATAC is
pbmc3k cells (1642). Because MaxFuse never uses pairing, the run code is unchanged
from the pbmc3k version -- it just consumes the 4 h5ads built by 00_prep_maxfuse_input.R,
whose obs_names already carry the modality suffixes (RNA -> _rna, ATAC -> _atac).

Inputs come from 00_prep_maxfuse_input.R (Signac), 4 h5ads:
  maxfuse_rna_shared.h5ad  : RNA lognorm on shared genes   (rows = <bc>_rna)
  maxfuse_atac_shared.h5ad : ATAC gene-activity, same genes (rows = <bc>_atac)
  maxfuse_rna_active.h5ad  : RNA scaled HVG                 (rows = <bc>_rna)
  maxfuse_atac_active.h5ad : ATAC LSI                       (rows = <bc>_atac)
shared_arr1/2 must share the SAME feature columns; active arrays are modality-specific.

API note: parameters below are MaxFuse's tutorial defaults (the benchmark runs every
method at defaults), copied verbatim from the pbmc3k 01_run_maxfuse.py.
"""
import os
import timeit
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import maxfuse as mf

ROOT = "/path/to/multiomeBench"
INP  = f"{ROOT}/pbmc_parse/MaxFuse/input"
OUT  = f"{ROOT}/pbmc_parse/MaxFuse"
SEED = 420
np.random.seed(SEED)
start = timeit.default_timer()

def arr(h5):
    a = sc.read_h5ad(h5)
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    return np.asarray(X, dtype=np.float32), a.obs_names

s1, rna_bc  = arr(f"{INP}/maxfuse_rna_shared.h5ad")
s2, atac_bc = arr(f"{INP}/maxfuse_atac_shared.h5ad")
a1, _       = arr(f"{INP}/maxfuse_rna_active.h5ad")
a2, _       = arr(f"{INP}/maxfuse_atac_active.h5ad")
assert s1.shape[1] == s2.shape[1], "shared arrays must have the same feature columns"
print(f"RNA {s1.shape}/{a1.shape}  ATAC {s2.shape}/{a2.shape}")

# ---- MaxFuse Fusor workflow (tutorial defaults) ----
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

# ---- joint latent.csv: RNA cells (_rna) + ATAC cells (_atac) in one CCA space ----
joint = pd.concat([pd.DataFrame(np.asarray(rna_cca),  index=rna_bc),
                   pd.DataFrame(np.asarray(atac_cca), index=atac_bc)])
os.makedirs(OUT, exist_ok=True)
joint.to_csv(os.path.join(OUT, "latent.csv"))
stop = timeit.default_timer()
with open(os.path.join(OUT, "maxfuse_runtime.txt"), "w") as fh:
    fh.write(str(stop - start))
print(f"wrote latent.csv  shape={joint.shape}  time={stop-start:.1f}s")
