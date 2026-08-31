#!/usr/bin/env python
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Builder for the ArchR gene-activity VARIANT of the portal atac_gene h5.
#
# Portal does not call Signac::GeneActivity; it reads a precomputed CellRanger-format
# atac_gene .h5 (cells x genes, feature_types == 'Gene Expression') via
# sc.read_10x_h5(..., gex_only=False). For the R2.1 gene-activity sensitivity test we
# replace that precomputed Signac gene-activity matrix with the precomputed ArchR
# GeneScoreMatrix exported at .../geneactivity_archr/archr/export/.
#
# The ArchR export is genes x cells (MatrixMarket .mtx + genes.csv + barcodes.csv).
# Here we read it, transpose to CELLS x GENES, and write a 10x CellRanger v3 .h5
# with every feature tagged 'Gene Expression' so that sc.read_10x_h5(gex_only=False)
# returns the whole matrix. Barcodes are kept raw (<bc>-1), exactly as ArchR exported
# them; the original portal atac_gene h5 used the same raw barcodes and portal treats
# the ATAC barcodes only as cell labels (it adds an "_atac" suffix downstream).
#
# Run this (in the portal conda env, which has scanpy/scipy/h5py) BEFORE 01_run_portal.py.

import os
import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

export_dir = "/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/archr/export"
out_h5 = "/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/portal/atac_gene_archr.h5"

mtx_path = os.path.join(export_dir, "archr_gene_scores.mtx")
genes_path = os.path.join(export_dir, "genes.csv")
barcodes_path = os.path.join(export_dir, "barcodes.csv")

# ArchR export: genes x cells
print("reading", mtx_path)
m = scipy.io.mmread(mtx_path)            # genes x cells (sparse COO)
genes = pd.read_csv(genes_path)["gene"].astype(str).tolist()
barcodes = pd.read_csv(barcodes_path)["barcode"].astype(str).tolist()
print("loaded matrix genes x cells:", m.shape, "| n_genes:", len(genes), "| n_cells:", len(barcodes))
assert m.shape[0] == len(genes), "gene dim mismatch"
assert m.shape[1] == len(barcodes), "cell dim mismatch"

# --- restrict to the 1642 TEST-ATAC barcodes (the exact cell set the other methods embed), 0-filling
#     any test cell ArchR's minFrags QC dropped. Without this, portal embeds ALL 2748 ArchR cells
#     (train+test) instead of the 1642 test split, inflating its ATAC count vs every other method. ---
test_bc_path = os.path.join(os.path.dirname(out_h5), "test_atac_barcodes.txt")
test_bc = pd.read_csv(test_bc_path, header=None)[0].astype(str).tolist()
pos = {b: i for i, b in enumerate(barcodes)}
present = [b for b in test_bc if b in pos]
missing = [b for b in test_bc if b not in pos]
m = sp.csc_matrix(m)
sub = m[:, [pos[b] for b in present]]                       # genes x present
if missing:
    zeros = sp.csc_matrix((m.shape[0], len(missing)), dtype=m.dtype)
    m = sp.hstack([sub, zeros], format="csc")               # present cells, then 0-filled missing
    barcodes = present + missing
else:
    m, barcodes = sub, present
print(f"restricted to {len(barcodes)} test-ATAC cells ({len(present)} from ArchR, {len(missing)} filled 0)")
assert m.shape[1] == len(barcodes)

# The 10x CellRanger h5 schema stores features(rows) x barcodes(cols) in CSC layout
# (matrix/shape == [n_features, n_barcodes]); sc.read_10x_h5 then transposes on read so
# the resulting AnnData is CELLS x GENES. The ArchR export m is already genes(features) x
# cells, i.e. exactly the features x barcodes orientation the h5 wants -- so we store it
# directly (no transpose needed) and scanpy delivers cells x genes downstream.
mat = sp.csc_matrix(m)                    # genes(features) x cells, CSC over feature rows
mat = mat.astype(np.float32)             # ArchR gene scores are real-valued
mat.eliminate_zeros()

n_features, n_barcodes = mat.shape
print("writing 10x h5: features x barcodes =", (n_features, n_barcodes))

feature_ids = np.array(genes, dtype="S")
feature_names = np.array(genes, dtype="S")
feature_types = np.array(["Gene Expression"] * n_features, dtype="S")
genome = np.array(["GRCh38"] * n_features, dtype="S")
bc = np.array(barcodes, dtype="S")

os.makedirs(os.path.dirname(out_h5), exist_ok=True)
if os.path.exists(out_h5):
    os.remove(out_h5)

with h5py.File(out_h5, "w") as f:
    grp = f.create_group("matrix")
    grp.create_dataset("data", data=mat.data.astype(np.float32))
    grp.create_dataset("indices", data=mat.indices.astype(np.int64))
    grp.create_dataset("indptr", data=mat.indptr.astype(np.int64))
    grp.create_dataset("shape", data=np.array([n_features, n_barcodes], dtype=np.int32))
    grp.create_dataset("barcodes", data=bc)
    fea = grp.create_group("features")
    fea.create_dataset("id", data=feature_ids)
    fea.create_dataset("name", data=feature_names)
    fea.create_dataset("feature_type", data=feature_types)
    fea.create_dataset("genome", data=genome)
    fea.create_dataset("_all_tag_keys", data=np.array([b"genome"]))

print("wrote", out_h5)

# Sanity check: re-read with scanpy exactly as 01_run_portal.py will.
try:
    import scanpy as sc
    a = sc.read_10x_h5(out_h5, gex_only=False)
    print("scanpy re-read OK:", a.shape, "| feature_types unique:", a.var["feature_types"].unique().tolist())
except Exception as e:
    print("WARNING: scanpy re-read check failed:", repr(e))
