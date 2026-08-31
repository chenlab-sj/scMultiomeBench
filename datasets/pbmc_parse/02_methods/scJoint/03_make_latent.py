#!/usr/bin/env python3
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# scJoint produces per-modality embeddings (output/*_scjoint_embeddings.txt), not a single
# latent.csv. This concatenates them (RNA then ATAC), indexes by the barcode lists (which
# already carry the _rna/_atac suffix), and writes latent.csv in the benchmark convention.
import os
import numpy as np

SC = "/path/to/multiomeBench/pbmc_parse/scJoint"
BC = "/path/to/multiomeBench/pbmc_parse/input_azimuth/scJoint"

REP  = os.environ.get("REP", "")          # "" = base run (seed 1); "rep2" = second-seed run
OUTD = os.path.join(SC, REP)
rna  = np.loadtxt(os.path.join(OUTD, "output/rna_scjoint_embeddings.txt"))
atac = np.loadtxt(os.path.join(OUTD, "output/atac_gene_scjoint_embeddings.txt"))
rbc  = [l.strip() for l in open(os.path.join(BC, "rna_bc.txt"))]
abc  = [l.strip() for l in open(os.path.join(BC, "atac_bc.txt"))]
assert rna.shape[0] == len(rbc), (rna.shape, len(rbc))
assert atac.shape[0] == len(abc), (atac.shape, len(abc))

emb = np.vstack([rna, atac]); bcs = rbc + abc
ncol = emb.shape[1]
with open(os.path.join(OUTD, "latent.csv"), "w") as f:
    f.write("," + ",".join(f"latent_{i}" for i in range(ncol)) + "\n")
    for bc, row in zip(bcs, emb):
        f.write(bc + "," + ",".join(repr(float(x)) for x in row) + "\n")
print(f"wrote latent.csv: {emb.shape[0]} x {ncol} (rna={len(rbc)}, atac={len(abc)})")
