#!/usr/bin/env python3
"""
Stage all 14 RMS method latents (rep1/rep2/rep3) normalized to label.csv's <bc>-1_rna / <bc>-1_atac,
so FigS4A (00_compute_metrics.py) and Fig4 (02_reproduce_metrics.py) can read uniform staged/<method>/rep{1,2,3}.csv.

Per-method barcode quirks (verified against label.csv):
  clean   : already <bc>-1_rna / <bc>-1_atac (+ maybe extra train/non-test rows the metrics ignore)
  cobolt  : "test_atac~<bc>-1_atac" etc. -> strip "^.*~"
  bindsc  : R make.unique -> "<bc>.1" (rna), "<bc>.1.1" (atac)
  seurat  : "<bc>-1_1" (rna), "<bc>-1_2" (atac)
  scjoint : space-delimited embeddings with NO ids; rna_bc.txt / atac_bc.txt give per-row barcodes
  scVI    : already fixed to scvi_latent.csv by scVI/fix_scvi_barcodes.py

Re-run after simba + scButterfly finish (their reps are SKIPPED with a notice until present).
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))            # .../RMS/benchmark
SCRIPT = os.path.normpath(os.path.join(HERE, "..", "Mast607", "script"))
STAGE = os.path.join(HERE, "staged")
SCJOINT_BC = "scJoint"   # rna_bc.txt / atac_bc.txt (shared across reps)

# display name -> (kind, [rep1, rep2, rep3] paths relative to SCRIPT).
# NOTE: Fig4 uses rep1/rep2/rep4 for the three seed-fragility methods (scVI, MIDAS, scDART) -- the 3rd
# slot below points to their rep4 run (seed 13). rep3/rep5 exist too; rep4 is included so the triplicate
# captures scDART's seed sensitivity (rep1/2/3 was the lucky combo: 0.99 -> 0.84 with rep4).
METHODS = {
    "scVI":             ("clean",  ["scVI/scvi_latent.csv", "scVI/res_Mast607A/rep2/scvi_latent.csv", "scVI/res_Mast607A/rep3/scvi_latent.csv", "scVI/res_Mast607A/rep4/scvi_latent.csv", "scVI/res_Mast607A/rep5/scvi_latent.csv"]),
    "Cobolt":           ("cobolt", ["cobolt/res_Mast607A/latent.csv", "cobolt/res_Mast607A/rep2/latent.csv", "cobolt/res_Mast607A/rep3/latent.csv", "cobolt/res_Mast607A/rep4/latent.csv", "cobolt/res_Mast607A/rep5/latent.csv"]),
    "Portal":           ("clean",  ["portal/lat_df.csv", "portal/rep2/lat_df.csv", "portal/rep3/lat_df.csv", "portal/rep4/lat_df.csv", "portal/rep5/lat_df.csv"]),
    "scglue":           ("clean",  ["scglue/latent.csv", "scglue/rep2/latent.csv", "scglue/rep3/latent.csv"]),
    "scglue(multiome)": ("clean",  ["scglue_withpair/scglue_latent.csv", "scglue_withpair/rep2/scglue_latent.csv", "scglue_withpair/rep3/scglue_latent.csv"]),
    "BindSC":           ("bindsc", ["bindsc/res_Mast607A/coembed_coor.csv", "bindsc/res_Mast607A/rep2/coembed_coor.csv", "bindsc/res_Mast607A/rep3/coembed_coor.csv", "bindsc/res_Mast607A/rep4/coembed_coor.csv", "bindsc/res_Mast607A/rep5/coembed_coor.csv"]),
    "Seurat(CCA)":      ("seurat", ["seurat/res_Mast607A/coembed_coor.csv", "seurat/res_Mast607A/rep2/coembed_coor.csv", "seurat/res_Mast607A/rep3/coembed_coor.csv"]),
    "MaxFuse":          ("clean",  ["maxfuse/latent.csv", "maxfuse/rep2/latent.csv", "maxfuse/rep3/latent.csv"]),
    "MIDAS":            ("clean",  ["midas/latent.csv", "midas/rep2/latent.csv", "midas/rep3/latent.csv", "midas/rep4/latent.csv", "midas/rep5/latent.csv"]),
    "scDART":           ("clean",  ["scDART/latent.csv", "scDART/rep2/latent.csv", "scDART/rep3/latent.csv", "scDART/rep4/latent.csv", "scDART/rep5/latent.csv"]),
    "scBridge":         ("clean",  ["scBridge/latent.csv", "scBridge/rep2/latent.csv", "scBridge/rep3/latent.csv"]),
    "simba":            ("clean",  ["simba/latent.csv", "simba/rep2/latent.csv", "simba/rep3/latent.csv"]),
    "scButterfly":      ("clean",  ["scbutterfly/latent.csv", "scbutterfly/rep4/latent.csv", "scbutterfly/rep5/latent.csv"]),  # forced-seed triplet (rep2/3 byte-identical -> artificial 1.0; rep4/5 FORCE_SEED 100/200 -> ~0.959)
    "scJoint":          ("scjoint", ["scJoint/output", "scJoint/rep2/output", "scJoint/rep3/output", "scJoint/rep4/output", "scJoint/rep5/output", "scJoint/rep6/output", "scJoint/rep7/output"]),
}


def norm_clean(df):
    df.index = df.index.astype(str).str.strip('"')
    return df


def norm_cobolt(df):
    df.index = df.index.astype(str).str.strip('"').str.replace(r"^.*~", "", regex=True)
    return df


def norm_seurat(df):
    df.index = (df.index.astype(str).str.strip('"')
                .str.replace(r"_1$", "_rna", regex=True)
                .str.replace(r"_2$", "_atac", regex=True))
    return df


def norm_bindsc(df):
    def fix(b):
        b = b.strip('"')
        if b.endswith(".1.1"):
            return b[:-4] + "-1_atac"
        if b.endswith(".1"):
            return b[:-2] + "-1_rna"
        return b
    df.index = [fix(b) for b in df.index.astype(str)]
    return df


def load_scjoint(rel_outdir):
    out = os.path.join(SCRIPT, rel_outdir)
    rna = pd.read_csv(os.path.join(out, "rna_scjoint_embeddings.txt"), sep=r"\s+", header=None)
    atac = pd.read_csv(os.path.join(out, "atac_gene_scjoint_embeddings.txt"), sep=r"\s+", header=None)
    rbc = [l.strip() for l in open(os.path.join(SCRIPT, SCJOINT_BC, "rna_bc.txt"))]
    abc = [l.strip() for l in open(os.path.join(SCRIPT, SCJOINT_BC, "atac_bc.txt"))]
    rna.index = [f"{b}_rna" for b in rbc]        # bc files are <bc>-1 -> <bc>-1_rna
    atac.index = [f"{b}_atac" for b in abc]
    df = pd.concat([rna, atac])
    df.columns = [f"latent_{i}" for i in range(df.shape[1])]
    return df


NORM = {"clean": norm_clean, "cobolt": norm_cobolt, "seurat": norm_seurat, "bindsc": norm_bindsc}


def main():
    n_ok, n_miss = 0, 0
    for method, (kind, paths) in METHODS.items():
        outdir = os.path.join(STAGE, method)
        os.makedirs(outdir, exist_ok=True)
        for i, p in enumerate(paths, 1):
            try:
                if kind == "scjoint":
                    if not os.path.exists(os.path.join(SCRIPT, p, "rna_scjoint_embeddings.txt")):
                        print(f"  MISSING {method} rep{i}: {p}"); n_miss += 1; continue
                    df = load_scjoint(p)
                else:
                    full = os.path.join(SCRIPT, p)
                    if not os.path.exists(full):
                        print(f"  MISSING {method} rep{i}: {full}"); n_miss += 1; continue
                    df = NORM[kind](pd.read_csv(full, index_col=0))
            except Exception as e:
                print(f"  ERROR  {method} rep{i}: {e}"); n_miss += 1; continue
            df.to_csv(os.path.join(outdir, f"rep{i}.csv"))
            print(f"  {method:18s} rep{i}: {df.shape[0]:6d} cells -> staged/{method}/rep{i}.csv")
            n_ok += 1
    print(f"DONE: staged {n_ok} latents, {n_miss} missing/errored "
          f"(re-run after simba + scButterfly finish to fill the gaps)")


if __name__ == "__main__":
    main()
