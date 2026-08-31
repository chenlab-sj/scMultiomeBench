#!/usr/bin/env python3
"""
Repair the two barcode-format bugs that made compute_metrics drop Seurat(CCA) & bindsc and let portal
embed the wrong cell set. Run on the Mac (pandas only).

1. Seurat(CCA) coembed_coor.csv index is <bc>_1 / <bc>_2  (Seurat merge convention) -> rename to
   <bc>_rna / <bc>_atac and write lat_df.csv.
2. bindsc coembed_coor.csv index is R-mangled  <bc>.1 (RNA) / <bc>.1.1 (ATAC, make.unique) -> recover
   the -1 barcode and write <bc>_rna / <bc>_atac to lat_df.csv.
   (compute_metrics keys cells off the <stem>_rna / <stem>_atac suffix; the latent VALUES are untouched,
    so this stays a clean Signac-vs-ArchR comparison.)
3. Extract the canonical 1642 test-ATAC barcodes (from a method that scored correctly, maxfuse) and
   cross-check against scJoint's atac_bc.txt -> portal/test_atac_barcodes.txt, the cell set portal's
   ArchR h5 must be restricted to (it currently embeds all 2748 ArchR cells).
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# labeled test barcodes (for sanity-checking the renames) -- 1531 _rna + 1531 _atac
lab = pd.read_csv(os.path.join(HERE, "metrics", "label.csv"), index_col=0)
test_bc = set(lab[lab["modality"] != "train multiomics"].index)


def relabel(in_csv, out_csv, fn, tag):
    df = pd.read_csv(in_csv, index_col=0)
    new = [fn(str(i)) for i in df.index]
    df.index = new
    df.to_csv(out_csv)
    nr = sum(b.endswith("_rna") for b in new)
    na = sum(b.endswith("_atac") for b in new)
    ov = len(set(new) & test_bc)
    print(f"  {tag:11s} {len(new)} rows -> _rna={nr} _atac={na} | overlap w/ labeled test = {ov} "
          f"(expect ~{len(test_bc)//2*2 if False else 3062}); wrote {os.path.relpath(out_csv, HERE)}")
    return set(new)


def seurat_fn(i):
    if i.endswith("_1"):
        return i[:-2] + "_rna"
    if i.endswith("_2"):
        return i[:-2] + "_atac"
    return i


def bindsc_fn(i):
    # ATAC first: <bc>.1.1 = (<bc>-1 mangled to <bc>.1) + make.unique ".1"
    if i.endswith(".1.1"):
        stem = i[:-2]              # drop make.unique ".1" -> "<bc>.1"
        return stem[:-2] + "-1_atac"
    if i.endswith(".1"):          # RNA: <bc>.1  (-1 mangled to .1)
        return i[:-2] + "-1_rna"
    return i


print("== relabel Seurat(CCA) & bindsc latents ==")
relabel(os.path.join(HERE, "Seurat_CCA", "res_pbmc3k_testall", "coembed_coor.csv"),
        os.path.join(HERE, "Seurat_CCA", "res_pbmc3k_testall", "lat_df.csv"), seurat_fn, "Seurat(CCA)")
relabel(os.path.join(HERE, "bindsc", "res_pbmc3k", "coembed_coor.csv"),
        os.path.join(HERE, "bindsc", "res_pbmc3k", "lat_df.csv"), bindsc_fn, "bindsc")

print("\n== canonical 1642 test-ATAC barcodes for portal ==")
mf = pd.read_csv(os.path.join(HERE, "maxfuse", "latent.csv"), index_col=0)
mf_atac = sorted(i[:-5] for i in mf.index if str(i).endswith("_atac"))
sj = sorted(pd.read_csv(os.path.join(HERE, "scJoint", "scJoint_need", "atac_bc.txt"),
                        header=None)[0].astype(str))
same = set(mf_atac) == set(sj)
print(f"  maxfuse _atac = {len(mf_atac)} cells; scJoint atac_bc = {len(sj)}; identical set = {same}")
out = os.path.join(HERE, "portal", "test_atac_barcodes.txt")
pd.Series(mf_atac).to_csv(out, index=False, header=False)
print(f"  wrote {os.path.relpath(out, HERE)} ({len(mf_atac)} barcodes)")
