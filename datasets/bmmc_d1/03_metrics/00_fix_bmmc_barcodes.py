#!/usr/bin/env python3
"""
Restore the gem-group '-1' that the mmread read-helper strips from BMMC latent barcodes, so
<bc>_rnaN -> <bc>-1_rnaN and the latent joins celltype.csv (which uses <bc>-1_rnaN).

The helper is a copy of scvi.data.read_10x_multiome, which drops the '-1' even though barcodes.tsv
carries it. maxfuse is unaffected (it used the pre-built h5ads that keep '-1'); midas/midas2/scbutterfly
read the raw 10x and need this. Idempotent: barcodes that already end in '-1' are left alone.

Run it after each deep method finishes:  python 00_fix_bmmc_barcodes.py            (does all 4)
                                          python 00_fix_bmmc_barcodes.py midas midas2
"""
import sys
import os
import re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def add_dash1(b):
    m = re.search(r'_(?:rna|atac)\d$', b)      # the _rna1.._atac3 suffix
    if not m:
        return b
    stem, suf = b[:m.start()], b[m.start():]
    return (stem if stem.endswith('-1') else stem + '-1') + suf


def main():
    methods = sys.argv[1:] or ["midas", "midas2", "scbutterfly", "maxfuse"]
    lab = set(pd.read_csv(os.path.join(HERE, "celltype.csv"), index_col=0).index.astype(str))
    for m in methods:
        f = os.path.join(HERE, m, "latent.csv")
        if not os.path.exists(f):
            print(f"  {m:12s}: no latent yet")
            continue
        d = pd.read_csv(f, index_col=0)
        before = len(set(d.index.astype(str)) & lab)
        d.index = [add_dash1(str(b)) for b in d.index]
        after = len(set(d.index) & lab)
        if after > before:
            d.to_csv(f)
            print(f"  {m:12s}: join {before} -> {after}/25207  (fixed + written)")
        else:
            print(f"  {m:12s}: join {after}/25207  (already ok, no change)")


if __name__ == "__main__":
    main()
