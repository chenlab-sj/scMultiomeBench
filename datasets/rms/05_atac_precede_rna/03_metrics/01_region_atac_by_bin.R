#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## W7 step 2 (Signac / seurat4): per-cell ATAC accessibility over the MYOD1 & MYOG loci, from the Mast607A
## atac_fragments, depth-normalized. Joined to the pseudotime/bin table from 00_dpt_pseudotime.py.
##   Rscript 01_region_atac_by_bin.R
suppressMessages({library(Signac); library(GenomicRanges); library(Matrix); library(data.table)})

HERE <- "/path/to/multiomeBench/RMS/benchmark/atac_precede_rna"
OUTS <- "/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs"
FRAG <- file.path(OUTS, "atac_fragments.tsv.gz")
PBM  <- file.path(OUTS, "per_barcode_metrics.csv")

tab   <- read.csv(file.path(HERE, "pseudotime_table.csv"), stringsAsFactors = FALSE)  # has barcode_base, bin, cell_type
cells <- tab$barcode_base

## hg19/GRCh37 gene-body+promoter windows (this data is hg19 — BAM chr1=249,250,621; MYOG frags hg19>>hg38).
## seqnames WITHOUT 'chr' to match this fragments file (col1 = '1','11'). ~20 kb LOCUS-level windows to reduce
## scATAC sparsity. Myogenic panel: MYF5/MYOD1 (early MRFs) -> MYOG (late MRF) -> MYH3/ACTA1 (terminal structural).
## hg19 coords (cellranger hg19 GTF). Compact genes -> gene-body+promoter window (good ATAC signal).
## MEF2C (186 kb) is huge -> ~20 kb promoter/TSS window (TSS 5:88,199,922, -); ACTN2 (78 kb) -> ~30 kb TSS window.
## Small genes: MYOD1 11:17,741,115-17,743,678 (+); MYOG 1:203,052,260-203,055,164 (-); MYL1 2:211,154,874-
## 211,179,914 (-); TNNT3 11:1,940,792-1,959,936 (+); TNNT2 1:201,328,136-201,346,890 (-).
genes <- data.frame(
  name  = c("MYOD1",  "MYOG",     "MYH3",   "MEF2C"),
  chr   = c("11",     "1",        "17",     "5"),
  start = c(17725000, 203045000,  10528000, 88190000),   # MEF2C = TSS +/- ~10 kb (promoter; 186 kb gene)
  end   = c(17745000, 203062000,  10565000, 88210000),
  stringsAsFactors = FALSE)
regions <- GRanges(seqnames = genes$chr, ranges = IRanges(start = genes$start, end = genes$end))
frag <- CreateFragmentObject(path = FRAG, cells = cells)
fm   <- FeatureMatrix(fragments = frag, features = regions, cells = cells)   # rows follow `genes` order
cat("FeatureMatrix dim (features x cells):", dim(fm)[1], "x", dim(fm)[2], "\n")

present <- intersect(cells, colnames(fm))
fm <- fm[, present, drop = FALSE]

## per-cell total ATAC fragments for depth normalization
pbm <- fread(PBM)
bc_col  <- if ("barcode" %in% names(pbm)) "barcode" else names(pbm)[1]
cand    <- c("atac_fragments", "passed_filters", "atac_peak_region_fragments", "atac_raw_fragments")
tot_col <- cand[cand %in% names(pbm)][1]
cat("depth column:", tot_col, " | barcode column:", bc_col, "\n")
depth <- setNames(as.numeric(pbm[[tot_col]]), pbm[[bc_col]])
d <- depth[present]; d[is.na(d) | d <= 0] <- NA

res <- data.frame(barcode_base = present, stringsAsFactors = FALSE)
for (i in seq_len(nrow(genes))) {                                 # one <gene>_atac col per locus, in `genes` order
  res[[paste0(genes$name[i], "_atac")]] <- as.numeric(fm[i, ]) / d * 1e4   # fragments-in-region per 10k total
}
res <- merge(tab[, c("barcode_base", "bin", "cell_type", "dpt_pseudotime")], res, by = "barcode_base")
write.csv(res, file.path(HERE, "region_atac_percell.csv"), row.names = FALSE)
cat("wrote region_atac_percell.csv:", nrow(res), "cells\n")
for (g in genes$name)
  cat("  mean", g, "atac =", round(mean(res[[paste0(g, "_atac")]], na.rm = TRUE), 3), "\n")
