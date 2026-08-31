#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## Fig5B new-method pileups: compute the MYOD1 + FOXO1 region ATAC coverage tracks (per predicted cell type)
## for the 3 NEW methods (MaxFuse/MIDAS/scButterfly), EXACTLY as the published ataclabel2peak.R did for the
## old methods (same fragments, same plot_region_pileups params). Also re-emits the "annotation" (true-label)
## track as a consistency reference. Output: new_methods_{myod1,foxo1}value.csv (position, Mesoderm, Myoblast,
## Myocyte, pipeline) -- spliced onto the published value CSVs for the Fig5B heatmap.  Run in seurat4 env.
source("/path/to/multiomeBench/common/SJRHB013758_X2/plot_predpeak/plot_regionpeak_fun.R")
suppressMessages({library(Seurat); library(Signac); library(dplyr)})
mypalette <- c("#1F77B4", "#FF7F0E", "#2CA02C")

ROOT <- "/path/to/multiomeBench"
FIG  <- file.path(ROOT, "RMS/benchmark/figS4a")                 # knn_pred_label__*.csv + label.csv
OUT  <- file.path(ROOT, "RMS/benchmark/fig5b")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fpath <- "/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/atac_fragments.tsv.gz"  # matches published

NEW <- c("MaxFuse", "MIDAS", "scButterfly")

## combine the 3 new methods' predicted labels (barcodes like AAAC...-1_atac)
KNN.label <- NULL
for (m in NEW) {
  d <- read.csv(file.path(FIG, paste0("knn_pred_label__", m, ".csv")), row.names = 1, check.names = FALSE)
  colnames(d) <- m
  KNN.label <- if (is.null(KNN.label)) d else cbind(KNN.label, d[rownames(KNN.label), , drop = FALSE])
}
rownames(KNN.label) <- gsub("_atac", "", rownames(KNN.label))
bc <- rownames(KNN.label)

## annotation (true cell type) for the reference track
annot <- read.csv(file.path(FIG, "label.csv"), row.names = 1)
rownames(annot) <- gsub("_atac", "", rownames(annot))
KNN.label$annotation <- annot[bc, "cell_type"]

## per-cell total ATAC reads (for depth normalisation inside plot_region_pileups)
atac.frag <- CountFragments(fragments = fpath) %>% dplyr::filter(CB %in% bc) %>% dplyr::select(CB, reads_count)
rownames(atac.frag) <- atac.frag$CB
reads <- atac.frag[bc, "reads_count"]

regions <- list(myod1 = c("11", 17725000, 17745000), foxo1 = c("13", 41100000, 41300000))
for (rn in names(regions)) {
  r <- regions[[rn]]
  peak_list <- list()
  for (m in c(NEW, "annotation")) {
    KNN_atac <- data.frame(cell = bc, group = KNN.label[[m]], total_reads = reads)
    message("  ", rn, " / ", m)
    p.tmp <- plot_region_pileups(chrom = r[1], start = as.numeric(r[2]), end = as.numeric(r[3]),
                                 bam_bed_file = fpath, grouping = KNN_atac,
                                 smoothing_window = 40, x_steps = 100000, gene_annot = FALSE,
                                 col_palette = mypalette, fixed_ylim = 10)
    pd <- p.tmp[[1]]; pd$pipeline <- m; peak_list[[m]] <- pd
  }
  out <- do.call(rbind, peak_list)
  write.csv(out, file.path(OUT, paste0("new_methods_", rn, "value.csv")), row.names = FALSE, quote = FALSE)
  message("wrote new_methods_", rn, "value.csv")
}
cat("DONE\n")
