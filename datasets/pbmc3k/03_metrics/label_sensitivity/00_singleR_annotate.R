# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# R4/R1.6 label-circularity rebuttal — independent RNA-only re-annotation of pbmc3k with SingleR.
#
# LABEL PROVENANCE (corrected 2026-07-16): the pbmc3k/pbmc10k "ground-truth" labels are from the **10X
# Genomics R&D team** -- Methods: "Cell-type labels from the 10X Genomics R&D team were used as surrogate
# ground truths". They are a THIRD-PARTY annotation, NOT derived here with Seurat WNN or any other
# benchmarked method. (An earlier version of this file wrongly called them "WNN labels" -- Seurat(WNN) is a
# *benchmarked method* in this study, which is a different thing. Do not reintroduce that claim.)
#
# R4 asks whether labels "computationally derived from paired multiomics" could favour integration methods
# sharing the label-deriver's bias. SingleR is a reference-CORRELATION annotator: it shares no algorithmic
# basis with any benchmarked integration method (no anchors/CCA/graph) and uses RNA ONLY. If method
# rankings are unchanged when scored against SingleR labels, the rankings are not an artifact of the
# label source -- whatever pipeline 10X used to produce them.
#
# This script: load pbmc3k TEST GEX -> SingleR(ref = MonacoImmuneData, label.fine) -> map fine labels to
# the benchmark's 9-type scheme -> write a new label.csv (test cells relabelled; train rows untouched,
# they aren't scored) + a confusion matrix vs the 10X labels + the per-cell calls + the mapping table.
#
# NOTE celldex fetches MonacoImmuneData from ExperimentHub (internet). If compute nodes lack internet,
# run once on a login node:  R -e 'celldex::MonacoImmuneData()'  to cache to ~/.cache/R/ExperimentHub.
suppressMessages({
  library(Seurat); library(SingleR); library(celldex); library(SummarizedExperiment)
})

ROOT      <- "/path/to/multiomeBench"
TEST_H5   <- file.path(ROOT, "pbmc/pbmc3k/Data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5")
ORIG_LABEL<- file.path(ROOT, "pbmc/pbmc3k/benchmark/fig2b/label.csv")          # 10X R&D ground-truth labels
OUT       <- file.path(ROOT, "pbmc/pbmc3k/benchmark/label_sensitivity")

# ---- MonacoImmune label.fine -> benchmark 9-type scheme. The mapping, and the rationale for every
#      debatable call, live in fine_label_map.R -- shared with 02_remap_labels.R so the two cannot drift.
#      EDIT THE MAPPING THERE, not here. This script prints any fine label it can't map, so nothing is
#      silently dropped. ----
source(file.path(OUT, "fine_label_map.R"))                  # -> MAP

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 1) load TEST GEX, lognormalize for SingleR ----
counts <- Read10X_h5(TEST_H5)
rna    <- if (is.list(counts)) counts[["Gene Expression"]] else counts
seu    <- CreateSeuratObject(rna)
seu    <- NormalizeData(seu, verbose = FALSE)
logm   <- GetAssayData(seu, slot = "data")                 # genes x cells, lognorm
cat(sprintf("test GEX: %d genes x %d cells\n", nrow(logm), ncol(logm)))

# ---- 2) reference + SingleR ----
ref  <- celldex::MonacoImmuneData()
pred <- SingleR(test = logm, ref = ref, labels = ref$label.fine)
fine <- pred$labels
names(fine) <- colnames(logm)

# ---- 3) map fine -> 9-type; flag unmapped ----
unmapped <- setdiff(unique(fine), names(MAP))
if (length(unmapped)) cat("WARNING unmapped fine labels (-> 'other'):\n  ", paste(unmapped, collapse=", "), "\n")
mapped <- MAP[fine]; mapped[is.na(mapped)] <- "other"; names(mapped) <- names(fine)
cat("\nSingleR mapped-label distribution:\n"); print(table(mapped))

# per-cell calls (transparency)
write.csv(data.frame(barcode = names(fine), singleR_fine = fine, singleR_celltype = mapped,
                     score = apply(pred$scores, 1, max)),
          file.path(OUT, "singleR_per_cell.csv"), row.names = FALSE)
write.csv(data.frame(fine = names(MAP), benchmark = unname(MAP)),
          file.path(OUT, "celltype_map.csv"), row.names = FALSE)

# ---- 4) new label.csv: test rows get SingleR label; train rows unchanged ----
lab      <- read.csv(ORIG_LABEL, stringsAsFactors = FALSE, check.names = FALSE)
is_test  <- lab$modality != "train multiomics"
cell_bc  <- sub("_(rna|atac)$", "", lab$Barcode)
lab_new  <- lab
hit      <- is_test & cell_bc %in% names(mapped)
lab_new$cell_type[hit] <- mapped[cell_bc[hit]]
miss <- sum(is_test & !(cell_bc %in% names(mapped)))
if (miss) cat(sprintf("\n%d test rows had no SingleR call (kept original 10X label)\n", miss))
write.csv(lab_new, file.path(OUT, "label_singleR.csv"), row.names = FALSE)

# ---- 5) confusion 10X vs SingleR on the RNA-side test cells (one row per cell) ----
rna_rows <- lab$modality == "test scRNA"
truth <- lab$cell_type[rna_rows]; names(truth) <- cell_bc[rna_rows]
common <- intersect(names(truth), names(mapped))
conf <- table(TenX = truth[common], SingleR = mapped[common])
write.csv(as.data.frame.matrix(conf), file.path(OUT, "confusion_10x_vs_singleR.csv"))
agree <- mean(truth[common] == mapped[common])
cat(sprintf("\n=== agreement (10X-annotated ground truth vs SingleR, %d test cells): %.1f%% ===\n",
            length(common), 100*agree))
cat("per-type agreement:\n")
for (t in sort(unique(truth[common]))) {
  idx <- names(truth[common])[truth[common] == t]
  cat(sprintf("  %-20s %5.1f%% (n=%d)\n", t, 100*mean(mapped[idx] == t), length(idx)))
}
cat("\nwrote: label_singleR.csv, singleR_per_cell.csv, celltype_map.csv, confusion_10x_vs_singleR.csv\n")
