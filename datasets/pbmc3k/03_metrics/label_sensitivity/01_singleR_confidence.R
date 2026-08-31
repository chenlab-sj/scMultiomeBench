#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Re-run SingleR on the pbmc3k test GEX to obtain per-cell CONFIDENCE, which 00_singleR_annotate.R did not save:
#   pruned.labels (NA = low-confidence / ambiguous call, SingleR's own outlier-based pruning) and delta.next.
# Local-capable on the Mac: SingleR + celldex are installed; MonacoImmuneData is fetched from ExperimentHub
# (needs internet once, then cached). VERIFIES the fine labels reproduce the cached singleR_per_cell.csv
# before writing, so the confidence flags are consistent with the existing (cached) analysis.
#   Rscript 01_singleR_confidence.R   ->  singleR_confidence.csv
suppressMessages({library(Seurat); library(SingleR); library(celldex); library(SummarizedExperiment)})

ROOT <- if (dir.exists("/path/to/project")) "/path/to/multiomeBench" else
                                          "/path/to/multiomeBench"
TEST_H5 <- file.path(ROOT, "pbmc/pbmc3k/Data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5")
OUT     <- file.path(ROOT, "pbmc/pbmc3k/benchmark/label_sensitivity")

counts <- Read10X_h5(TEST_H5)
rna    <- if (is.list(counts)) counts[["Gene Expression"]] else counts
seu    <- NormalizeData(CreateSeuratObject(rna), verbose = FALSE)
logm   <- GetAssayData(seu, slot = "data")
cat(sprintf("test GEX: %d genes x %d cells\n", nrow(logm), ncol(logm)))

ref  <- celldex::MonacoImmuneData()
pred <- SingleR(test = logm, ref = ref, labels = ref$label.fine)

df <- data.frame(barcode      = rownames(pred),
                 singleR_fine = pred$labels,
                 pruned_label = pred$pruned.labels,          # NA => pruned (low confidence)
                 is_confident = !is.na(pred$pruned.labels),  # TRUE => SingleR kept it (confident)
                 delta_next   = pred$delta.next,             # gap to next-best label
                 score        = apply(pred$scores, 1, max),
                 stringsAsFactors = FALSE)

# ---- reproducibility check vs the cached run (fine labels must match, else the confidence is for a
#      different labelling and cannot be spliced in) ----
cache <- read.csv(file.path(OUT, "singleR_per_cell.csv"), stringsAsFactors = FALSE)
m <- merge(cache[, c("barcode", "singleR_fine")], df[, c("barcode", "singleR_fine")],
           by = "barcode", suffixes = c("_cache", "_rerun"))
agree <- mean(m$singleR_fine_cache == m$singleR_fine_rerun)
cat(sprintf("fine-label reproducibility vs cache: %.2f%% (%d/%d match)\n",
            100*agree, sum(m$singleR_fine_cache == m$singleR_fine_rerun), nrow(m)))
if (agree < 0.99) cat("!! WARNING: re-run diverges from the cached labels -- do NOT splice; investigate.\n")

write.csv(df, file.path(OUT, "singleR_confidence.csv"), row.names = FALSE)
cat(sprintf("confident (non-pruned): %d / %d (%.1f%%)\n",
            sum(df$is_confident), nrow(df), 100*mean(df$is_confident)))
cat("wrote singleR_confidence.csv\n")
