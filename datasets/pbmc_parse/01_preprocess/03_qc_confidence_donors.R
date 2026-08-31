#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# =============================================================================
# 03_qc_confidence_donors.R
# Reads the saved parse_<DONOR>_annotated.rds files and summarizes, per donor:
#   (1) single-cell QC  -- genes/cell, UMIs/cell, %mito
#   (2) Azimuth annotation CONFIDENCE -- per-cell L2 score, per cell type
# to finalize the donor choice beyond composition, and to feed the R1.6(a)
# annotation-uncertainty analysis (the per-cell score IS the uncertainty measure).
# No re-run of Azimuth; just loads the .rds. Needs Seurat (for readRDS of the object).
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })

ROOT    <- "/path/to/multiomeBench"
ANN_DIR <- file.path(ROOT, "pbmc_parse", "azimuth")
DONORS  <- c("Donor_1", "Donor_2", "Donor_3", "Donor_4")   # <-- EDIT to your sample names
CONF_CUT <- 0.5                                            # Azimuth L2 score => "high confidence"
OUT     <- file.path(ANN_DIR, "qc_confidence"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ----------------------------- gather metadata from each .rds -------------
md_all <- list()
for (d in DONORS) {
  f <- file.path(ANN_DIR, sprintf("parse_%s_annotated.rds", d))
  if (!file.exists(f)) { message("skip (missing): ", f); next }
  m <- readRDS(f)@meta.data
  sc_col <- intersect(c("predicted.celltype.l2.score", "azimuth_score"), colnames(m))[1]
  md_all[[d]] <- data.frame(donor = d,
                            nFeature = m$nFeature_RNA, nCount = m$nCount_RNA,
                            pct_mt = m$percent.mt, celltype = m$celltype,
                            score = m[[sc_col]], stringsAsFactors = FALSE)
}
MD <- bind_rows(md_all)

## ----------------------------- per-donor QC + confidence summary ----------
summ <- MD %>% group_by(donor) %>% summarise(
  n_cells       = n(),
  median_genes  = median(nFeature),
  median_umi    = median(nCount),
  median_pct_mt = round(median(pct_mt), 2),
  median_score  = round(median(score), 3),
  pct_high_conf = round(mean(score > CONF_CUT) * 100, 1),
  .groups = "drop"
) %>% arrange(desc(median_score))
print(summ); write.csv(summ, file.path(OUT, "donor_qc_confidence_summary.csv"), row.names = FALSE)

## ----------------------------- per-cell-type confidence -------------------
# which types are RELIABLY called (the hard ones -- naive CD4/CD8, CD16 Mono --
# will score lower; that low-confidence subset is exactly the R1.6a story).
ct_conf <- MD %>% group_by(donor, celltype) %>%
  summarise(n = n(), median_score = round(median(score), 3),
            pct_high = round(mean(score > CONF_CUT) * 100, 1), .groups = "drop")
write.csv(ct_conf, file.path(OUT, "celltype_confidence_by_donor.csv"), row.names = FALSE)

## ----------------------------- plots --------------------------------------
ggsave(file.path(OUT, "qc_genes_by_donor.png"),
  ggplot(MD, aes(donor, nFeature, fill = donor)) + geom_violin(scale = "width") +
    scale_y_log10() + labs(y = "genes / cell", x = NULL, title = "QC: genes per cell") + theme_bw(),
  width = 7, height = 5, dpi = 200)
ggsave(file.path(OUT, "qc_pctmt_by_donor.png"),
  ggplot(MD, aes(donor, pct_mt, fill = donor)) + geom_violin(scale = "width") +
    labs(y = "% mito", x = NULL, title = "QC: mitochondrial %") + theme_bw(),
  width = 7, height = 5, dpi = 200)
ggsave(file.path(OUT, "azimuth_score_by_donor.png"),
  ggplot(MD, aes(donor, score, fill = donor)) + geom_violin(scale = "width") +
    labs(y = "Azimuth L2 score", x = NULL, title = "Annotation confidence per donor") + theme_bw(),
  width = 7, height = 5, dpi = 200)
ggsave(file.path(OUT, "azimuth_score_by_celltype.png"),
  ggplot(MD, aes(celltype, score, fill = donor)) + geom_boxplot(outlier.size = 0.3) +
    coord_flip() + labs(y = "Azimuth L2 score", x = NULL, title = "Confidence per cell type") + theme_bw(),
  width = 9, height = 6, dpi = 200)
message("wrote ", OUT)
