#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# =============================================================================
# 01_azimuth_annotate_parse.R
# Reproduce the Parse "Evercode WT Mini v3 PBMC" Figure-1 annotation for ONE donor:
#   load Parse data -> QC -> Seurat embedding -> Azimuth (pbmcref) cell typing
#   -> harmonize to the benchmark's 6-type taxonomy -> export label CSV + UMAPs.
#
# Why Azimuth (not transfer from our pbmc3k): Azimuth uses an EXTERNAL reference
# (Hao 2021), so the Parse RNA labels are independent of our 10x reference that
# defines the ATAC labels -> keeps the cross-platform test genuinely independent
# (answers the R4 circularity concern), and matches the published Parse Fig 1.
#
# Run PER DONOR (change DONOR), compare the printed compositions, pick the donor
# whose 6-type makeup + cell count is closest to pbmc3k, use that donor's CSV.
#
# ENV: needs Seurat (v4/5) + Azimuth + SeuratData(pbmcref). RunAzimuth downloads
# the pbmcref reference on first use -> run on a node WITH internet (login/interactive),
# not a compute node (same no-network trap we hit with BABEL).
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)     # provides RunAzimuth() + the "pbmcref" reference
  library(dplyr)
  library(ggplot2)
})

set.seed(1234)

# Seurat v5 stores assays as "Assay5" (layered), which the v4-era Azimuth rejects with
# 'invalid class "Assay5" object: Layers must be two-dimensional objects'. Force the
# classic v3 assay so CreateSeuratObject + RunAzimuth get the structure Azimuth expects.
options(Seurat.object.assay.version = "v3")

## ----------------------------- 0. paths / args ----------------------------
ROOT      <- "/path/to/multiomeBench"
PARSE_DIR <- file.path(ROOT, "pbmc_parse", "data")   # <-- EDIT to the new 4-donor ReadParseBio dir
# DONOR from the command line (e.g. `Rscript 01_azimuth_annotate_parse.R Donor_2`); default Donor_1
.args     <- commandArgs(trailingOnly = TRUE)
DONOR     <- if (length(.args) >= 1) .args[1] else "Donor_1"
OUT_DIR   <- file.path(ROOT, "pbmc_parse", "azimuth")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ----------------------------- 1. load Parse, subset donor ----------------
mat <- ReadParseBio(PARSE_DIR)
rownames(mat)[rownames(mat) == ""] <- "unknown"
meta <- read.csv(file.path(PARSE_DIR, "cell_metadata.csv"), row.names = 1)

obj <- CreateSeuratObject(counts = mat, meta.data = meta, names.field = 0)
stopifnot("sample" %in% colnames(obj@meta.data))
obj <- subset(obj, subset = sample == DONOR)
message(sprintf("[%s] loaded %d cells x %d genes", DONOR, ncol(obj), nrow(obj)))

## ----------------------------- 2. QC (mirror 00_data_prep.R) -----------------
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nCount_RNA < 50000 &
                            nCount_RNA > 500 & percent.mt < 20)
message(sprintf("[%s] after QC: %d cells", DONOR, ncol(obj)))

## ----------------------------- 3. Seurat embedding (for a Fig1-like UMAP) --
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)   # the object's OWN umap (reduction = "umap")

## ----------------------------- 4. Azimuth typing (Parse Fig 1 method) -----
# Adds: predicted.celltype.l1/l2, predicted.celltype.l1/l2.score, and a ref.umap.
# Safety: if anything upstream left a v5 Assay5, convert it back to v3 for Azimuth.
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
obj <- RunAzimuth(obj, reference = "pbmcref")

## ----------------------------- 5. harmonize Azimuth L2 -> pbmc3k taxonomy --
# Maps Azimuth L2 -> the EXACT pbmc3k labels (from pbmc3k_celltype_annotation.csv):
#   CD14 Monocytes, CD16 Monocytes, B cells, Naive CD4 T cells, Naive CD8 T cells,
#   Memory T cells, NK/Effector T cells, Myeloid DC, Plasmacytoid DC.
# The ONE judgment call left: whether cytotoxic CD8 TEM / CD4 CTL belong with
# "NK/Effector T cells" (default here, matching the class name) or "Memory T cells".
L2_TO_TYPE <- c(
  # --- Monocytes (CD16 kept SEPARATE) ---
  "CD14 Mono" = "CD14 Monocytes",
  "CD16 Mono" = "CD16 Monocytes",
  # --- B ---
  "B naive" = "B cells", "B intermediate" = "B cells", "B memory" = "B cells",
  "Plasmablast" = "B cells",
  # --- Naive T ---
  "CD4 Naive" = "Naive CD4 T cells",
  "CD8 Naive" = "Naive CD8 T cells",
  # --- Memory / non-cytotoxic T ---
  "CD4 TCM" = "Memory T cells", "CD4 TEM" = "Memory T cells", "CD8 TCM" = "Memory T cells",
  "Treg" = "Memory T cells", "MAIT" = "Memory T cells",
  "CD4 Proliferating" = "Memory T cells", "CD8 Proliferating" = "Memory T cells",
  # --- NK + cytotoxic/effector T (CD8 TEM/CD4 CTL here; move to Memory T if your
  #     MOFA+ annotation grouped them there instead) ---
  "NK" = "NK/Effector T cells", "NK_CD56bright" = "NK/Effector T cells",
  "NK Proliferating" = "NK/Effector T cells",
  "CD8 TEM" = "NK/Effector T cells", "CD4 CTL" = "NK/Effector T cells",
  # --- Dendritic cells (real classes in your set, not "Other") ---
  "cDC1" = "Myeloid DC", "cDC2" = "Myeloid DC", "ASDC" = "Myeloid DC",
  "pDC" = "Plasmacytoid DC"
  # everything else (Eryth/HSPC/Platelet/gdT/dnT/ILC) -> "Other"
)

l2 <- as.character(obj$predicted.celltype.l2)
obj$celltype <- unname(L2_TO_TYPE[l2])
obj$celltype[is.na(obj$celltype)] <- "Other"   # not one of the pbmc3k types

message("[composition] harmonized 6-type makeup for ", DONOR, ":")
print(sort(table(obj$celltype), decreasing = TRUE))

## ----------------------------- 6. export label CSV ------------------------
# `celltype` = harmonized 6-type label to feed the benchmark.
# Keep the Azimuth fine label + score too: the score is a per-cell annotation-
# confidence (useful for the R1.6a annotation-uncertainty analysis later).
labels <- data.frame(
  Barcode       = colnames(obj),
  celltype      = obj$celltype,
  azimuth_l2    = obj$predicted.celltype.l2,
  azimuth_l1    = obj$predicted.celltype.l1,
  azimuth_score = obj$predicted.celltype.l2.score,
  row.names = NULL, stringsAsFactors = FALSE
)
label_csv <- file.path(OUT_DIR, sprintf("parse_%s_azimuth_labels.csv", DONOR))
write.csv(labels, label_csv, row.names = FALSE)
message("wrote labels -> ", label_csv)

# composition table (for the donor-selection comparison across the 4 donors)
comp <- as.data.frame(table(obj$celltype)); colnames(comp) <- c("celltype", "n")
comp$frac <- round(comp$n / sum(comp$n), 4)
write.csv(comp, file.path(OUT_DIR, sprintf("parse_%s_composition.csv", DONOR)), row.names = FALSE)

## ----------------------------- 7. Fig1-like UMAPs -------------------------
p1 <- DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE,
              repel = TRUE) + ggtitle(sprintf("Parse %s - harmonized 6 types", DONOR))
p2 <- DimPlot(obj, reduction = "umap", group.by = "predicted.celltype.l2",
              label = TRUE, repel = TRUE) + ggtitle(sprintf("Parse %s - Azimuth L2", DONOR)) +
      theme(legend.text = element_text(size = 6))
ggsave(file.path(OUT_DIR, sprintf("parse_%s_umap_6type.png", DONOR)), p1, width = 8, height = 6, dpi = 200)
ggsave(file.path(OUT_DIR, sprintf("parse_%s_umap_azimuthL2.png", DONOR)), p2, width = 9, height = 6, dpi = 200)
message("wrote UMAPs -> ", OUT_DIR)

saveRDS(obj, file.path(OUT_DIR, sprintf("parse_%s_annotated.rds", DONOR)))
message("done.")
