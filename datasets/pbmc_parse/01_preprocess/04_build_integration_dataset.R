#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# =============================================================================
# 04_build_integration_dataset.R
# Build a COMPLETE, self-contained integration-input dir (input_azimuth/) for the
# cross-platform run:  RNA = Parse donor (Azimuth-annotated, QC'd, "Other" dropped),
# ATAC = pbmc3k multiome test ATAC (unchanged -> copied in).
#
# Uses the SAME filenames as the original input/ so every method script needs only
# input/ -> input_azimuth/. Two real changes vs 00_data_prep.R:
#   (1) labels from the INDEPENDENT Azimuth annotation (not the circular pbmc3k transfer)
#   (2) cells whose harmonized type is "Other" are removed.
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(Matrix); library(reticulate) })
use_python("~/.conda/envs/seurat4/bin/python"); library(anndata)
source("/path/to/multiomeBench/common/0_create10x_datafmt.R")   # sub_bc_matrix / sub_h5 (as 00_data_prep.R)

ROOT       <- "/path/to/multiomeBench"
DONOR      <- "Donor_1"
ANN_RDS    <- file.path(ROOT, "pbmc_parse", "azimuth", sprintf("parse_%s_annotated.rds", DONOR))
OLD_INPUT  <- file.path(ROOT, "pbmc_parse", "input")                 # unchanged ATAC/train files + old labels live here
OUT        <- file.path(ROOT, "pbmc_parse", "input_azimuth"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PARSE_DATA <- file.path(ROOT, "pbmc_parse", "data")
TRAIN_FEATURES <- "/path/to/data/pbmc3k/filtered_feature_bc_matrix/features.tsv"

## ---------------- 1. RNA: Azimuth donor (QC'd), drop "Other" -----------------
rna <- readRDS(ANN_RDS); DefaultAssay(rna) <- "RNA"
stopifnot("celltype" %in% colnames(rna@meta.data))
rna <- subset(rna, subset = celltype != "Other")
message(sprintf("[RNA] %s: %d cells after dropping 'Other'", DONOR, ncol(rna)))
print(sort(table(rna$celltype), decreasing = TRUE))

## ---------------- 2. gene alignment: Parse genes overlapping 10x -------------
train_feature <- read.csv(TRAIN_FEATURES, sep = "\t", header = FALSE)
gene_feature  <- subset(train_feature, V3 == "Gene Expression"); rownames(gene_feature) <- gene_feature$V1
parse_feature <- read.csv(file.path(PARSE_DATA, "all_genes.csv"))
stopifnot(nrow(parse_feature) == nrow(rna))                 # gene order must line up (min.cells=0 -> no gene dropped)
keep <- !is.na(match(parse_feature$gene_id, gene_feature$V1))
parse_feature_10x     <- parse_feature[keep, ]
common.features       <- rownames(rna)[keep]
parse_feature_10x.fmt <- gene_feature[parse_feature_10x$gene_id, ]; rownames(parse_feature_10x.fmt) <- NULL
rna_sub <- subset(rna, features = common.features)
message(sprintf("[RNA] kept %d genes overlapping 10x", length(common.features)))

## ---------------- 3. labels: Azimuth RNA + (unchanged) pbmc3k ATAC -----------
old_lab   <- read.csv(file.path(OLD_INPUT, "muti_celltype.csv"), stringsAsFactors = FALSE)
atac_rows <- old_lab[grepl("scATAC", old_lab$modality), c("bc", "atac_celltype", "modality")]
rna_rows  <- data.frame(bc = colnames(rna_sub), rna_celltype = rna_sub$celltype, modality = "parse scRNA")
write.csv(bind_rows(rna_rows, atac_rows), file.path(OUT, "muti_celltype.csv"), row.names = FALSE)

## ---------------- 4. RNA exports (OLD filenames so scripts only change the dir) ----
rna_counts <- GetAssayData(rna_sub, assay = "RNA", slot = "counts")
write_h5ad(AnnData(X = t(as.matrix(rna_counts))), file.path(OUT, "pbmc_parsed1_rna.h5ad"))

rownames(rna_counts) <- parse_feature_10x.fmt$V2            # Parse symbols -> 10x symbols
write.table(parse_feature_10x.fmt, file.path(OUT, "tmp.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
sub_bc_matrix(barcodes = colnames(rna_counts), feature.mtx = file.path(OUT, "tmp.tsv"),
              count.mtx = rna_counts, output_dir = file.path(OUT, "pbmc_parse_d1/"))
sub_h5(barcodes = colnames(rna_counts), feature_mtx = file.path(OUT, "tmp.tsv"), count_mtx = rna_counts,
       path = file.path(OUT, "pbmc_parse_d1.h5"), test = TRUE, module = "Gene Expression")

rna_sub <- rna_sub %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
write_h5ad(AnnData(X = as.matrix(Embeddings(rna_sub, "pca"))), file.path(OUT, "pbmc_parse_testrna_embed.h5ad"))

## ---------------- 5. RData for the R methods (new all.rna + unchanged all.atac) ----
# R methods (Seurat3/bindSC/scBridge/scDART region2gene/scJoint) load pbmc_parse_seurat.RData.
# all.atac is unchanged -> take it from the old RData; all.rna is the new filtered Parse RNA.
ne <- new.env(); load(file.path(OLD_INPUT, "pbmc_parse_seurat.RData"), envir = ne)
all.atac <- ne$all.atac
all.rna  <- CreateSeuratObject(counts = GetAssayData(rna_sub, "RNA", "counts"))   # clean object, no Azimuth cruft
save(all.rna, all.atac, file = file.path(OUT, "pbmc_parse_seurat.RData"))

## ---------------- 6. copy the UNCHANGED ATAC / train inputs into input_azimuth ----
unchanged <- c("pbmc3k_test_atac.h5ad", "pbmc3k_test_atac_gene.h5ad", "pbmc3k_test_atac.h5",
               "pbmc_parse_testatac_embed.h5ad",
               "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_train.h5")
for (f in unchanged) file.copy(file.path(OLD_INPUT, f), file.path(OUT, f), overwrite = TRUE)
for (d in c("pbmc3k_test_atac", "pbmc3k_filtered_feature_bc_matrix_train"))
  file.copy(file.path(OLD_INPUT, d), OUT, recursive = TRUE, overwrite = TRUE)

message("DONE. Complete self-contained input dir -> ", OUT)
