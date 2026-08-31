# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## 00_prep_maxfuse_input.R
## Build MaxFuse inputs for pbmc3k using the SAME Signac recipe as bindsc.R:
##   RNA  : lognorm expression (shared genes) + scaled HVG (active)
##   ATAC : GeneActivity scores (shared genes) + LSI embedding (active)
## Test cells are a modality split of the same cells -> RNA cells get '_rna',
## ATAC cells get '_atac' so the joint latent.csv matches benchmark_metrics.py.
## Run in the existing `bindsc` conda env (has Seurat/Signac). See 00_prep_maxfuse_input_sub.sh.

# seurat4 env (like 00_scBridge_s0.R / 00_scjoint_s1.R): its python has the `anndata`
# module that reticulate needs for write_h5ad. Without this binding, write_h5ad fails.
library(reticulate); use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat); library(Signac); library(Matrix)
library(EnsDb.Hsapiens.v86); library(anndata)

## ---- inputs (confirm the fragments path; it lives on clusterhome, not in Data/) ----
DATA   <- "/path/to/multiomeBench/pbmc/pbmc3k/Data"
TEST_H5 <- file.path(DATA, "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5")
FRAG    <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"  # <-- confirm
OUT     <- "/path/to/multiomeBench/pbmc/pbmc3k/scripts/maxfuse/input"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

x <- Read10X_h5(TEST_H5)
all.rna  <- CreateSeuratObject(counts = x$`Gene Expression`, assay = "RNA")
all.atac <- CreateSeuratObject(
  counts = CreateChromatinAssay(counts = x$Peaks, sep = c(":", "-"), fragments = FRAG),
  assay = "peaks")

ann <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(ann) <- paste0("chr", seqlevels(ann)); genome(ann) <- "hg38"
Annotation(all.atac) <- ann

## RNA: lognorm + HVG + scale + PCA
all.rna <- NormalizeData(all.rna); all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna);     all.rna <- RunPCA(all.rna)

## ATAC: LSI (active) + GeneActivity on the RNA HVGs (shared)
all.atac <- RunTFIDF(all.atac); all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)
ga <- GeneActivity(all.atac, features = VariableFeatures(all.rna))
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = ga)
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- NormalizeData(all.atac)

## shared gene set present in BOTH modalities
gene.use <- intersect(VariableFeatures(all.rna), rownames(all.atac[["ACTIVITY"]]))
hvg      <- VariableFeatures(all.rna)
cat("shared genes:", length(gene.use), " | RNA HVG:", length(hvg), "\n")

## barcodes: same cells, split by modality
rbc <- paste0(colnames(all.rna),  "_rna")
abc <- paste0(colnames(all.atac), "_atac")

save_h5ad <- function(mat, bc, path) {            # mat: features x cells -> cells x features
  m <- t(as.matrix(mat)); rownames(m) <- bc
  write_h5ad(AnnData(X = m), path)
}
# Seurat v5 (Assay5): layers via [["assay"]]$layer, NOT the old @data slot.
save_h5ad(all.rna[["RNA"]]$data[gene.use, ],          rbc, file.path(OUT, "maxfuse_rna_shared.h5ad"))
save_h5ad(all.atac[["ACTIVITY"]]$data[gene.use, ],    abc, file.path(OUT, "maxfuse_atac_shared.h5ad"))
save_h5ad(all.rna[["RNA"]]$scale.data[hvg, ],         rbc, file.path(OUT, "maxfuse_rna_active.h5ad"))

al <- Embeddings(all.atac, reduction = "lsi"); rownames(al) <- abc   # cells x LSI (active)
write_h5ad(AnnData(X = as.matrix(al)), file.path(OUT, "maxfuse_atac_active.h5ad"))

cat("wrote 4 MaxFuse input h5ads to", OUT, "\n")
