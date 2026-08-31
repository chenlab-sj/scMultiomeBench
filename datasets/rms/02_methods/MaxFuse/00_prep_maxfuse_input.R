# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## 00_prep_maxfuse_input.R  (RMS Mast607A)
## Build MaxFuse inputs with the SAME Signac recipe as the pbmc3k MaxFuse prep: read the test multiome
## h5 (Gene Expression + common Peaks) + the ATAC fragments, then write 4 h5ads:
##   maxfuse_rna_shared.h5ad  : RNA lognorm on shared genes   (rows = <bc>_rna)
##   maxfuse_atac_shared.h5ad : ATAC gene-activity, same genes (rows = <bc>_atac)
##   maxfuse_rna_active.h5ad  : RNA scaled HVG                 (rows = <bc>_rna)
##   maxfuse_atac_active.h5ad : ATAC LSI                       (rows = <bc>_atac)
## Run in seurat4 env. See 00_prep_maxfuse_input_sub.sh.
## VERIFY first run: peak sep (CreateChromatinAssay sep=c(":","-")) matches the commonpeaks.h5 peak names;
## the fragments path resolves; GeneActivity returns a non-empty matrix.

library(reticulate); use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat); library(Signac); library(Matrix)
library(EnsDb.Hsapiens.v86); library(anndata)

DATA <- "/path/to/multiomeBench/RMS/Mast607/data/Mast607A_TB19_22652"
TEST_H5 <- file.path(DATA, "filtered_feature_bc_matrix_Mast607A_TB19_22652_commonpeaks.h5")  # RNA + common peaks
FRAG    <- file.path(DATA, "hg38", "outs", "atac_fragments.tsv.gz")
OUT     <- "/path/to/multiomeBench/RMS/Mast607/script/maxfuse/input"
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

gene.use <- intersect(VariableFeatures(all.rna), rownames(all.atac[["ACTIVITY"]]))
hvg      <- VariableFeatures(all.rna)
cat("shared genes:", length(gene.use), " | RNA HVG:", length(hvg), "\n")

rbc <- paste0(colnames(all.rna),  "_rna")
abc <- paste0(colnames(all.atac), "_atac")

save_h5ad <- function(mat, bc, path) {            # features x cells -> cells x features
  m <- t(as.matrix(mat)); rownames(m) <- bc
  write_h5ad(AnnData(X = m), path)
}
# Seurat v5 (Assay5): layers via [["assay"]]$layer
save_h5ad(all.rna[["RNA"]]$data[gene.use, ],          rbc, file.path(OUT, "maxfuse_rna_shared.h5ad"))
save_h5ad(all.atac[["ACTIVITY"]]$data[gene.use, ],    abc, file.path(OUT, "maxfuse_atac_shared.h5ad"))
save_h5ad(all.rna[["RNA"]]$scale.data[hvg, ],         rbc, file.path(OUT, "maxfuse_rna_active.h5ad"))

al <- Embeddings(all.atac, reduction = "lsi"); rownames(al) <- abc
write_h5ad(AnnData(X = as.matrix(al)), file.path(OUT, "maxfuse_atac_active.h5ad"))

cat("wrote 4 MaxFuse input h5ads to", OUT, "\n")
