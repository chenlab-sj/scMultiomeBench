# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## 00_prep_maxfuse_input.R  (BRCA HT243B1-S1H4)
## Build MaxFuse inputs with the SAME Signac recipe as the pbmc3k MaxFuse prep, but loading
## BRCA's pre-built Seurat object (HT243_S1H4_seurat.RData -> all.rna, all.atac) the way
## BRCA/HT243B1-S1H4/bindSC/bindsc.R does, instead of Read10X_h5. Outputs 4 h5ads:
##   maxfuse_rna_shared.h5ad  : RNA lognorm on shared genes   (rows = <bc>_rna)
##   maxfuse_atac_shared.h5ad : ATAC gene-activity, same genes (rows = <bc>_atac)
##   maxfuse_rna_active.h5ad  : RNA scaled HVG                 (rows = <bc>_rna)
##   maxfuse_atac_active.h5ad : ATAC LSI                       (rows = <bc>_atac)
## Run in seurat4 env (Seurat/Signac + reticulate anndata for write_h5ad). See 00_prep_maxfuse_input_sub.sh.

library(reticulate); use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat); library(Signac); library(Matrix)
library(EnsDb.Hsapiens.v86); library(anndata)

## ---- inputs: BRCA's pre-built Seurat object (same source as bindsc.R) ----
DATA <- "/path/to/multiomeBench/BRCA/HT243B1-S1H4"
load(file.path(DATA, "HT243_S1H4_seurat.RData"))   # -> all.rna, all.atac (ATAC assay named "ATAC")
OUT  <- file.path(DATA, "maxfuse", "input")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## test cells are a modality split of the same cells -> suffix so the joint latent matches the benchmark
all.rna  <- RenameCells(all.rna,  new.names = paste0(Cells(all.rna),  "_rna"))
all.atac <- RenameCells(all.atac, new.names = paste0(Cells(all.atac), "_atac"))

## RNA: lognorm + HVG + scale + PCA
all.rna <- NormalizeData(all.rna); all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna);     all.rna <- RunPCA(all.rna)

## ATAC: annotations + LSI (active) + GeneActivity on the RNA HVGs (shared)
ann <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(ann) <- paste0("chr", seqlevels(ann)); genome(ann) <- "hg38"
Annotation(all.atac) <- ann
all.atac <- RunTFIDF(all.atac); all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)
ga <- GeneActivity(all.atac, features = VariableFeatures(all.rna))   # default assay = "ATAC" ChromatinAssay
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = ga)
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- NormalizeData(all.atac)

## shared gene set present in BOTH modalities
gene.use <- intersect(VariableFeatures(all.rna), rownames(all.atac[["ACTIVITY"]]))
hvg      <- VariableFeatures(all.rna)
cat("shared genes:", length(gene.use), " | RNA HVG:", length(hvg), "\n")

## cells already carry the _rna/_atac suffix (RenameCells), so matrix colnames -> h5ad rownames
save_h5ad <- function(mat, path) {                 # mat: features x cells -> cells x features
  m <- t(as.matrix(mat))                           # rownames(m) come from colnames(mat) = suffixed barcodes
  write_h5ad(AnnData(X = m), path)
}
# Seurat v5 (Assay5): layers via [["assay"]]$layer (matches bindsc.R's all.rna$RNA$data etc.)
save_h5ad(all.rna[["RNA"]]$data[gene.use, ],       file.path(OUT, "maxfuse_rna_shared.h5ad"))
save_h5ad(all.atac[["ACTIVITY"]]$data[gene.use, ], file.path(OUT, "maxfuse_atac_shared.h5ad"))
save_h5ad(all.rna[["RNA"]]$scale.data[hvg, ],      file.path(OUT, "maxfuse_rna_active.h5ad"))

al <- Embeddings(all.atac, reduction = "lsi")      # cells x LSI (active); rownames already suffixed
write_h5ad(AnnData(X = as.matrix(al)), file.path(OUT, "maxfuse_atac_active.h5ad"))

cat("wrote 4 MaxFuse input h5ads to", OUT, "\n")
