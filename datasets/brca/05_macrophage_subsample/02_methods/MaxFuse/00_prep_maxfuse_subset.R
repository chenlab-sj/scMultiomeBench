# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## 00_prep_maxfuse_subset.R -- per-(MACROSUB, SUB) MaxFuse inputs for the HT243 macrophage-subsample analysis.
## SAME Signac recipe as HT243B1-S1H4/maxfuse/prep_maxfuse_input.R, but loads the per-sub Seurat object
## (HT243_S1H4_sub{SUB}.RData -> all.rna_sub, all.atac_sub, as the subset bindsc.R does) and writes the 4
## h5ads to {MACROSUB}/maxfuse/sub{SUB}/input/. Env: MACROSUB (HT243_S1H4_macrosub..5) + SUB (1..5).
## Run in seurat4 env (Seurat/Signac + reticulate anndata). Isolated per (macrosub,sub) -> parallel-safe.

library(reticulate); use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat); library(Signac); library(Matrix)
library(EnsDb.Hsapiens.v86); library(anndata)

SUBS <- "/path/to/multiomeBench/BRCA/HT243-S1H4_subsample"
MACROSUB <- Sys.getenv("MACROSUB")
SUB      <- Sys.getenv("SUB")
stopifnot(nzchar(MACROSUB), nzchar(SUB))

## ---- inputs: the per-sub Seurat object (same objects the subset bindsc.R uses) ----
load(file.path(SUBS, MACROSUB, paste0("HT243_S1H4_sub", SUB, ".RData")))   # -> all.rna_sub, all.atac_sub
all.rna  <- all.rna_sub
all.atac <- all.atac_sub
OUT <- file.path(SUBS, MACROSUB, "maxfuse", paste0("sub", SUB), "input")
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
ga <- GeneActivity(all.atac, features = VariableFeatures(all.rna))
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = ga)
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- NormalizeData(all.atac)

## shared gene set present in BOTH modalities
gene.use <- intersect(VariableFeatures(all.rna), rownames(all.atac[["ACTIVITY"]]))
hvg      <- VariableFeatures(all.rna)
cat(MACROSUB, "sub", SUB, "shared genes:", length(gene.use), " | RNA HVG:", length(hvg), "\n")

save_h5ad <- function(mat, path) {
  m <- t(as.matrix(mat))
  write_h5ad(AnnData(X = m), path)
}
save_h5ad(all.rna[["RNA"]]$data[gene.use, ],       file.path(OUT, "maxfuse_rna_shared.h5ad"))
save_h5ad(all.atac[["ACTIVITY"]]$data[gene.use, ], file.path(OUT, "maxfuse_atac_shared.h5ad"))
save_h5ad(all.rna[["RNA"]]$scale.data[hvg, ],      file.path(OUT, "maxfuse_rna_active.h5ad"))

al <- Embeddings(all.atac, reduction = "lsi")
write_h5ad(AnnData(X = as.matrix(al)), file.path(OUT, "maxfuse_atac_active.h5ad"))

cat("wrote 4 MaxFuse input h5ads to", OUT, "\n")
