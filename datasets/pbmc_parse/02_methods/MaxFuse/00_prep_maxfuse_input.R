# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## 00_prep_maxfuse_input.R  (pbmc_parse cross-platform)
## Build MaxFuse inputs using the SAME Signac recipe as the pbmc3k maxfuse prep:
##   RNA  : lognorm expression (shared genes) + scaled HVG (active)
##   ATAC : GeneActivity scores (shared genes) + LSI embedding (active)
##
## CROSS-PLATFORM DIFFERENCE vs pbmc3k: there RNA and ATAC were a modality split
## of the SAME multiome cells (built from one 10x h5 + fragments). HERE the RNA is
## independent Parse cells (1481) and the ATAC is pbmc3k cells (1642) -- unpaired,
## different barcodes. MaxFuse is a diagonal/unpaired method, so this is exactly what
## it expects; no pairing is used. We therefore load the two pre-built cross-platform
## Seurat objects from pbmc_parse_seurat.RData (the same RData the pbmc_parse bindSC
## adaptation uses): `all.rna` (Parse, RNA assay) and `all.atac` (pbmc3k, `peaks`
## assay + fragments-backed annotation so GeneActivity works). RNA cells get '_rna',
## ATAC cells get '_atac' so the joint latent.csv matches the other pbmc_parse methods.
## Run in the existing `seurat4` conda env (has Seurat/Signac + reticulate->anndata).

# seurat4 env (like 00_scBridge_s0.R / 00_scjoint_s1.R): its python has the `anndata`
# module that reticulate needs for write_h5ad. Without this binding, write_h5ad fails.
library(reticulate); use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat); library(Signac); library(Matrix)
library(EnsDb.Hsapiens.v86); library(anndata)

## ---- inputs: cross-platform Seurat objects from input_azimuth ----
DATA <- "/path/to/multiomeBench/pbmc_parse/input_azimuth"
OUT  <- "/path/to/multiomeBench/pbmc_parse/MaxFuse/input"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## pbmc_parse_seurat.RData provides `all.rna` (Parse RNA) and `all.atac` (pbmc3k,
## `peaks` assay + fragments annotation). Same object set the pbmc_parse bindSC uses.
load(file.path(DATA, "pbmc_parse_seurat.RData"))

## RNA: lognorm + HVG + scale + PCA  (identical to pbmc3k maxfuse prep)
all.rna <- NormalizeData(all.rna); all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna);     all.rna <- RunPCA(all.rna)

## ATAC: ensure annotation (cross-platform RData may not carry it) then LSI (active) +
## GeneActivity on the RNA HVGs (shared). Identical recipe to pbmc3k maxfuse prep.
DefaultAssay(all.atac) <- "peaks"
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
cat("shared genes:", length(gene.use), " | RNA HVG:", length(hvg), "\n")

## barcodes: cross-platform cells, suffixed by modality (RNA cells != ATAC cells)
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
