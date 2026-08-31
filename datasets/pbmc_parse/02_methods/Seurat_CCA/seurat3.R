# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyr)
t1 <- Sys.time()
data.dir <- "/path/to/multiomeBench/pbmc_parse/input_azimuth/"
load(paste0(data.dir,"pbmc_parse_seurat.RData"))
n_lat <-30
all.rna <- RenameCells(all.rna,  new.names = paste0(Cells(all.rna),"_rna")) 
all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna)
all.rna <- RunUMAP(all.rna, dims = 1:n_lat)


all.atac <-RenameCells(all.atac,  new.names = paste0(Cells(all.atac),"_atac"))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)
all.atac <- RunUMAP(all.atac, reduction = "lsi", dims = 2:n_lat, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
gene.activities <- GeneActivity(all.atac,features = VariableFeatures(all.rna))
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- NormalizeData(all.atac)
all.atac <- ScaleData(all.atac, features = rownames(all.atac))

#############################################################
## integration
#############################################################

transfer.anchors <- FindTransferAnchors(reference = all.rna, query = all.atac, features = VariableFeatures(object = all.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
genes.use <- VariableFeatures(all.rna)
refdata <- GetAssayData(all.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = all.atac[["lsi"]],
                           dims = 2:n_lat)
all.atac$RNA <- imputation
coembed <- merge(x = all.rna, y = all.atac)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
t2 <- Sys.time()
df_umap = as.data.frame(coembed@reductions$pca@cell.embeddings[,1:n_lat])
colnames(df_umap) = paste0("latent_",1:ncol(df_umap))


outdir <- '/path/to/multiomeBench/pbmc_parse/Seurat3/'


write.csv(df_umap,paste0(outdir,"latent.csv"))

write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)


