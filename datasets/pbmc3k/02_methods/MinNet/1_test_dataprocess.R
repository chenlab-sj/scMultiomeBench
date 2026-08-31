# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(ggplot2)
library(patchwork)
library(rhdf5)
library(EnsDb.Hsapiens.v86)
library(Signac)
library(Matrix)
rm(list=ls())

set.seed(14)

## according the pre-process from author

bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
out.dir <-"/path/to/tools/MinNet/pbmc3k/"

inputdata.10x<-Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks



# Create RNA Seurat object ----
rna <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = "RNA")

rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))

rna <- FindNeighbors(rna, dims = 1:15)
rna <- FindClusters(rna, resolution = 0.1)
rna <- RunUMAP(rna, dims = 1:15)

#DimPlot(rna, reduction = "umap",pt.size = 0.5)



# Create ATAC Seurat object ----
#                                            seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments= fpath
)
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# convert to UCSC style
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"
Annotation(atac) <- annotation

gene.activities <- GeneActivity(atac)
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(atac) <- "ACTIVITY"

atac <- FindVariableFeatures(atac)
atac <- NormalizeData(atac)
atac <- ScaleData(atac)


DefaultAssay(atac) <- "ATAC"
VariableFeatures(atac) <- names(which(Matrix::rowSums(atac) > 50))
#atac <- RunSLSI(atac, n = 50, scale.max = NULL)
#atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30)
#DimPlot(atac, reduction = "umap") #,group.by = 'predicted.id')

#saveRDS(atac,'atac.rds')

# save for SiaNN ----
#rna <- readRDS('rna.rds')
rna <- RenameCells(rna,  new.names = paste0(Cells(rna),"_rna")) 
rna_mtx <- t(rna$RNA$counts)
rna_mtx <- as.matrix(rna_mtx)
h5createFile(paste0(out.dir,"rna.test.h5"))
h5write(rna_mtx, paste0(out.dir,"rna.test.h5"), 'RNA')
rna_meta <- rna@meta.data
write.table(rna_meta,paste0(out.dir, 'rna.meta.csv'),quote = F, sep='\t')
rna_gene <- rownames(rna)
write.table(rna_gene, paste0(out.dir,'rna_gene_name.txt'),quote=F,row.names=F,col.names=F)

#atac <- readRDS('atac.rds')
atac <- RenameCells(atac,  new.names = paste0(Cells(atac),"_atac")) 
atac_mtx <- Matrix::t(atac$ACTIVITY$counts)
atac_mtx <- as.matrix(atac_mtx)
h5createFile(paste0(out.dir,'atac.test.h5'))
h5write(atac_mtx, paste0(out.dir,'atac.test.h5'), 'ATAC')
atac_meta <- atac@meta.data
write.table(atac_meta, paste0(out.dir,'atac.meta.csv'),quote = F, sep='\t')
DefaultAssay(atac) <- 'ACTIVITY'

atac_gene <- rownames(atac)
write.table(atac_gene, paste0(out.dir,'atac_gene_name.txt'),quote=F,row.names=F,col.names=F)




















