# NOTE: paths below are placeholders. See config/config.R and the README.
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr) 
library(Matrix)
library(anndata)
library(tidyr)



data.dir <- "/path/to/data/HTAN/HT137B1-S1H7/"
load(paste0(data.dir,"HT137B1-S1H7_seurat.RData"))
n_lat <-30

all.rna <- RenameCells(all.rna, new.names = paste0(Cells(all.rna),"_rna")) 
all.rna <- NormalizeData(all.rna)

all.atac <-RenameCells(all.atac,  new.names = paste0(Cells(all.atac),"_atac"))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations
gene.activities <- GeneActivity(all.atac,features = VariableFeatures(all.rna))
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")
rna_count <- all.rna$RNA$counts
common.genes <- intersect(row.names(acces.counts),row.names(rna_count))
acces.counts<-as.matrix(acces.counts)
acces.counts.sel <-acces.counts[common.genes,]
acces.counts.sel <- as(acces.counts.sel, "sparseMatrix")

atac.adata <- AnnData(X= t(acces.counts.sel))
atac.adata$var$gene_ids <- common.genes
outdir <- paste0(data.dir,"scBridge/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
write_h5ad(atac.adata, paste0(outdir,"scBridge_atac.h5ad"))


rna.counts<-as.matrix(rna_count)
rna.counts.sel<-rna.counts[common.genes,] 
rna.counts.sel <- as(rna.counts.sel, "sparseMatrix")

rna_celltype.file <- paste0(data.dir,"rna_celltype.csv")
rna_celltype <- read.csv(rna_celltype.file,row.names = 1)
rna.adata <- AnnData(X= t(rna.counts.sel))
test_celltype<- rna_celltype[colnames(rna.counts.sel),]$cell_type
test_celltype <- ifelse(is.na(test_celltype), "unknown",test_celltype)

rna.adata <- AnnData(X= t(rna.counts.sel))
rna.adata$obs$CellType <-test_celltype
rna.adata$var$gene_ids <- common.genes
write_h5ad(rna.adata, paste0(outdir,"scBridge_rna.h5ad"))
