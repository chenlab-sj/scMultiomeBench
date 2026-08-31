# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(Matrix)
library(anndata)
library(tidyr)



data.dir <- '/path/to/data/BMMC_d1/'
load(paste0(data.dir,"all_seurat.RData" ))
rna3 <- RenameCells(rna3, new.names = paste0(Cells(rna3),"_rna3"))
all.rna <- rna3
all.rna <- NormalizeData(all.rna)
n_lat <-30

atac3 <- RenameCells(atac3, new.names = paste0(Cells(atac3),"_atac3"))
all.atac <- atac3
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations

all.rna <- FindVariableFeatures(all.rna)
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
outdir <- paste0("/path/to/multiomeBench/BMMC_d1/s1d1_paired/scBridge/")
if (!file.exists(outdir)) {
  dir.create(outdir, recursive=TRUE)
}
write_h5ad(atac.adata, paste0(outdir,"scBridge_atac.h5ad"))


rna.counts<-as.matrix(rna_count)
rna.counts.sel<-rna.counts[common.genes,] 
rna.counts.sel <- as(rna.counts.sel, "sparseMatrix")

label_file<-"/path/to/data/BMMC_d1/celltype.csv"
label_df <- read.csv(label_file,row.names = 1)
rna_celltype <- label_df


rna.adata <- AnnData(X= t(rna.counts.sel))

test_celltype<- rna_celltype[colnames(rna.counts.sel),]$cell_type
print(unique(test_celltype))
test_celltype <- ifelse(is.na(test_celltype), "unknown",test_celltype)

#rna.adata <- AnnData(X= t(rna.counts.sel))
rna.adata$obs$CellType <-test_celltype
rna.adata$var$gene_ids <- common.genes
write_h5ad(rna.adata, paste0(outdir,"scBridge_rna.h5ad"))
