# NOTE: paths below are placeholders. See config/config.R and the README.
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyr)

source("/path/to/tools/scJoint/data_to_h5.R")
data.dir <- "/path/to/data/HTAN/HT137B1-S1H7/"
load(paste0(data.dir,"HT137B1-S1H7_seurat.RData"))

all.rna <- RenameCells(all.rna, new.names = paste0(Cells(all.rna),"_rna")) 
all.rna <- NormalizeData(all.rna)

all.atac <-RenameCells(all.atac,  new.names = paste0(Cells(all.atac),"_atac"))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations

gene.activities <- GeneActivity(all.atac)
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"

acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")
rna_count <- all.rna$RNA$counts
common.genes <- intersect(row.names(acces.counts),row.names(rna_count))
acces.counts<-as.matrix(acces.counts)
acces.counts.sel <-acces.counts[common.genes,]

acces_counts.lst<-list()
acces_counts.lst[[1]]<-acces.counts.sel
write_h5_scJoint(acces_counts.lst,"atac_gene_scjoint.h5")

rna.counts<-as.matrix(rna_count)
rna.counts.sel<-rna.counts[common.genes,] 
rna_counts.lst<-list()
rna_counts.lst[[1]]<-rna.counts.sel
write_h5_scJoint(rna_counts.lst,"rna_scjoint.h5")

rna_celltype.file <- paste0(data.dir,"rna_celltype.csv")
rna_celltype <- read.csv(rna_celltype.file,row.names = 1)
test_celltype<- rna_celltype[colnames(rna.counts.sel),]$cell_type
test_celltype <- ifelse(is.na(test_celltype), "unknown",test_celltype)
test_celltype.lst<-list()
test_celltype.lst[[1]]<-test_celltype
write_csv_scJoint(test_celltype.lst,"testrna_celltype_scjoint.csv")

outdir <- paste0(data.dir,"scjoint/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
rna_bc.file = paste0(outdir,"rna_bc.txt")
writeLines(Cells(all.rna), rna_bc.file)
atac_bc.file = paste0(outdir,"atac_bc.txt")
writeLines(Cells(all.atac), atac_bc.file)
unique(test_celltype)
