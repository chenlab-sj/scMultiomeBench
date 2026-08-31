# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyr)

source("/path/to/tools/scJoint/data_to_h5.R")
data.dir <- '/path/to/data/BMMC_d1/'
load(paste0(data.dir,"all_seurat.RData" ))
rna3 <- RenameCells(rna3, new.names = paste0(Cells(rna3),"_rna3"))
all.rna <- rna3
all.rna <-JoinLayers(all.rna)
all.rna <- NormalizeData(all.rna)

atac2 <- RenameCells(atac2, new.names = paste0(Cells(atac2),"_atac2"))
all.atac<- atac2
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
print(length(common.genes))
acces.counts<-as.matrix(acces.counts)
acces.counts.sel <-acces.counts[common.genes,]

outdir <- "/path/to/multiomeBench/BMMC_d1/crosssite/scjoint/"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
setwd(outdir)
writeLines(as.character(length(common.genes)), "input_size.txt")   # for 02_scjoint_patch_sub.sh (config input_size)

acces_counts.lst<-list()
acces_counts.lst[[1]]<-acces.counts.sel
write_h5_scJoint(acces_counts.lst,"atac_gene_scjoint.h5")

rna.counts<-as.matrix(rna_count)
rna.counts.sel<-rna.counts[common.genes,] 
rna_counts.lst<-list()
rna_counts.lst[[1]]<-rna.counts.sel
write_h5_scJoint(rna_counts.lst,"rna_scjoint.h5")

label_file<-"/path/to/data/BMMC_d1/celltype.csv"
label_df <- read.csv(label_file,row.names = 1)
rna_celltype <- label_df


test_celltype<- rna_celltype[colnames(rna.counts.sel),]$cell_type
print(unique(test_celltype))
test_celltype <- ifelse(is.na(test_celltype), "unknown",test_celltype)
test_celltype.lst<-list()
test_celltype.lst[[1]]<-test_celltype
write_csv_scJoint(test_celltype.lst,"testrna_celltype_scjoint.csv")

rna_bc.file = paste0(outdir,"rna_bc.txt")
writeLines(Cells(all.rna), rna_bc.file)
atac_bc.file = paste0(outdir,"atac_bc.txt")
writeLines(Cells(all.atac), atac_bc.file)
unique(test_celltype)
