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
rna1 <- RenameCells(rna1, new.names = paste0(Cells(rna1),"_rna1")) 
rna2 <- RenameCells(rna2, new.names = paste0(Cells(rna2),"_rna2")) 
rna3 <- RenameCells(rna3, new.names = paste0(Cells(rna3),"_rna3")) 

rna1 <- NormalizeData(rna1)
rna2 <- NormalizeData(rna2)
rna3 <- NormalizeData(rna3)

atac1 <- RenameCells(atac1 , new.names = paste0(Cells(atac1),"_atac1")) 
atac2 <- RenameCells(atac2, new.names = paste0(Cells(atac2),"_atac2")) 
atac3 <- RenameCells(atac3, new.names = paste0(Cells(atac3),"_atac3"))

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(atac1) <- annotations
Annotation(atac2) <- annotations
Annotation(atac3) <- annotations

gene.activities1 <- GeneActivity(atac1)
atac1[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities1)
DefaultAssay(atac1) <- "ACTIVITY"
acces.counts1 <- GetAssayData(object = atac1, layer  = "ACTIVITY", slot = "counts")
rna_count1 <- rna1$RNA$counts
common.genes1 <- intersect(row.names(acces.counts1),row.names(rna_count1))


gene.activities2 <- GeneActivity(atac2)
atac2[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities2)
DefaultAssay(atac2) <- "ACTIVITY"
acces.counts2 <- GetAssayData(object = atac2, layer  = "ACTIVITY", slot = "counts")
rna_count2 <- rna2$RNA$counts
common.genes2 <- intersect(row.names(acces.counts2),row.names(rna_count2))

gene.activities3 <- GeneActivity(atac3)
atac3[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities3)
DefaultAssay(atac3) <- "ACTIVITY"
acces.counts3 <- GetAssayData(object = atac3, layer  = "ACTIVITY", slot = "counts")
rna_count3 <- rna3$RNA$counts
common.genes3 <- intersect(row.names(acces.counts3),row.names(rna_count3))

common.genes <- intersect(common.genes1,common.genes2)
common.genes<- intersect(common.genes,common.genes3)
print(length(common.genes))
acces.counts1<-as.matrix(acces.counts1)
acces.counts1.sel <-acces.counts1[common.genes,]
acces.counts2<-as.matrix(acces.counts2)
acces.counts2.sel <-acces.counts2[common.genes,]
acces.counts3<-as.matrix(acces.counts3)
acces.counts3.sel <-acces.counts3[common.genes,]

acces_counts.lst<-list()
acces_counts.lst[[1]]<-acces.counts1.sel
acces_counts.lst[[2]]<-acces.counts2.sel
acces_counts.lst[[3]]<-acces.counts3.sel
write_h5_scJoint(acces_counts.lst,c("atac1_gene_scjoint.h5",
                                    "atac2_gene_scjoint.h5",
                                    "atac3_gene_scjoint.h5"))

rna.counts1<-as.matrix(rna_count1)
rna.counts1.sel<-rna.counts1[common.genes,]
rna.counts2<-as.matrix(rna_count2)
rna.counts2.sel<-rna.counts2[common.genes,]
rna.counts3<-as.matrix(rna_count3)
rna.counts3.sel<-rna.counts3[common.genes,]


rna_counts.lst<-list()
rna_counts.lst[[1]]<-rna.counts1.sel
rna_counts.lst[[2]]<-rna.counts2.sel
rna_counts.lst[[3]]<-rna.counts3.sel
write_h5_scJoint(rna_counts.lst,c("rna1_scjoint.h5",
                                  "rna2_scjoint.h5",
                                  "rna3_scjoint.h5"))

label_file<-"/path/to/data/BMMC_d1/celltype.csv"
label_df <- read.csv(label_file,row.names = 1)
rna_celltype <- label_df


test_celltype1<- rna_celltype[colnames(rna.counts1.sel),]$cell_type
print(unique(test_celltype1))
test_celltype1 <- ifelse(is.na(test_celltype1), "unknown",test_celltype1)

test_celltype2<- rna_celltype[colnames(rna.counts2.sel),]$cell_type
print(unique(test_celltype2))
test_celltype2 <- ifelse(is.na(test_celltype2), "unknown",test_celltype2)

test_celltype3<- rna_celltype[colnames(rna.counts3.sel),]$cell_type
print(unique(test_celltype3))
test_celltype3 <- ifelse(is.na(test_celltype3), "unknown",test_celltype3)

test_celltype.lst<-list()
test_celltype.lst[[1]]<-test_celltype1
test_celltype.lst[[2]]<-test_celltype2
test_celltype.lst[[3]]<-test_celltype3
write_csv_scJoint(test_celltype.lst,c("test1rna_celltype_scjoint.csv",
                                      "test2rna_celltype_scjoint.csv",
                                      "test3rna_celltype_scjoint.csv"))

outdir <- paste0(data.dir,"scjoint2/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
rna1_bc.file = paste0(outdir,"rna1_bc.txt")
writeLines(Cells(rna1), rna1_bc.file)
atac1_bc.file = paste0(outdir,"atac1_bc.txt")
writeLines(Cells(atac1), atac1_bc.file)

rna2_bc.file = paste0(outdir,"rna2_bc.txt")
writeLines(Cells(rna2), rna2_bc.file)
atac2_bc.file = paste0(outdir,"atac2_bc.txt")
writeLines(Cells(atac2), atac2_bc.file)

rna3_bc.file = paste0(outdir,"rna3_bc.txt")
writeLines(Cells(rna3), rna3_bc.file)
atac3_bc.file = paste0(outdir,"atac3_bc.txt")
writeLines(Cells(atac3), atac3_bc.file)
unique(test_celltype)

