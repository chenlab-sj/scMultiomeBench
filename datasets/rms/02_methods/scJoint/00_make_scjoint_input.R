# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v75)
library(tidyr)
data.dir <- "/path/to/data/RMS/"

bcmtx.h5 <- "/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5"
fpath <-'/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/atac_fragments.tsv.gz'

source("/path/to/tools/scJoint/data_to_h5.R")
source("/path/to/tools/scJoint/data_to_h5.R")
inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_count  <- inputdata.10x$`Gene Expression`
atac_count  <- inputdata.10x$Peaks

chrom_assay <- CreateChromatinAssay(
  counts = atac_count,
  sep = c(":", "-"),
  fragments= fpath
)
all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# convert to UCSC style
#seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg19"
Annotation(all.atac) <- annotation

# create RNA seurat object
all.rna <- CreateSeuratObject(counts = rna_count,assay = "RNA")
all.rna <- NormalizeData(all.rna)

gene.activities <- GeneActivity(all.atac)
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")

common.genes <- intersect(row.names(acces.counts),row.names(rna_count))
acces.counts<-as.matrix(acces.counts)
acces.counts.sel <-acces.counts[common.genes,]
colnames(acces.counts.sel)<-paste0(colnames(acces.counts.sel),"_atac")

acces_counts.lst<-list()
acces_counts.lst[[1]]<-acces.counts.sel
write_h5_scJoint(acces_counts.lst,"atac_gene_scjoint.h5")

## common genes between rna and atac
rna.counts<-as.matrix(rna_count)
rna.counts.sel<-rna.counts[common.genes,] 
rna_counts.lst<-list()
colnames(rna.counts.sel)<- paste0(colnames(rna.counts.sel),"_rna")
rna_counts.lst[[1]]<-rna.counts.sel
write_h5_scJoint(rna_counts.lst,"rna_scjoint.h5")

########################################################################

##### cell type label
########################################################################
## celltype list file

label_file<-"//path/to/data/RMS/Mast607A_TB19_22652/Mast607A_SJRHB013758_X2_singleR.csv"
label_df <- read.csv(label_file,row.names = 1)
row.names(label_df)<- paste0(row.names(label_df),"_rna")
test_label <- label_df[colnames(rna.counts.sel),]
test_celltype.lst<-list()
test_celltype.lst[[1]]<-test_label$labels
# csv_list<-list()
# csv_list[[1]]<-"/path/to/tools/scJoint/RMS/SJRHB10468_X1/testrna_celltype_scjoint.csv"
write_csv_scJoint(test_celltype.lst,"testrna_celltype_scjoint.csv")
#write.csv(test_celltype.lst[[1]],"testrna_celltype_scjoint.csv")
out.dir <- "/path/to/tools/scJoint/RMS/Mast607A_TB19_22652/"
rna_bc.file = paste0(out.dir,"rna_bc.txt")
writeLines(colnames(rna_count), rna_bc.file)
atac_bc.file = paste0(out.dir,"atac_bc.txt")
writeLines(colnames(atac_count), atac_bc.file)
unique(test_label$labels)

