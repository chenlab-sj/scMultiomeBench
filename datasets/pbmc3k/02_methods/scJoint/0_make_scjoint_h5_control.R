# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v75)
library(tidyr)
data.dir <- "/path/to/data/RMS/"
rna_count.mtx <- paste0(data.dir,"SJRHB010468_X1_scRNA/GSM5390487_1763612_DYE2651_matrix.mtx")
rna_bc.tsv <- paste0(data.dir,"SJRHB010468_X1_scRNA/GSM5390487_1763612_DYE2651_barcodes.tsv")
rna.feature.tsv <-paste0(data.dir,"SJRHB010468_X1_scRNA/GSM5390487_1763612_DYE2651_features.tsv")
atac_count.mtx <- paste0(data.dir, "SJRHB010468_X1_scATAC/GSM5390508_1761672_DYE2649_matrix.mtx")
atac_peak.bed <- paste0(data.dir,"SJRHB010468_X1_scATAC/GSM5390508_1761672_DYE2649_peaks.bed")
atac_bc.tsv <- paste0(data.dir,"SJRHB010468_X1_scATAC/GSM5390508_1761672_DYE2649_barcodes.tsv")
fpath <- paste0(data.dir,"SJRHB010468_X1_scATAC/GSM5390508_1761672_DYE2649_fragments.tsv.gz")

source("/path/to/tools/scJoint/data_to_h5.R")

rna_count <-readMM(rna_count.mtx)
rna_bc <- read.table(rna_bc.tsv, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rna_feature <- read.table(rna.feature.tsv, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
rna_feature <-separate(rna_feature, V2, into = c("ref","gene"), sep = "_", remove = FALSE)
rna_count <- rna_count[rna_feature$ref=="hg19",]
rna_hg19gene <- rna_feature[rna_feature$ref=="hg19",]$gene
row.names(rna_count)<- rna_hg19gene
colnames(rna_count)<-rna_bc$V1
rna_count <- rna_count[!duplicated(rna_hg19gene),]


atac_count <- readMM(atac_count.mtx)
atac_feature <- read.table(atac_peak.bed , header = FALSE, sep = "\t", stringsAsFactors = FALSE)
atac_bc <- read.table(atac_bc.tsv, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
atac_feature$peak <- paste0(atac_feature$V1,":",atac_feature$V2,"-",atac_feature$V3)
rownames(atac_count)<-atac_feature$peak
colnames(atac_count)<-atac_bc$V1

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
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
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
label_file<-"/path/to/data/RMS/SJRHB010468_X1_meta.csv"
label_df <- read.csv(label_file,row.names = 1)
test_label <- label_df[colnames(rna.counts.sel),]
test_celltype.lst<-list()
test_celltype.lst[[1]]<-test_label$cell_type
csv_list<-list()
csv_list[[1]]<-"/path/to/tools/scJoint/RMS/SJRHB10468_X1/testrna_celltype_scjoint.csv"
write_csv_scJoint(test_celltype.lst,"testrna_celltype_scjoint.csv")
write.csv(test_celltype.lst[[1]],"testrna_celltype_scjoint.csv")
out.dir <- "/path/to/tools/scJoint/RMS/SJRHB10468_X1/"
rna_bc.file = paste0(out.dir,"rna_bc.txt")
writeLines(colnames(rna_count), rna_bc.file)
atac_bc.file = paste0(out.dir,"atac_bc.txt")
writeLines(colnames(atac_count), atac_bc.file)
unique(test_label$cell_type)

