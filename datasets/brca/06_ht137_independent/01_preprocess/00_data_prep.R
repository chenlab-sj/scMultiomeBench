# NOTE: paths below are placeholders. See config/config.R and the README.
## data download
## https://humantumoratlas.org/explore?selectedFilters=%5B%7B%22group%22%3A%22AtlasName%22%2C%22value%22%3A%22HTAN+WUSTL%22%7D%5D
## search by piece_ID
# module load conda3/202311
#conda activate myenv
## eg.synapse get syn52176727

## generate tbi for fragement tsv.gz
# https://github.com/stuart-lab/signac/issues/242
# gzip -d <fragments>
#   bgzip <fragments>
#   tabix -p bed <fragments>
data.dir <- "/path/to/data/HTAN/HT137B1-S1H7/"
rna_data <- readRDS(paste0(data.dir,"syn53214894/HT137B1-S1H7.rds"))
atac_data <- readRDS(paste0(data.dir,"syn53215781/HT137B1-S1H7.rds"))
fpath <- paste0(data.dir,"HT137B1-S1H7-fragments.tsv.gz")
rna_count <- rna_data$RNA$counts
colnames(rna_count)<-gsub(".*_", "", colnames(rna_count))
atac_count <- atac_data$peaksMACS2$counts
colnames(atac_count)<-gsub(".*_", "", colnames(atac_count))
rna_count.filter <- rna_count[which(rowSums(rna_count)!=0),]
atac_count.filter <- atac_count[which(rowSums(atac_count)!=0),]
#############################################################
## prepare for R
all.rna <- CreateSeuratObject(counts = rna_count.filter,assay = "RNA", project = "scRNA")

chrom_assay <- CreateChromatinAssay(
  counts = atac_count.filter,
  sep = c(":", "-"),
  fragments= fpath
)

all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
)

save(all.rna, all.atac, file = paste0(data.dir,"HT137B1-S1H7_seurat.RData"))

#############################################################
## save celltype
rna_celltype <- data.frame(
  bc =paste0(rna_data$Original_barcode,"_rna"),
  cell_type =rna_data$cell_type
)

rna_celltype$modality = "test_scRNA"

atac_celltype <- data.frame(
  bc =rownames(atac_data@meta.data),
  cell_type =atac_data$cell_type
)

atac_celltype$modality = "test_scATAC"
atac_celltype$bc <- gsub(".*_", "", atac_celltype$bc)
atac_celltype$bc <- paste0(atac_celltype$bc,"_atac")
write.csv(rna_celltype, paste0(data.dir,"rna_celltype.csv"),row.names = FALSE)
write.csv(atac_celltype, paste0(data.dir,"atac_celltype.csv"),row.names = FALSE)
#################################################################

## write to 10x
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)

rna<- AnnData(
  X = t(as.matrix(rna_count.filter)),
)

atac<- AnnData(
  X = t(as.matrix(atac_count.filter)),
)
write_h5ad(rna, paste0(data.dir,"HT137B1-S1H7_rna.h5ad"))
write_h5ad(atac, paste0(data.dir,"HT137B1-S1H7_atac.h5ad"))

library(EnsDb.Hsapiens.v86)
library(Signac)
library(Seurat)


data.dir <- "/path/to/data/HTAN/HT137B1-S1H7/"
load(paste0(data.dir,"HT137B1-S1H7_seurat.RData"))
atac_counts  <- all.atac$ATAC$counts

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments= fpath
)
all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# convert to UCSC style
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"
Annotation(all.atac) <- annotation
gene.activities <- GeneActivity(all.atac)
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")

atac_gene<- AnnData(
  X = t(as.matrix(acces.counts)),
)
write_h5ad(atac_gene, paste0(data.dir,"HT137B1-S1H7_atac_gene.h5ad"))


########################################################################
atac_counts  <- all.atac$ATAC$counts

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments= fpath
)
all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# convert to UCSC style
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations
## normalization 
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)


# Perform standard analysis of each modality independently RNA analysis
all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna)

rna.embed <- Embeddings(all.rna, reduction = "pca")
rownames(rna.embed)<- paste0(rownames(rna.embed),"_rna")

atac.embed <- Embeddings(all.atac, reduction = "lsi")
rownames(atac.embed)<- paste0(rownames(atac.embed),"_atac")


rna.embed<- AnnData(
  X = as.matrix(rna.embed),
)

atac.embed<- AnnData(
  X =as.matrix(atac.embed),
)

write_h5ad(rna.embed, paste0(data.dir,"rna_embed.h5ad"))
write_h5ad(atac.embed, paste0(data.dir,"atac_embed.h5ad"))




#######################################################################################


library(Seurat)
library(Signac)
data.dir <- "/path/to/data/HTAN/HT235B1-S1H1/"
rna_data <- readRDS(paste0(data.dir,"syn53214703_snRNA/HT235B1-S1H1.rds"))
atac_data <- readRDS(paste0(data.dir,"syn53215775/HT235B1-S1H1.rds"))

# module load tabix/0.2.6 
# tabix -p bed HT235B1-S1H1-atac_fragments.tsv.gz
fpath <- paste0(data.dir,"syn53215775/HT235B1-S1H1-atac_fragments.tsv.gz")

rna_data<- RenameCells(rna_data, new.names = gsub(".*_", "", Cells(rna_data))) 
atac_data<- RenameCells(atac_data, new.names = gsub(".*_", "", Cells(atac_data))) 
atac_data$Original_barcode <- Cells(atac_data)
sel.bc <- intersect(Cells(rna_data),Cells(atac_data))
rna_data.sel <- subset(rna_data, subset = Original_barcode %in% sel.bc)
atac_data.sel <- subset(atac_data, subset = Original_barcode %in% sel.bc)

rna_count.sel <- rna_data.sel$RNA$counts
atac_count.sel <- atac_data.sel$peaksMACS2$counts


#################################################################
## for atac, find common peak with test data

train.peaks <- StringToGRanges(regions = rownames(atac_count.sel), sep = c(":","-"))
#######################################################################################


## load test data
# test.dir <- "/path/to/data/HTAN/HT163B1-S1H6/"
# test.fpath <- paste0(test.dir,"syn53215756/HT163B1-S1H6-fragments.tsv.gz")
# load(paste0(test.dir,"S1H6_seurat.RData"))


test.dir <- "/path/to/data/HTAN/HT137B1-S1H7/"
test.fpath <- paste0(test.dir,"HT137B1-S1H7-fragments.tsv.gz")
load(paste0(test.dir,"HT137B1-S1H7_seurat.RData"))

test_atac_count <- all.atac$ATAC$counts
test.peaks <- StringToGRanges(regions = rownames(test_atac_count), sep = c(":","-"))
combined.peaks <- reduce(x = c(train.peaks,test.peaks))
## creat fragment obj

train.frags <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_count.sel)
)

test.frags <- CreateFragmentObject(
  path = test.fpath,
  cells = colnames(test_atac_count)
)

train.frag_counts <- FeatureMatrix(
  fragments = train.frags,
  features = combined.peaks,
  cells = colnames(atac_count.sel)
)

test.frag_counts <- FeatureMatrix(
  fragments = test.frags,
  features = combined.peaks,
  cells = colnames(test_atac_count)
)

#######################################################################
source("/path/to/multiomeBench/common/pbmc10k/0_create10x_datafmt.R")
library(dplyr)
train_features.file <- paste0(data.dir,"syn53215775/HT235B1-S1H1-features.tsv")
train_feature <- read.csv(train_features.file, sep = "\t", header = FALSE)

train_rna_feature<-train_feature[which(train_feature$V2 %in% rownames(rna_count.sel)),]
train_peaks_feature<-  strsplit(rownames(train.frag_counts), "-")
train_peaks_feature<- as.data.frame(do.call(rbind, train_peaks_feature))
train_peaks_feature$peak <- paste0(train_peaks_feature$V1,":",train_peaks_feature$V2,"-",train_peaks_feature$V3)
train_peaks_feature$peak2 <- train_peaks_feature$peak
train_peaks_feature$type <-"Peaks"
train_peaks_feature <- train_peaks_feature %>%
  dplyr:: select(peak, peak2,type,V1, V2, V3)
library(data.table)
train_feature <- rbindlist(list (train_rna_feature, train_peaks_feature), use.names = FALSE)
write.table(train_feature, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
rownames(train.frag_counts)<-train_peaks_feature$peak

sub_bc_matrix(barcodes = colnames(train.frag_counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = rbind(rna_count.sel,train.frag_counts),
              output_dir <- "train_forHT137B1-S1H7/"
)

train<-sub_h5(barcodes = colnames(train.frag_counts),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count.sel,train.frag_counts),
              path = "train_forHT137B1-S1H7.h5",
              test = TRUE)

## for test side



write.table(train_peaks_feature, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
rownames(test.frag_counts) <- sub("-", ":", rownames(test.frag_counts))
sub_bc_matrix(barcodes = colnames(test.frag_counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = test.frag_counts,
              output_dir <- "HT137B1-S1H7_commonpeak/",
              module = 'Peaks')




test<-sub_h5(barcodes = colnames(test.frag_counts),
              feature_mtx = "tmp.tsv", 
              count_mtx = test.frag_counts,
              path = "HT137B1-S1H7_commonpeak.h5",
              test = TRUE,
              module = 'Peaks')


testrna_feature <-train_rna_feature[which(rowSums(rna_count)!=0),]
write.table(testrna_feature, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
sub_bc_matrix(barcodes = colnames(rna_count.filter),
              feature.mtx = "tmp.tsv", 
              count.mtx = rna_count.filter,
              output_dir = "HT137B1-S1H7_rna/",
              module = 'Gene Expression')


test<-sub_h5(barcodes = colnames(rna_count.filter),
             feature_mtx = "tmp.tsv", 
             count_mtx = rna_count.filter,
             path = "HT137B1-S1H7_rna.h5",
             test = TRUE,
             module = 'Gene Expression')