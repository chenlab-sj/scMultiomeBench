# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## HT243_S1H4_adj
## load HT243_S1H4_celltype propotion

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
data.dir <- "/path/to/data/HTAN/HT243B1-S1H4/"
data2.dir <- "/path/to/data/HTAN/HT243_S1H4_adj/"
fpath <- fpath <- paste0(data.dir,"HT243B1-S1H4-atac_fragments.tsv.gz")
celltype <- read.csv(paste0(train.dir,'muti_celltype.csv'))
tumor_rm <- celltype %>%
  filter(rna_celltype == "Tumor")%>%
  sample_frac(0.7)
bc_rm <- tumor_rm$bc

load(paste0(data.dir,"HT243_S1H4_seurat.RData"))
rna_count<-all.rna$RNA$counts
rna_count.sel <- rna_count[,!colnames(rna_count) %in% bc_rm]
atac_count<-all.atac$ATAC$counts
atac_count.sel <- atac_count[,!colnames(atac_count) %in% bc_rm]

#############################################################
## prepare for R
all.rna <- CreateSeuratObject(counts = rna_count.sel,assay = "RNA", project = "scRNA")

chrom_assay <- CreateChromatinAssay(
  counts = atac_count.sel,
  sep = c(":", "-"),
  fragments= fpath
)

all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
)

save(all.rna, all.atac, file = paste0(data2.dir,"HT243_S1H4_adj_seurat.RData"))

#############################################################
## write to 10x
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)

rna<- AnnData(
  X = t(as.matrix(rna_count.sel)),
)

atac<- AnnData(
  X = t(as.matrix(atac_count.sel)),
)
write_h5ad(rna, paste0(data2.dir,"HT243_S1H4adj_rna.h5ad"))
write_h5ad(atac, paste0(data2.dir,"HT243_S1H4adj_atac.h5ad"))
celltype_adj <- celltype %>%
  dplyr::filter(! bc %in% bc_rm )

write.csv(celltype_adj, paste0(data2.dir,"muti_celltype.csv"),row.names = FALSE)


############################################
## from training data
newtest_dir <-'/path/to/data/HTAN/HT263B1-S1H1_train_HT243_S1H4/'
source("/path/to/multiomeBench/common/pbmc10k/0_create10x_datafmt.R")
library(dplyr)
train.dir <-"/path/to/data/HTAN/HT263B1-S1H1/"
train_celltype <- read.csv(paste0(train.dir,"muti_celltype.csv"))

train.dir <- "/path/to/data/HTAN/HT263B1-S1H1/"
rna_data <- readRDS(paste0(train.dir,"syn53214683/HT263B1-S1H1.rds"))
atac_data <- readRDS(paste0(train.dir,"syn53215774/HT263B1-S1H1.rds"))
# module load tabix/0.2.6 
# tabix -p bed HT263B1-S1H1-atac_fragments.tsv.gz
fpath <- paste0(train.dir,"HT263B1-S1H1-atac_fragments.tsv.gz")
rna_data<- RenameCells(rna_data, new.names = gsub(".*_", "", Cells(rna_data))) 
atac_data<- RenameCells(atac_data, new.names = gsub(".*_", "", Cells(atac_data))) 
atac_data$Original_barcode <- Cells(atac_data)
sel.bc <- intersect(Cells(rna_ata),Cells(atac_ata))
rna_data.sel <- subset(rna_data, subset = Original_barcode %in% sel.bc)
atac_data.sel <- subset(atac_data, subset = Original_barcode %in% sel.bc)

rna_count.sel <- rna_data.sel$RNA$counts
atac_count.sel <- atac_data.sel$peaksMACS2$counts
atac_count.sel <- atac_count.sel[which(rowSums(atac_count.sel)!=0),]

#################################################################
## for atac, find common peak with test data

train.peaks <- StringToGRanges(regions = rownames(atac_count.sel), sep = c(":","-"))
#######################################################################################
## load test data
# test.dir <- "/path/to/data/HTAN/HT163B1-S1H6/"
# test.fpath <- paste0(test.dir,"syn53215756/HT163B1-S1H6-fragments.tsv.gz")
# load(paste0(test.dir,"S1H6_seurat.RData"))


test.dir <- "/path/to/data/HTAN/HT243B1-S1H4/"
test.fpath <- paste0(test.dir,"HT243B1-S1H4-atac_fragments.tsv.gz")
load(paste0(test.dir,"HT243_S1H4_seurat.RData"))

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
train_peaks_feature<-  strsplit(rownames(train.frag_counts), "[-:]")

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
              output_dir <- paste0(newtest_dir,"train_forHT243B1-S1H4/")
)

train<-sub_h5(barcodes = colnames(train.frag_counts),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count.sel,train.frag_counts),
              path = paste0(newtest_dir,"train_forHT243B1-S1H4.h5"),
              test = TRUE)
###################################################################################

write.table(train_peaks_feature, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
rownames(test.frag_counts)<-train_peaks_feature$peak
sub_bc_matrix(barcodes = colnames(test.frag_counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = test.frag_counts,
              module = "Peaks",
              output_dir <- paste0(newtest_dir,"HT243B1-S1H4_commonpeak/")
)

atac_test<-sub_h5(barcodes = colnames(test.frag_counts),
                  feature_mtx = "tmp.tsv", 
                  count_mtx = test.frag_counts,
                  path = paste0(newtest_dir,"HT243B1-S1H4_commonpeaks.h5"),
                  module = "Peaks",
                  test = TRUE)

################################################################################
################################################################################
test_rna_feature<-train_feature[which(train_feature$V2 %in% rownames(all.rna$RNA$counts)),]
write.table(test_rna_feature, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
sub_bc_matrix(barcodes = colnames(all.rna$RNA$counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = all.rna$RNA$counts,
              module = 'Gene Expression',
              output_dir <- paste0(newtest_dir,"HT243B1-S1H4_rna/")
)

atac_test<-sub_h5(barcodes = colnames(all.rna$RNA$counts),
                  feature_mtx = "tmp.tsv", 
                  count_mtx =  all.rna$RNA$counts,
                  path = paste0(newtest_dir,"HT243B1-S1H4_rna.h5"),
                  module = 'Gene Expression',
                  test = TRUE)