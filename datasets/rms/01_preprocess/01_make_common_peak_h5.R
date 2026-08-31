# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# match feature tsv
library(tidyr)
library(dplyr)
library(Matrix)
library(rhdf5)

data.dir <- "/path/to/data/RMS/"
ref_feature_tsv<- paste0(data.dir,"Mast39_TB12_1442/Mast39_TB12_1442/outs/filtered_feature_bc_matrix/features.tsv")
ref_feature <- read.table(ref_feature_tsv, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

##find common peaks between reference and test and make new data seperately
library(Seurat)
library(Signac)
test.h5<-  paste0(data.dir,"Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5") 
train.h5 <-  paste0(data.dir,"Mast39_TB12_1442/Mast39_TB12_1442/outs/filtered_feature_bc_matrix.h5") 
train2.h5 <- paste0(data.dir,"MAST213F_TB15_5705/MAST213F_TB15_5705/outs/filtered_feature_bc_matrix.h5") 

train_fpath <- paste0(data.dir,"Mast39_TB12_1442/Mast39_TB12_1442/outs/atac_fragments.tsv.gz")
train2_fpath <- paste0(data.dir,"MAST213F_TB15_5705/MAST213F_TB15_5705/outs/atac_fragments.tsv.gz")
test_fpath <-paste0(data.dir,"Mast607A_TB19_22652/Mast607A_TB19_22652/outs/atac_fragments.tsv.gz")
inputdata.10x <- Read10X_h5(train.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks

inputdata2.10x <- Read10X_h5(train2.h5)
rna2_counts  <- inputdata2.10x$`Gene Expression`
atac2_counts  <- inputdata2.10x$Peaks

testdata.10x<-Read10X_h5(test.h5)
test_rna_count<-testdata.10x$`Gene Expression`
test_atac_count<-testdata.10x$Peaks
  
train.peaks <- StringToGRanges(regions = rownames(atac_counts), sep = c(":","-"))
train2.peaks <- StringToGRanges(regions = rownames(atac2_counts), sep = c(":","-"))
test.peaks <- StringToGRanges(regions = rownames(test_atac_count), sep = c(":","-"))

# train.peaks2 <- StringToGRanges(regions = paste0("chr",rownames(atac_counts)), sep = c(":","-"))
# train2.peaks2 <- StringToGRanges(regions = paste0("chr",rownames(atac2_counts)), sep = c(":","-"))
# test.peaks2 <- StringToGRanges(regions = paste0("chr",rownames(test_atac_count)), sep = c(":","-"))

combined.peaks <- reduce(x = c(train.peaks,train2.peaks,test.peaks))
#combined.peaks2 <- reduce(x = c(train.peaks2,train2.peaks2,test.peaks2))
# Newcounts based on unions

# create fragment objects
train.frags <- CreateFragmentObject(
  path = train_fpath,
  cells = colnames(atac_counts)
)

train2.frags <- CreateFragmentObject(
  path = train2_fpath,
  cells = colnames(atac2_counts)
)

test.frags <- CreateFragmentObject(
  path = test_fpath,
  cells = colnames(test_atac_count)
)

train.frag_counts <- FeatureMatrix(
  fragments = train.frags,
  features = combined.peaks,
  cells = colnames(atac_counts)
)

train2.frag_counts <- FeatureMatrix(
  fragments = train2.frags,
  features = combined.peaks,
  cells = colnames(atac2_counts)
)

test.frag_counts <- FeatureMatrix(
  fragments = test.frags,
  features = combined.peaks,
  cells = colnames(test_atac_count)
)

colnames(train.frag_counts)<- paste0(colnames(train.frag_counts),"_Mast39")
colnames(train2.frag_counts)<- paste0(colnames(train2.frag_counts),"_Mast213F")
train.frag_counts <- cbind(train.frag_counts,train2.frag_counts)
## make peak feature for reference
# train_peaks_feature2<-  strsplit(rownames(train.frag_counts), "-")
# train_peaks_feature2<- as.data.frame(do.call(rbind, train_peaks_feature2))
# train_peaks_feature2$peak <- paste0("chr",train_peaks_feature2$V1,":",train_peaks_feature2$V2,"-",train_peaks_feature2$V3)
# train_peaks_feature2$peak2 <- train_peaks_feature2$peak
# train_peaks_feature2$chr <- paste0("chr", train_peaks_feature2$V1)
# train_peaks_feature2$type <-"Peaks"
# train_peaks_feature2.sel <- train_peaks_feature2 %>%
#   dplyr:: select(peak, peak2,type,chr, V2, V3)
# 
# train_rna_feature <- ref_feature[ref_feature$V3 == "Gene Expression",]
library(data.table)
# train_feature <- rbindlist(list (train_rna_feature, train_peaks_feature2.sel), use.names = FALSE)
# write.table(train_feature, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
#rownames(train.frag_counts)<-train_peaks_feature2.sel$peak

colnames(rna_counts)<- paste0(colnames(rna_counts),"_Mast39")
colnames(rna2_counts)<- paste0(colnames(rna2_counts),"_Mast213F")
train.rna_counts <- cbind(rna_counts, rna2_counts)
source("/path/to/multiomeBench/common/0_create10x_datafmt.R")


########################################################################################
##### remake peak feature for test
test_peaks_feature2<-  strsplit(rownames(test.frag_counts), "-")
test_peaks_feature2<- as.data.frame(do.call(rbind, test_peaks_feature2))
test_peaks_feature2$peak <- paste0(test_peaks_feature2$V1,":",test_peaks_feature2$V2,"-",test_peaks_feature2$V3)
test_peaks_feature2$peak2 <- test_peaks_feature2$peak

test_peaks_feature2$type <-"Peaks"
test_peaks_feature2.sel <- test_peaks_feature2 %>%
  dplyr:: select(peak, peak2,type,V1, V2, V3)
train_rna_feature <- ref_feature[ref_feature$V3 == "Gene Expression",]
train_feature <- rbindlist(list (train_rna_feature, test_peaks_feature2.sel), use.names = FALSE)
write.table(train_feature, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
rownames(test.frag_counts)<-test_peaks_feature2.sel$peak
rownames(train.frag_counts)<-test_peaks_feature2.sel$peak
sub_bc_matrix(barcodes = colnames(train.frag_counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = rbind(train.rna_counts,train.frag_counts),
              
              output_dir <- "Mast39_Mast213F_filtered_feature_bc_matrix_Mast607A/"
)

train<-sub_h5(barcodes = colnames(train.frag_counts),
              feature_mtx = "tmp.tsv", 
              count_mtx =  rbind(train.rna_counts,train.frag_counts),
              path = "Mast39_Mast213F_filtered_feature_bc_matrix_Mast607A.h5",
              test = TRUE)


sub_bc_matrix(barcodes = colnames(test.frag_counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = rbind(test_rna_count,test.frag_counts),
              output_dir <-  "filtered_feature_bc_matrix_Mast607A_TB19_22652_commonpeaks/")

test<-sub_h5(barcodes = colnames(test_rna_count),
             feature_mtx = "tmp.tsv", 
             count_mtx =  rbind(test_rna_count,test.frag_counts),
             path = "filtered_feature_bc_matrix_Mast607A_TB19_22652_commonpeaks.h5",
             test = TRUE)


##########################################################
library(Signac)
library(EnsDb.Hsapiens.v86)
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)
library(EnsDb.Hsapiens.v75)

## input args
bcmtx.h5 <- "/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5"
fpath <-'/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/atac_fragments.tsv.gz'
data.dir <- '/path/to/data/RMS/Mast607A_TB19_22652/'
n_lat <-30


inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks

## prepare for R
all.rna <- CreateSeuratObject(counts = rna_counts,assay = "RNA", project = "scRNA")
#colnames(rna_counts)<-paste0(colnames(rna_counts),"_rna")
#colnames(atac_counts)<-paste0(colnames(atac_counts),"_atac")
# create atac object
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments= fpath
)

all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

save(all.rna, all.atac, file = paste0(data.dir,"Mast607A_seurat.RData"))

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# convert to UCSC style
#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg37"
Annotation(all.atac) <- annotations
## normalization 
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)

########################################################################

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

