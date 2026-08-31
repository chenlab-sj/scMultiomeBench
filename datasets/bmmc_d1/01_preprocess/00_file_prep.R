# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
data.dir1 <- "/path/to/data/BMMC/NCBI_sra/s2d1/"
fpath1 <-"/path/to/data/BMMC/NCBI_sra/s2d1/outs/atac_fragments.tsv.gz"
bcmx1 <- Read10X_h5(paste0(data.dir1,"outs/filtered_feature_bc_matrix.h5"))
rna_count1  <- bcmx1$`Gene Expression`
atac_count1  <- bcmx1$Peaks

########################################################

data.dir2 <- "/path/to/data/BMMC/NCBI_sra/s4d1/"
fpath2 <-"/path/to/data/BMMC/NCBI_sra/s4d1/outs/atac_fragments.tsv.gz"
bcmx2 <- Read10X_h5(paste0(data.dir2,"outs/filtered_feature_bc_matrix.h5"))
rna_count2  <- bcmx2$`Gene Expression`
atac_count2  <- bcmx2$Peaks


####################################################################
data.dir3 <- "/path/to/data/BMMC/NCBI_sra/s1d1/"
fpath3 <-"/path/to/data/BMMC/NCBI_sra/s1d1/outs/atac_fragments.tsv.gz"
bcmx3 <- Read10X_h5(paste0(data.dir3,"outs/filtered_feature_bc_matrix.h5"))
rna_count3  <- bcmx3$`Gene Expression`
atac_count3  <- bcmx3$Peaks


####################################################################
## for atac, find common peak between test1 and test2

peaks1 <- StringToGRanges(regions = rownames(atac_count1), sep = c(":","-"))
peaks2 <- StringToGRanges(regions = rownames(atac_count2), sep = c(":","-"))
peaks3 <- StringToGRanges(regions = rownames(atac_count3), sep = c(":","-"))

combined.peaks <- reduce(x = c(peaks1,peaks2, peaks3))

# create fragment objects
frags1 <- CreateFragmentObject(
  path = fpath1,
  cells = colnames(atac_count1)
)
frags2 <- CreateFragmentObject(
  path = fpath2,
  cells = colnames(atac_count2)
)
frags3 <- CreateFragmentObject(
  path = fpath3,
  cells = colnames(atac_count3)
)
frag_counts1 <- FeatureMatrix(
  fragments = frags1,
  features = combined.peaks,
  cells = colnames(atac_count1)
)
frag_counts2 <- FeatureMatrix(
  fragments = frags2,
  features = combined.peaks,
  cells = colnames(atac_count2)
)
frag_counts3 <- FeatureMatrix(
  fragments = frags3,
  features = combined.peaks,
  cells = colnames(atac_count3)
)




# create atac object

chrom_assay1 <- CreateChromatinAssay(
  counts = frag_counts1,
  sep = c(":", "-"),
  fragments= fpath1
)

chrom_assay2 <- CreateChromatinAssay(
  counts = frag_counts2,
  sep = c(":", "-"),
  fragments= fpath2
)

chrom_assay3 <- CreateChromatinAssay(
  counts = frag_counts3,
  sep = c(":", "-"),
  fragments= fpath3 
)




atac1 <- CreateSeuratObject(
  counts = chrom_assay1,
  assay = "ATAC",
  project = "s2d1 scATAC"
)

atac2 <- CreateSeuratObject(
  counts = chrom_assay2,
  assay = "ATAC",
  project = "s4d1 scATAC"
)

atac3 <- CreateSeuratObject(
  counts = chrom_assay3,
  assay = "ATAC",
  project = "s1d1 scATAC"
)


rna1 <- CreateSeuratObject(counts = rna_count1,assay = "RNA", project = "s2d1 scRNA")
rna2 <- CreateSeuratObject(counts = rna_count2,assay = "RNA", project = "s4d1 scRNA")
rna3 <- CreateSeuratObject(counts = rna_count3,assay = "RNA", project = "s1d1 scRNA")

outdir <- '/path/to/data/BMMC_d1/'
save(rna1, atac1, rna2, atac2, rna3,atac3, file = paste0(outdir,"all_seurat.RData" ))
#######################################################################################

chrom_assay1 <- CreateChromatinAssay(
  counts = frag_counts1,
  sep = c(":", "-"),
  fragments= fpath1
)

chrom_assay2 <- CreateChromatinAssay(
  counts = frag_counts2,
  sep = c(":", "-"),
  fragments= fpath2
)

chrom_assay3 <- CreateChromatinAssay(
  counts = frag_counts3,
  sep = c(":", "-"),
  fragments= fpath3 
)

atac1 <- CreateSeuratObject(
  counts = chrom_assay1,
  assay = "ATAC",
  project = "s2d1 scATAC"
)

atac2 <- CreateSeuratObject(
  counts = chrom_assay2,
  assay = "ATAC",
  project = "s4d1 scATAC"
)

atac3 <- CreateSeuratObject(
  counts = chrom_assay3,
  assay = "ATAC",
  project = "s1d1 scATAC"
)





atac1 <- RenameCells(atac1 , new.names = paste0(Cells(atac1),"_atac1")) 
atac2 <- RenameCells(atac2, new.names = paste0(Cells(atac2),"_atac2"))
atac3 <- RenameCells(atac3, new.names = paste0(Cells(atac3),"_atac3")) 
all.atac<- merge(x = atac1,y = c(atac2,atac3))




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
write_h5ad(atac_gene, paste0(outdir,"test_atac_gene.h5ad"))








colnames(rna_count1)<- paste0(colnames(rna_count1),"_rna1")
colnames(rna_count2)<- paste0(colnames(rna_count2),"_rna2")
colnames(rna_count3)<- paste0(colnames(rna_count3),"_rna3")
test_rna <- cbind(rna_count1,rna_count2,rna_count3)
colnames(frag_counts1)<- paste0(colnames(frag_counts1),"_atac1")
colnames(frag_counts2)<- paste0(colnames(frag_counts2),"_atac2")
colnames(frag_counts3)<- paste0(colnames(frag_counts3),"_atac3")
test_atac <- cbind(frag_counts1,frag_counts2,frag_counts3)
## write to 10x
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)

test_rna<- AnnData(
  X = t(as.matrix(test_rna)),
)

test_atac<- AnnData(
  X = t(as.matrix(test_atac)),
)
write_h5ad(test_rna, paste0(outdir,"test_rna.h5ad"))
write_h5ad(test_atac, paste0(outdir,"test_atac.h5ad"))


chrom_assay1 <- CreateChromatinAssay(
  counts = frag_counts1,
  fragments= frags1 
)

chrom_assay2 <- CreateChromatinAssay(
  counts = frag_counts2,
  fragments= frags2 
)

chrom_assay3 <- CreateChromatinAssay(
  counts = frag_counts3,
  fragments= frags3 
)










#######################################################################################

## create multiomics h5 file
source("/path/to/multiomeBench/common/0_create10x_datafmt.R")
library(dplyr)
train_features.file <- "/path/to/data/BMMC/NCBI_sra/s1d2/outs/filtered_feature_bc_matrix/features.tsv"
train_feature <- read.csv(train_features.file, sep = "\t", header = FALSE)
rna1_feature<-train_feature[which(train_feature$V2 %in% rownames(rna_count1)),]
peaks_feature<-  strsplit(rownames(frag_counts1), "-")
peaks_feature<- as.data.frame(do.call(rbind, peaks_feature))
peaks_feature$peak <- paste0(peaks_feature$V1,":",peaks_feature$V2,"-",peaks_feature$V3)
peaks_feature$peak2 <- peaks_feature$peak
peaks_feature$type <-"Peaks"
peaks_feature <- peaks_feature %>%
  dplyr:: select(peak, peak2,type,V1, V2, V3)

library(data.table)
feature1 <- rbindlist(list (rna1_feature, peaks_feature), use.names = FALSE)
write.table(feature1, file = "tmp.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
rownames(frag_counts1)<-peaks_feature$peak

test1<-sub_h5(barcodes = colnames(frag_counts1),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count1,frag_counts1),
              path = paste0(outdir,"test1_s2d1.h5"),
              test = TRUE)
rownames(frag_counts2)<-peaks_feature$peak
test2<-sub_h5(barcodes = colnames(frag_counts2),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count2,frag_counts2),
              path = paste0(outdir,"test2_s4d1.h5"),
              test = TRUE)
rownames(frag_counts3)<-peaks_feature$peak
test3<-sub_h5(barcodes = colnames(frag_counts3),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count3,frag_counts3),
              path = paste0(outdir,"test3_s1d1.h5"),
              test = TRUE)
#######################################################################################












#################################################################################
## load train data
data.dir <- "/path/to/data/BMMC/NCBI_sra/s1d2/"
fpath <-"/path/to/data/BMMC/NCBI_sra/s1d2/outs/atac_fragments.tsv.gz"
bcmx <- Read10X_h5(paste0(data.dir,"outs/filtered_feature_bc_matrix.h5"))
rna_count  <- bcmx$`Gene Expression`
atac_count  <- bcmx$Peaks
train.peaks <- StringToGRanges(regions = rownames(atac_count), sep = c(":","-"))

combined.peaks <- reduce(x = c(train.peaks,peaks1,peaks2,peaks3))


train.frags <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_count)
)

frags1 <- CreateFragmentObject(
  path = fpath1,
  cells = colnames(atac_count1)
)
frags2 <- CreateFragmentObject(
  path = fpath2,
  cells = colnames(atac_count2)
)
frags3 <- CreateFragmentObject(
  path = fpath3,
  cells = colnames(atac_count3)
)

train.frag_counts <- FeatureMatrix(
  fragments = train.frags,
  features = combined.peaks,
  cells = colnames(atac_count)
)

frag_counts1 <- FeatureMatrix(
  fragments = frags1,
  features = combined.peaks,
  cells = colnames(atac_count1)
)
frag_counts2 <- FeatureMatrix(
  fragments = frags2,
  features = combined.peaks,
  cells = colnames(atac_count2)
)

frag_counts3 <- FeatureMatrix(
  fragments = frags3,
  features = combined.peaks,
  cells = colnames(atac_count3)
)
#######################################################################
source("/path/to/multiomeBench/common/0_create10x_datafmt.R")
library(dplyr)

train_features.file <- paste0(data.dir,"outs/filtered_feature_bc_matrix/features.tsv")
train_feature <- read.csv(train_features.file, sep = "\t", header = FALSE)

train_rna_feature<-train_feature[which(train_feature$V2 %in% rownames(rna_count)),]
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
              count.mtx = rbind(rna_count,train.frag_counts),
              output_dir <- paste0(outdir,"s1d2_train/")
)

train<-sub_h5(barcodes = colnames(train.frag_counts),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count,train.frag_counts),
              path = paste0(outdir,"s1d2_train.h5"),
              test = TRUE)

###################################################################################
colnames(rna_count1)<- paste0(colnames(rna_count1),"_rna1")
colnames(rna_count2)<- paste0(colnames(rna_count2),"_rna2")
colnames(rna_count3)<- paste0(colnames(rna_count3),"_rna3")
test_rna <- cbind(rna_count1,rna_count2,rna_count3)
colnames(frag_counts1)<- paste0(colnames(frag_counts1),"_atac1")
colnames(frag_counts2)<- paste0(colnames(frag_counts2),"_atac2")
colnames(frag_counts3)<- paste0(colnames(frag_counts3),"_atac3")
test_atac <- cbind(frag_counts1,frag_counts2,frag_counts3)
## write to 10x
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)

test_rna<- AnnData(
  X = t(as.matrix(test_rna)),
)

test_atac<- AnnData(
  X = t(as.matrix(test_atac)),
)
write_h5ad(test_rna, paste0(outdir,"test_rna.h5ad"))
write_h5ad(test_atac, paste0(outdir,"test_atac.h5ad"))

train_rna<- AnnData(
  X = t(as.matrix(rna_count)),
)

write_h5ad(train_rna, paste0(outdir,"train_rna.h5ad"))

train_atac<- AnnData(
  X = t(as.matrix(train.frag_counts)),
)

write_h5ad(train_rna, paste0(outdir,"train_atac.h5ad"))

###################################################################################
rownames(frag_counts1)<-train_peaks_feature$peak
sub_bc_matrix(barcodes = colnames(frag_counts1),
              feature.mtx = "tmp.tsv", 
              count.mtx = rbind(rna_count1,frag_counts1),
              output_dir <- paste0(outdir,"train_test1_s2d1/")
)
test1<-sub_h5(barcodes = colnames(frag_counts1),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count1,frag_counts1),
              path = paste0(outdir,"train_test1_s2d1.h5"),
              test = TRUE)


rownames(frag_counts2)<-train_peaks_feature$peak
sub_bc_matrix(barcodes = colnames(frag_counts2),
              feature.mtx = "tmp.tsv", 
              count.mtx = rbind(rna_count2,frag_counts2),
              output_dir <- paste0(outdir,"train_test2_s4d1/")
)
test2<-sub_h5(barcodes = colnames(frag_counts2),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count2,frag_counts2),
              path = paste0(outdir,"train_test2_s4d1.h5"),
              test = TRUE)

rownames(frag_counts3)<-train_peaks_feature$peak
sub_bc_matrix(barcodes = colnames(frag_counts3),
              feature.mtx = "tmp.tsv", 
              count.mtx = rbind(rna_count3,frag_counts3),
              output_dir <- paste0(outdir,"train_test3_s1d1/")
)
test3<-sub_h5(barcodes = colnames(frag_counts3),
              feature_mtx = "tmp.tsv", 
              count_mtx = rbind(rna_count3,frag_counts3),
              path = paste0(outdir,"train_test3_s1d1.h5"),
              test = TRUE)

###################################################################################
#label_file prep

meta_file<-"/path/to/data/BMMC/BMMC_meta.csv"
meta<- read.csv(meta_file,row.names = 1)
meta.sel <- meta%>%
  filter(batch %in% c(s1d2,s1d1,s2d1,s4d1))%>%
  select(batch, cell_type)