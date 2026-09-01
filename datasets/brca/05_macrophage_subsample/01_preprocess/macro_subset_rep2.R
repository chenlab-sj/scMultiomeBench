# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# NOTE: "HT514B1-S1H3" in the paths below is a historical folder name on the
# original cluster that held COPIES of HT243B1-S1H4 files (see the HT243-named
# filenames inside it). The data are HT243B1-S1H4 -- map these paths to
# $DATA_ROOT/brca/HT243B1-S1H4/. No HT514 data is used anywhere in this benchmark.
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)

## original data
data.dir <-  "/path/to/data/HTAN/HT243B1-S1H4/"
outdir <- '/path/to/data/HTAN/HT243_S1H4_macrosub2/' ## update
## from celltype bc, subset cells
set.seed(220) ## update
celltype <- read.csv(paste0(data.dir,'muti_celltype.csv'))
celltype_macro <- subset(celltype, rna_celltype == "Macrophages")
celltype_else <- subset(celltype, rna_celltype != "Macrophages")
## subset 
celltype_subset1<- celltype_macro %>%
  group_by(rna_celltype)%>%
  sample_frac(1/2)

celltype_subset2 <- celltype_macro %>%
  group_by(rna_celltype)%>%
  sample_frac(1/2^2)

celltype_subset3 <- celltype_macro %>%
  group_by(rna_celltype)%>%
  sample_frac(1/2^3)

celltype_subset4 <- celltype_macro %>%
  group_by(rna_celltype)%>%
  sample_frac(1/2^4)

celltype_subset5 <- celltype_macro %>%
  group_by(rna_celltype)%>%
  sample_frac(1/2^5)

celltype_subset1_bc <- c(celltype_else$bc,celltype_subset1$bc)
celltype_subset2_bc <- c(celltype_else$bc,celltype_subset2$bc)
celltype_subset3_bc <- c(celltype_else$bc,celltype_subset3$bc)
celltype_subset4_bc <- c(celltype_else$bc,celltype_subset4$bc)
celltype_subset5_bc <- c(celltype_else$bc,celltype_subset5$bc)

bc_list <- list(celltype_subset1_bc, celltype_subset2_bc,celltype_subset3_bc,
             celltype_subset4_bc,celltype_subset5_bc)
save(bc_list, file = paste0(outdir,"macro_subset_bc.RData"))

print("subset R data generated")
## write to 10x
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)


library(EnsDb.Hsapiens.v86)
library(Signac)
library(Seurat)
load(paste0(data.dir,"HT243_S1H4_seurat.RData"))
atac_counts  <- all.atac$ATAC$counts
fpath <- paste0(data.dir,"HT243B1-S1H4-atac_fragments.tsv.gz")
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
load(paste0(data.dir,"HT243_S1H4_seurat.RData"))
for (i in 1:length(bc_list)){
  subset_bc <- bc_list[[i]]
  print(i)
  celltype_sub <- subset(celltype, bc %in% subset_bc)
  print(dim(celltype_sub))
  write.csv(celltype_sub, paste0(outdir, "multi_celltype_sub",i,".csv"))
  print("celltype file make")
  
  all.rna_sub <- subset(all.rna, cells = subset_bc)
  all.atac_sub <- subset(all.atac, cells = subset_bc)
  print(all.rna_sub)
  save(all.rna_sub,all.atac_sub, file = paste0(outdir,"HT243_S1H4_sub",i,".RData"))
  print("subset RData make")
  
  rna_count.sel <- all.rna_sub$RNA$counts
  atac_count.sel <- all.atac_sub$ATAC$counts
  print(dim(rna_count.sel))
  print(dim(atac_count.sel))
  rna<- AnnData(
    X = t(as.matrix(rna_count.sel)),
  )
  
  atac<- AnnData(
    X = t(as.matrix(atac_count.sel)),
  )
  write_h5ad(rna, paste0(outdir,"HT243_S1H4_rna_sub",i,"_rna.h5ad"))
  write_h5ad(atac, paste0(outdir,"HT243_S1H4_atac_sub",i,"_atac.h5ad"))
  print("h5ad create")
  
  acces.counts_sub <- acces.counts[,subset_bc]
  print(dim(acces.counts_sub))
  atac_gene<- AnnData(
    X = t(as.matrix(acces.counts_sub)),
  )
  write_h5ad(atac_gene, paste0(outdir,"S1H4_sub",i,"_atac_gene.h5ad"))
  
  print("subset h5ad atac_gene file generated")
  
  ## 
  ## normalization 
  all.atac_sub <- RunTFIDF(all.atac_sub)
  all.atac_sub <- FindTopFeatures(all.atac_sub, min.cutoff = "q0")
  all.atac_sub <- RunSVD(all.atac_sub)
  
  
  # Perform standard analysis of each modality independently RNA analysis
  all.rna_sub <- NormalizeData(all.rna_sub)
  all.rna_sub <- FindVariableFeatures(all.rna_sub)
  all.rna_sub <- ScaleData(all.rna_sub)
  all.rna_sub <- RunPCA(all.rna_sub)
  
  rna.embed <- Embeddings(all.rna_sub, reduction = "pca")
  rownames(rna.embed)<- paste0(rownames(rna.embed),"_rna")
  
  atac.embed <- Embeddings(all.atac_sub, reduction = "lsi")
  rownames(atac.embed)<- paste0(rownames(atac.embed),"_atac")
  
  print(dim(rna.embed))
  print(dim(atac.embed))
  rna.embed<- AnnData(
    X = as.matrix(rna.embed),
  )
  
  atac.embed<- AnnData(
    X =as.matrix(atac.embed),
  )
  
  write_h5ad(rna.embed, paste0(outdir,"rna_sub",i,"_embed.h5ad"))
  write_h5ad(atac.embed, paste0(outdir,"atac_sub",i,"_embed.h5ad"))
  print("embed write")
  
}








#############################################################
## for test with common peaks 
print("common peak start")
load(paste0(outdir,"macro_subset_bc.RData"))
source("/path/to/multiomeBench/common/0_create10x_datafmt.R")
common_peak<- Read10X_h5('/path/to/data/HTAN/HT514B1-S1H3/HT243B1-S1H4_commonpeak.h5')
test_rna <- Read10X_h5('/path/to/data/HTAN/HT514B1-S1H3//HT243B1-S1H4_rna.h5')

for (i in 1:length(bc_list)){
  subset_bc <- bc_list[[i]]
  print(i)
  common_peak_sub <- common_peak[,subset_bc]
  print(dim(common_peak_sub))
  sub_bc_matrix(barcodes = colnames(common_peak_sub),
                feature.mtx = "/path/to/data/HTAN/HT514B1-S1H3//HT243B1-S1H4_commonpeak/features.tsv", 
                count.mtx = common_peak_sub,
                output_dir <- paste0(outdir,"HT243B1-S1H4_commonpeak_sub",i,"/"),
                module = 'Peaks')
  
  test<-sub_h5(barcodes = colnames(common_peak_sub),
               feature_mtx ="/path/to/data/HTAN/HT514B1-S1H3//HT243B1-S1H4_commonpeak/features.tsv", 
               count_mtx = common_peak_sub,
               path = paste0(outdir,"HT243B1-S1H4_commonpeak_sub",i,".h5"),
               test = TRUE,
               module = 'Peaks')
  
  test_rna_sub <- test_rna[, subset_bc]
  print(dim(test_rna_sub))
  
  sub_bc_matrix(barcodes = colnames(test_rna_sub),
                feature.mtx = "/path/to/data/HTAN/HT514B1-S1H3/HT243B1-S1H4_rna/features.tsv", 
                count.mtx = test_rna_sub,
                output_dir =paste0(outdir,"HT243B1-S1H4_rna_sub",i,"/"),
                module = 'Gene Expression')
  
  test<-sub_h5(barcodes = colnames(test_rna_sub),
               feature_mtx = "/path/to/data/HTAN/HT514B1-S1H3/HT243B1-S1H4_rna/features.tsv", 
               count_mtx = test_rna_sub,
               path =paste0(outdir,"HT243B1-S1H4_rna_sub",i,".h5"),
               test = TRUE,
               module = 'Gene Expression')
  
  
}




