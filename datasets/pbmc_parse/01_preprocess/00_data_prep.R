# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(dplyr)
library(Matrix)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(anndata)
## scRNA
## parse pbmc data
# https://www.parsebiosciences.com/datasets/performance-of-evercode-wt-mini-v3-in-human-pbmcs/#download
## d1 
out.dir <-'/path/to/multiomeBench/pbmc_parse/input/'
parse.file <- '/path/to/multiomeBench/pbmc_parse/data/'
mat <- ReadParseBio(parse.file)
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"
# Read in cell meta data
cell_meta <- read.csv(paste0(parse.file, "/cell_metadata.csv"), row.names = 1)
# Create object
pbmc_parse <- CreateSeuratObject(mat,
                           names.field = 0, meta.data = cell_meta)
pbmc_parse[["percent.mt"]] <- PercentageFeatureSet(pbmc_parse, pattern = "^MT-")
pbmc_parse[["percent.rb"]] <- PercentageFeatureSet(pbmc_parse, pattern = "^RP[SL]")
pbmc_parse <- subset(pbmc_parse, subset = nFeature_RNA > 200 & nCount_RNA < 50000 & nCount_RNA > 500 & percent.mt < 20)

## celltype annotaion, with pbmc3k-10X as reference
## input args
bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
pbmc3k.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
cell_annot.file <- '/path/to/data/pbmc3k/pbmc3k_cellannot_1113.csv'
celltype_label.file <- "/path/to/data/pbmc3k/pbmc3k_celltype_annotation.csv" 
celltype_label <- read.csv(celltype_label.file)
## ref
pbmc3k <-Read10X_h5(pbmc3k.h5)
rna_counts  <- pbmc3k$`Gene Expression`
cell_annot <- read.csv(cell_annot.file)

bc <- intersect(celltype_label$Barcode, colnames(rna_counts))
rna_counts <- rna_counts[, bc, drop = FALSE]
celltype_label<- celltype_label[match(bc, celltype_label$Barcode), , drop = FALSE]


pbmc3k_rna <- CreateSeuratObject(
  counts = rna_counts,
  project = "pbmc3k_rna",
  min.cells = 0,
  min.features = 0
)
pbmc3k_rna$celltype <-celltype_label$`Annotated.Cell.Types`
keep <- !is.na(pbmc3k_rna$celltype) & pbmc3k_rna$celltype != ""
pbmc3k_rna <- subset(pbmc3k_rna, cells = colnames(pbmc3k_rna)[keep])

## find 10x overlap transcripts from parse
train_features.file <-"/path/to/data/pbmc3k/filtered_feature_bc_matrix/features.tsv"
train_feature <- read.csv(train_features.file, sep = "\t", header = FALSE)
parse_features.file <-'/path/to/multiomeBench/pbmc_parse/data/all_genes.csv'
parse_feature <- read.csv(parse_features.file)
gene_feature <- subset(train_feature, V3 == "Gene Expression")
idx_in_10x <- match( parse_feature$gene_id,gene_feature$V1,)
keep <- !is.na(idx_in_10x)
parse_feature_10x <- parse_feature[keep,]
common.features<-Features(pbmc_parse)[keep]


#common.features <- intersect(Features(pbmc3k_rna), Features(pbmc_parse))
pbmc3k_rna_sub <- subset(pbmc3k_rna, features = common.features)
pbmc_parse_sub <- subset(pbmc_parse, features = common.features)
## celltype label transfer
pbmc3k_rna_sub<-pbmc3k_rna_sub %>% 
  NormalizeData()%>%
  FindVariableFeatures(nfeatures = 2000)%>%
  ScaleData()%>%
  RunPCA()


## query
pbmc_parse_sub <-pbmc_parse_sub %>% 
  NormalizeData()%>%
  FindVariableFeatures( nfeatures = 2000)%>%
  ScaleData()%>%
  RunPCA()




pbmc.anchors <- FindTransferAnchors(reference = pbmc3k_rna_sub, query = pbmc_parse_sub, dims = 1:30,features = common.features,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = pbmc.anchors, refdata = pbmc3k_rna_sub$celltype, dims = 1:30)
pbmc_parse_sub <- AddMetaData(pbmc_parse_sub, metadata = predictions)
pbmc_parse.d1 <- subset(pbmc_parse_sub,subset= sample == "Donor_1")





## quick visulization to check
pbmc_parse.d1 <- pbmc_parse.d1 %>%
  NormalizeData()%>%
  FindVariableFeatures()%>%
  ScaleData()%>%
  RunPCA()%>%
  FindNeighbors()%>%
  RunUMAP(dims = 1:30, reduction = "pca")

pbmc_parse.d1 <-FindNeighbors(pbmc_parse.d1, graph.name = "test")
pbmc_parse.d1 <- FindClusters(pbmc_parse.d1,graph.name = "test", resolution = 0.5)
pbmc_parse.d1 <- FindSubCluster(pbmc_parse.d1, '0', graph.name =  "test", subcluster.name = "seurat_clustsub0",  resolution =0.5)
pbmc_parse.d1 <- FindSubCluster(pbmc_parse.d1, '4', graph.name =  "test", subcluster.name = "seurat_clustsub4",  resolution =1)
pbmc_parse.d1 <- FindSubCluster(pbmc_parse.d1, '6', graph.name =  "test", subcluster.name = "seurat_clustsub6",  resolution =1)

Idents(pbmc_parse.d1)<-"seurat_cluster"

png(paste0(out.dir, "pbmc_parsed1_check.png"), width=3200, height=1250, res= 400)
p<-DimPlot(pbmc_parse.d1, reduction = "umap", group.by = c( "seurat_clusters",'predicted.id'),pt.size = 0.1)
p
dev.off()
png(paste0(out.dir, "pbmc_parsed1_checksub0.png"), width=3200, height=1250, res= 400)
p<-DimPlot(pbmc_parse.d1, reduction = "umap", group.by = c( "seurat_clustsub0",'predicted.id'),pt.size = 0.1)
p
dev.off()
png(paste0(out.dir, "pbmc_parsed1_checksub4.png"), width=3200, height=1250, res= 400)
p<-DimPlot(pbmc_parse.d1, reduction = "umap", group.by = c( "seurat_clustsub4",'predicted.id'),pt.size = 0.1)

dev.off()
## also quick check marker
T_NK <-c("CD3D","CD3E","CD3G", ## T cell
         "CD4","IL7R","CCR7","CD8A","CD8B","GZMB","GNLY",
         'NCR1','KLRD1')

B_DC <-c("MS4A1","CD79A","CD79B","CLEC9A","CLEC10A")
Mono <- c("CD14",'LYZ',"S100A9","FCGR3A","MS4A7","LST1",'CX3CR1',"CD68","CD163","MRC1",'MSR1')

p1 <- FeaturePlot(pbmc_parse.d1, features = T_NK,  reduction = "umap", ncol = 4, pt.size = 0.1)
p2 <- FeaturePlot(pbmc_parse.d1, features = B_DC,  reduction = "umap", ncol = 5, pt.size = 0.1)
p3 <- FeaturePlot(pbmc_parse.d1, features = Mono, reduction = "umap", ncol = 4, pt.size = 0.1) 
png(paste0(out.dir, "pbmc_parse_d1_marker1.png"), width=5000, height=5500, res=400)
print(p1 / p3)
dev.off()
png(paste0(out.dir, "pbmc_parse_d1_marker2.png"), width=5000, height=100, res=400)
print(p2)
dev.off()

## tmp used the predicted.id as celltype
##  together with atac_test, save file
inputdata.10x <- Read10X_h5(bcmtx.h5)
atac_counts  <- inputdata.10x$Peaks

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
all.rna <- pbmc_parse.d1
save(all.rna, all.atac, file = paste0(out.dir,"pbmc_parse_seurat.RData"))
###############################################################
## save celltype
rna_celltype <- data.frame(
  bc =Cells(pbmc_parse.d1),
  rna_celltype =pbmc_parse.d1$predicted.id
)
rna_celltype$modality = "parse_d1 scRNA"

## add celltype infor to all.atac

rownames(celltype_label)<-celltype_label$Barcode

all.atac$cell_type<-celltype_label[Cells(all.atac),]$`Annotated.Cell.Types`


atac_celltype <- data.frame(
  bc =Cells(all.atac),
  atac_celltype =all.atac$cell_type
)
atac_celltype$modality = "pbmc3k_test scATAC"
celltype <- bind_rows(rna_celltype, atac_celltype)


write.csv(celltype, paste0(out.dir,"muti_celltype.csv"),row.names = FALSE)

######################################################################

## write to 10x
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)


rna_counts <- all.rna$RNA$counts
#rna_counts_common <-rna_counts[common.features,]
## only keep overlap with 10x?
rna<- AnnData(
  X = t(as.matrix(rna_counts)),
)


atac<- AnnData(
  X = t(as.matrix(atac_counts)),
)
write_h5ad(rna, paste0(out.dir,"pbmc_parsed1_rna.h5ad"))
write_h5ad(atac, paste0(out.dir,"pbmc3k_test_atac.h5ad"))
#write_h5ad(rna_common, paste0(out.dir,"pbmc_parsed1_rna_common.h5ad"))



library(EnsDb.Hsapiens.v86)
library(Signac)
library(Seurat)


load(paste0(out.dir,"pbmc_parse_seurat.RData"))
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
write_h5ad(atac_gene, paste0(out.dir,"pbmc3k_test_atac_gene.h5ad"))


########################################################
## save as 10X multiomic dir
source("/path/to/multiomeBench/common/0_create10x_datafmt.R")
train_features.file <-"/path/to/data/pbmc3k/filtered_feature_bc_matrix/features.tsv"
train_feature <- read.csv(train_features.file, sep = "\t", header = FALSE)
# parse_features.file <-'/path/to/multiomeBench/pbmc_parse/data/all_genes.csv'
# parse_feature <- read.csv(parse_features.file)

# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# # convert to UCSC style
# seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# genome(annotations) <- "hg38"
# 
# library(GenomicRanges)
# library(dplyr)
# 
# # gene_gr is a GRangesList (per gene)
# gene_gr <- suppressWarnings(
#   reduce(split(annotations, mcols(annotations)$gene_id))
# )
# 
# # get a single "gene span" per gene: min start, max end
# gene_df <- data.frame(
#   gene_id = names(gene_gr),
#   chr = vapply(gene_gr, function(gr) as.character(seqnames(gr))[1], character(1)),
#   start = vapply(gene_gr, function(gr) min(start(gr)), integer(1)),
#   end   = vapply(gene_gr, function(gr) max(end(gr)), integer(1)),
#   stringsAsFactors = FALSE
# )
# 
# # join to parse_feature and write 10x-like features.tsv
# parse_annot <- parse_feature %>%
#   mutate(gene_id = as.character(gene_id),
#          gene_name = as.character(gene_name)) %>%
#   left_join(gene_df, by = "gene_id")
gene_feature <- subset(train_feature, V3 == "Gene Expression")
rownames(gene_feature)<-gene_feature$V1
parse_feature_10x.fmt <- gene_feature[parse_feature_10x$gene_id,]

rownames(parse_feature_10x.fmt)<-NULL
# drop genes without coordinates
# parse_feature_10x <- parse_feature_10x[!is.na(parse_feature_10x$V4), , drop = FALSE]
# keep_idx <- !duplicated(genes_target)
# parse_feature_10x <- parse_feature_10x[keep_idx, , drop = FALSE]
# genes_target <- as.character(parse_feature_10x$V2)
# genes_keep <- genes_target[genes_target %in% rownames(rna_counts)]
# rna_counts_sub <- rna_counts[genes_keep, , drop = FALSE]
# parse_feature_10x <- subset(parse_feature_10x, V2 %in% genes_keep)
write.table(
  parse_feature_10x.fmt,
  file = "tmp.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


rna_counts<-all.rna$RNA$counts
rownames(rna_counts)<-parse_feature_10x.fmt$V2
sub_bc_matrix(barcodes = colnames(rna_counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = rna_counts,
              output_dir <- paste0(out.dir,"pbmc_parse_d1/")
)

pbmc_parse_d1_test<-sub_h5(barcodes = colnames(rna_counts),
             feature_mtx = "tmp.tsv", 
             count_mtx = rna_counts,
             path = paste0(out.dir,"pbmc_parse_d1.h5"),
             test = TRUE,
             module = 'Gene Expression')


## then save atac_counts
peaks_feature <- subset(train_feature, V3 == "Peaks")

write.table(
  peaks_feature,
  file = "tmp.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sub_bc_matrix(barcodes = colnames(atac_counts),
              feature.mtx = "tmp.tsv", 
              count.mtx = atac_counts,
              output_dir <- paste0(out.dir,"pbmc3k_test_atac/")
)

pbmc3k_test_atac_test<-sub_h5(barcodes = colnames(atac_counts),
                           feature_mtx = "tmp.tsv", 
                           count_mtx = atac_counts,
                           path = paste0(out.dir,"pbmc3k_test_atac.h5"),
                           test = TRUE,
                           module = 'Peaks')

######################################################################
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

all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna)

rna.embed <- Embeddings(all.rna, reduction = "pca")
atac.embed <- Embeddings(all.atac, reduction = "lsi")



rna.embed<- AnnData(
  X = as.matrix(rna.embed),
)

atac.embed<- AnnData(
  X =as.matrix(atac.embed),
)

write_h5ad(rna.embed, paste0(out.dir,"pbmc_parse_testrna_embed.h5ad"))
write_h5ad(atac.embed, paste0(out.dir,"pbmc_parse_testatac_embed.h5ad"))
