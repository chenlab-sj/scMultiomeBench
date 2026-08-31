# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
data.dir <- "/path/to/data/HTAN/HT243B1-S1H4/"
rna_data <- readRDS(paste0(data.dir,"syn53214720/HT243B1-S1H4.rds"))
atac_data <- readRDS(paste0(data.dir,"syn53215789/HT243B1-S1H4.rds"))
# module load tabix/0.2.6 
# tabix -p bed HT263B1-S1H1-atac_fragments.tsv.gz
fpath <- paste0(data.dir,"HT243B1-S1H4-atac_fragments.tsv.gz")
rna_data<- RenameCells(rna_data, new.names = gsub(".*_", "", Cells(rna_data))) 
atac_data<- RenameCells(atac_data, new.names = gsub(".*_", "", Cells(atac_data))) 
atac_data$Original_barcode <- Cells(atac_data)
sel.bc <- intersect(Cells(rna_data),Cells(atac_data))
rna_data.sel <- subset(rna_data, subset = Original_barcode %in% sel.bc)
atac_data.sel <- subset(atac_data, subset = Original_barcode %in% sel.bc)

rna_count.sel <- rna_data.sel$RNA$counts
atac_count.sel <- atac_data.sel$peaksMACS2$counts
atac_count.sel <- atac_count.sel[which(rowSums(atac_count.sel)!=0),]
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

save(all.rna, all.atac, file = paste0(data.dir,"HT243_S1H4_seurat.RData"))

#############################################################
## save celltype
rna_celltype1 <- data.frame(
  bc =rna_data1$Original_barcode,
  rna_celltype =rna_data1$cell_type
)
rna_celltype1$modality = "HT243B1-S1H4 scRNA"
atac_celltype1 <- data.frame(
  bc =atac_data1$Original_barcode,
  atac_celltype =atac_data1$cell_type
)
atac_celltype1$modality = "HT243B1-S1H4 scATAC"
celltype1 <- merge(rna_celltype1, atac_celltype1, by = "bc")


write.csv(celltype1, paste0(data.dir1,"muti_celltype.csv"),row.names = FALSE)

#################################################################

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
write_h5ad(rna, paste0(data.dir,"HT243_S1H4_rna.h5ad"))
write_h5ad(atac, paste0(data.dir,"HT243_S1H4_atac.h5ad"))

library(EnsDb.Hsapiens.v86)
library(Signac)
library(Seurat)


data.dir <- "/path/to/data/HTAN/HT243B1-S1H4/"
load(paste0(data.dir,"HT243_S1H4_seurat.RData"))
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
write_h5ad(atac_gene, paste0(data.dir,"S1H4_atac_gene.h5ad"))


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