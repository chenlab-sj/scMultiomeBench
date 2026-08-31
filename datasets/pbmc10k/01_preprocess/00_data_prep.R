# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## for pbmc10k


library(Signac)
library(EnsDb.Hsapiens.v86)
library(reticulate)
use_python("~/.conda/envs/seurat4/bin/python")
library(Seurat)
library(anndata)

## input args
bcmtx.h5 <- "/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
data.dir <- '/path/to/data/pbmc10k/'
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

save(all.rna, all.atac, file = paste0(data.dir,"pbmc10k_seurat.RData"))

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

write_h5ad(rna.embed, paste0(data.dir,"pbmc10k_testrna_embed.h5ad"))
write_h5ad(atac.embed, paste0(data.dir,"pbmc10k_testatac_embed.h5ad"))
