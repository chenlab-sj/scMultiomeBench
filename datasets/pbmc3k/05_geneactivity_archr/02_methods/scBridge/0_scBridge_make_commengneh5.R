# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
source("/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/archr_GeneActivity.R")
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
#library(DropletUtils)
library(Matrix)
library(reticulate); use_python("~/.conda/envs/seurat4/bin/python")   # bind reticulate to seurat4's python (has anndata) -> write_h5ad works
library(anndata)
#source("/path/to/tools/scJoint/data_to_h5.R")
bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
out.dir <- "/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/scBridge/"
label_file<-"/path/to/data/pbmc3k/pbmc3k_cellannot_1113.csv"

bcmtx_counts <- Read10X_h5(bcmtx.h5)
rna_counts <- bcmtx_counts$`Gene Expression`
atac_counts <- bcmtx_counts$Peaks


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

# create RNA seurat object
all.rna <- CreateSeuratObject(counts = rna_counts,assay = "RNA")
all.rna <- NormalizeData(all.rna)
gene.activities <- archr_GeneActivity(all.atac)
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")

common.genes <- intersect(row.names(acces.counts),row.names(rna_counts))
acces.counts<-as.matrix(acces.counts)
acces.counts.sel <-acces.counts[common.genes,]
acces.counts.sel <- as(acces.counts.sel, "sparseMatrix")
#write10xCounts(x = acces.counts.sel , path = paste0(out.dir,"atac_comgene.h5"))
colnames(acces.counts.sel)<-paste0(colnames(acces.counts.sel),"_atac")
atac.adata <- AnnData(X= t(acces.counts.sel))
atac.adata$var$gene_ids <- common.genes
write_h5ad(atac.adata, paste0(out.dir,"scBridge_atac.h5ad"))
# acces_counts.lst<-list()
# acces_counts.lst[[1]]<-acces.counts.sel
#write_h5_scJogene_idsint(acces_counts.lst, paste0(out.dir, "atac_comgene.h5"))
#write10xCounts(x = rna.counts.sel , path = paste0(out.dir,"rna_comgene.h5"))
## common genes between rna and atac
rna.counts<-as.matrix(rna_counts)
rna.counts.sel<-rna.counts[common.genes,]
rna.counts.sel <- as(rna.counts.sel, "sparseMatrix")
colnames(rna.counts.sel)<-paste0(colnames(rna.counts.sel),"_rna")
# rna_counts.lst<-list()
# rna_counts.lst[[1]]<-rna.counts.sel
#write_h5_scJoint(rna.counts.sel,paste0(out.dir,"rna_comgene.h5"))
#write10xCounts(x = rna.counts.sel , path = paste0(out.dir,"scB.h5"))
label_df <- read.csv(label_file,row.names = 1)
test_label <- label_df[colnames(rna.counts.sel),]
test_label$cell_type <- ifelse(is.na(test_label$cell_type), "unknown",test_label$cell_type)

rna.adata <- AnnData(X= t(rna.counts.sel))
rna.adata$obs$CellType <- test_label$cell_type
rna.adata$var$gene_ids <- common.genes
write_h5ad(rna.adata, paste0(out.dir,"scBridge_rna.h5ad"))
