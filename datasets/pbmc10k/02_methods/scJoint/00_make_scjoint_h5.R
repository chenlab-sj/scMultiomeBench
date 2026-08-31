# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
source("/path/to/tools/scJoint/data_to_h5.R")
bcmtx.h5<-"/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5"
fpath <- "/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

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
gene.activities <- GeneActivity(all.atac)
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")

common.genes <- intersect(row.names(acces.counts),row.names(rna_counts))
acces.counts<-as.matrix(acces.counts)
acces.counts.sel <-acces.counts[common.genes,]

acces_counts.lst<-list()
acces_counts.lst[[1]]<-acces.counts.sel
write_h5_scJoint(acces_counts.lst,"pbmc10k_control_atac_gene_scjoint.h5")

## common genes between rna and atac
rna.counts<-as.matrix(rna_counts)
rna.counts.sel<-rna.counts[common.genes,] 
rna_counts.lst<-list()
rna_counts.lst[[1]]<-rna.counts.sel
write_h5_scJoint(rna_counts.lst,"pbmc10k_control_rna_scjoint.h5")

########################################################################

##### cell type label
########################################################################
## celltype list file
label_file<-"/path/to/multiomeBench/common/1_pbmc10k_annot0818.csv"
label_df <- read.csv(label_file,row.names = 1)
test_label <- label_df[colnames(rna.counts.sel),]
test_celltype.lst<-list()
test_celltype.lst[[1]]<-test_label$celltype
write_csv_scJoint(test_celltype.lst,"pbmc10k_control_celltype_scjoint.csv")
out.dir <- "/path/to/tools/scJoint/pbmc10k/scJoint_control_need/"
rna_bc.file = paste0(out.dir,"rna_bc.txt")
writeLines(colnames(rna_counts), rna_bc.file)
rna_bc.file = paste0(out.dir,"rna_bc.txt")
writeLines(colnames(rna_counts), rna_bc.file)
