# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
source("/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/archr_GeneActivity.R")
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyr)

bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'

source("/path/to/tools/scJoint/data_to_h5.R")
inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_count  <- inputdata.10x$`Gene Expression`
atac_count  <- inputdata.10x$Peaks

chrom_assay <- CreateChromatinAssay(
  counts = atac_count,
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
all.rna <- CreateSeuratObject(counts = rna_count,assay = "RNA")
all.rna <- NormalizeData(all.rna)

# scJoint expects COUNT-like gene activity (its log-normalization + model are built for counts; the
# working Signac/cross-donor runs feed integer counts). ArchR gene scores are continuous (100% non-
# integer), which breaks scJoint's normalization -> modalities don't mix. Round to integer counts so
# the ATAC matrix matches what scJoint expects. (scJoint-only; the other 5 methods use raw scores.)
gene.activities <- round(archr_GeneActivity(all.atac))
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")

common.genes <- intersect(row.names(acces.counts),row.names(rna_count))
cat("ArchR common-gene count (input_size for 02_config.py):", length(common.genes), "\n")
acces.counts<-as.matrix(acces.counts)
acces.counts.sel <-acces.counts[common.genes,]
colnames(acces.counts.sel)<-paste0(colnames(acces.counts.sel),"_atac")

acces_counts.lst<-list()
acces_counts.lst[[1]]<-acces.counts.sel
write_h5_scJoint(acces_counts.lst,"/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/scJoint/scJoint_need/pbmc3k_atac_gene_scjoint.h5")

## common genes between rna and atac
rna.counts<-as.matrix(rna_count)
rna.counts.sel<-rna.counts[common.genes,]
rna_counts.lst<-list()
colnames(rna.counts.sel)<- paste0(colnames(rna.counts.sel),"_rna")
rna_counts.lst[[1]]<-rna.counts.sel
write_h5_scJoint(rna_counts.lst,"/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/scJoint/scJoint_need/pbmc3k_rna_scjoint.h5")

########################################################################

##### cell type label
########################################################################
## celltype list file
label_file<-"/path/to/data/pbmc3k/pbmc3k_cellannot_1113.csv"
label_df <- read.csv(label_file,row.names = 1)
test_label <- label_df[colnames(rna_count),]
test_celltype.lst<-list()
test_celltype.lst[[1]]<-test_label$cell_type
write_csv_scJoint(test_celltype.lst,"/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/scJoint/scJoint_need/pbmc3k_testrna_celltype_scjoint.csv")
out.dir <- "/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/scJoint/scJoint_need/"
rna_bc.file = paste0(out.dir,"rna_bc.txt")
writeLines(colnames(rna_count), rna_bc.file)
atac_bc.file = paste0(out.dir,"atac_bc.txt")
writeLines(colnames(atac_count), atac_bc.file)
unique(test_label$cell_type)
