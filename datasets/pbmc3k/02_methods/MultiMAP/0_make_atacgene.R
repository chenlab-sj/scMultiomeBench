# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## use seurat to covert atac_peak.h5 to atac_gene.h5
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(DropletUtils)
library(hdf5r)
## example: Rscript 00_atac_gene.R /path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5 /path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5
#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# # test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one (input file).h5", call.=FALSE)
# } else if (length(args)==4) {
#   print("load input rna.h5 and atac.h5")
#   # default output file
#   rna.h5<-args[1]
#   atac.h5 <-args[2]
#   rna_counts <- Read10X_h5(rna.h5)
#   atac_counts <- Read10X_h5(atac.h5)
#   fpath <- args[3]
#   outfile<-args[4]
# }else if (length(args)==3){
#   print("load input multi modality h5")
#   bcmtx.h5<-args[1]
#   inputdata.10x<-Read10X_h5(bcmtx.h5)
#   rna_counts  <- inputdata.10x$`Gene Expression`
#   atac_counts  <- inputdata.10x$Peaks
#   fpath <- args[2]
#   outfile<-args[3]
# }
# #bcmtx.h5<-"/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5"

bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
data_dir <- "/path/to/tools/multimap/pbmc3k/"

inputdata.10x<-Read10X_h5(bcmtx.h5)
 atac_counts  <- inputdata.10x$Peaks

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
# all.rna <- CreateSeuratObject(counts = rna_counts,assay = "RNA")
# all.rna <- NormalizeData(all.rna)
gene.activities <- GeneActivity(all.atac)
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
acces.counts <- GetAssayData(object = all.atac, layer  = "ACTIVITY", slot = "counts")
outfile <- paste0(data_dir, "atac_gene.h5")
write10xCounts(x = acces.counts, outfile)
