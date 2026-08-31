# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(rliger)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
################################################################
## liger tutorial: http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html
###################################################################
## input files and pram

# #######################################################################
###############################################
bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
n_lat <-30
project <- "pbmc3k"

outdir <- paste0("res_",project,"/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}

###############################################
t1 <- Sys.time()
inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks

#bc<-colnames(rna_counts)
#atac_counts<-atac_counts[,bc]
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

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# convert to UCSC style
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations

# ## create RNA seurat object
# all.rna <- CreateSeuratObject(counts = rna_counts,assay = "RNA")
# all.rna <- NormalizeData(all.rna)
# all.rna <- FindVariableFeatures(all.rna)
#genes <- row.names(rna_counts)
gene.activities <- GeneActivity(all.atac)
# add gene activities as a new assay
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
DefaultAssay(all.atac) <- "ACTIVITY"
activity.count <- all.atac@assays[["ACTIVITY"]]@counts
## renamed bc
colnames(rna_counts)<-paste0(colnames(rna_counts),"_rna")
colnames(activity.count)<-paste0(colnames(activity.count),"_atac")

#pbmc3k.data<-list(rna=rna_counts, atac=activity.count)  
## start rliger
liger.data <- list(atac = activity.count, rna = rna_counts)
int.data <- createLiger(liger.data ) ## liger require unpair data
## preprocessing before iNMF
int.data <- rliger::normalize(int.data)
int.data <- selectGenes(int.data, datasets.use = 2)
int.data <- scaleNotCenter(int.data)
## joint matrix factorization
int.data <- optimizeALS(int.data, k = 20)
int.data <- quantile_norm(int.data)
t2 <- Sys.time()

print("------ Saving integration result ------")
## export cell embeddings
df_umap<-int.data@H.norm
colnames(df_umap) = paste0("latent_",1:ncol(df_umap))


if (!file.exists(outdir)) {
  dir.create(outdir)
}

write.csv(df_umap,paste0(outdir,"coembed_coor.csv"))
write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
