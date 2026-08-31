# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(bindSC)
#####################################
bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
set.seed(40)
n_lat <-30
project <- "res_pbmc3k"
####################################
## prepare input from seurat
#rna_test.h5 <-'/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5'
#atac_test.h5 <-'/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5'
#fpath <-'/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
#n_lat <-30
#project <- "res_pbmc10k_test"

t1 <- Sys.time()
inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks

## create RNA seurat object
all.rna <- CreateSeuratObject(counts = rna_counts,assay = "RNA")
# Perform standard analysis of each modality independently RNA analysis
DefaultAssay(all.rna) <- "RNA"
all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna)
all.rna <- FindNeighbors(all.rna, dims = 1:n_lat, reduction = "pca")
all.rna <- FindClusters(all.rna, resolution = 0.5)
all.rna <- RunUMAP(all.rna, dims = 1:n_lat)

# create atac object
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments= fpath
)

all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
)
#granges(all.atac)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# convert to UCSC style
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations

gene.activities <- GeneActivity(all.atac, features = VariableFeatures(all.rna))
# add gene activities as a new assay
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- FindVariableFeatures(all.atac)

DefaultAssay(all.atac) <- "ATAC"
## normalization 
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)
all.atac <- FindNeighbors(all.atac, dims = 1:n_lat, reduction = "lsi")
all.atac <- FindClusters(all.atac, resolution = 0.5)
all.atac <- RunUMAP(all.atac, reduction = "lsi", dims = 2:n_lat, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
#gene.activities <- GeneActivity(all.atac)


DefaultAssay(all.atac) <- "ACTIVITY"
gene.use <- intersect(VariableFeatures(all.rna), 
                      VariableFeatures(all.atac))
#####################################################
## start bindSC
## X: gene expression matrix
## Y accessibility matrix
## gene activity matrix


X <- all.rna$RNA$data[gene.use,]
Y <- all.atac[["ATAC"]]$data[]
Z0 <- all.atac$ACTIVITY$data[gene.use,]
type <- c(rep("RNA", ncol(X)), rep("ATAC", ncol(X)))
a <- rowSums(as.matrix(Y))
out <- dimReduce(dt1 =  X, dt2 = Z0,  K = 30)

x <- out$dt1
z0 <- out$dt2
y  <- all.atac@reductions$lsi@cell.embeddings


res <- BiCCA( X = t(x) ,
              Y = t(y), 
              Z0 =t(z0), 
              X.clst = all.rna$seurat_clusters,
             Y.clst = all.atac$seurat_clusters,
              alpha = 0.5, 
              lambda = 0.5,
              K = 15,
              temp.path  = "out",
              num.iteration = 50,
              tolerance = 0.01,
              save = TRUE,
              parameter.optimize = FALSE,
              block.size = 0)

## export co-embeddings
df_umap <- as.data.frame(rbind(res$u, res$r))
colnames(df_umap) <- paste0("latent_",1:dim(df_umap)[2])
t2 <- Sys.time()

print("------ Saving integration result ------")
outdir <- paste0("res_",project,"/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}

write.csv(df_umap,paste0(outdir,"coembed_coor.csv"))

write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

############################################
# ## imputation of the full scRNA features
# Z_impu <- impuZ(X=all.rna$RNA$data, bicca = res)
# # whole range normalization, ran in plot_geneScoreChange 
# #Z_impu_norm<- (Z_impu-min(Z_impu))/(max(Z_impu)-min(Z_impu))
# write10xCounts(x = Z_impu, path = paste0(outdir,"imputed_rna.h5"))
