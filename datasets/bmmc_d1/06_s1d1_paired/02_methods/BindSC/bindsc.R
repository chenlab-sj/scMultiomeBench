# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(bindSC)
data.dir <- '/path/to/data/BMMC_d1/'
load(paste0(data.dir,"all_seurat.RData" ))
rna3 <- RenameCells(rna3, new.names = paste0(Cells(rna3),"_rna3"))
all.rna <- rna3

n_lat <-30
 
all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna)
all.rna <- FindNeighbors(all.rna, dims = 1:n_lat, reduction = "pca")
all.rna <- FindClusters(all.rna, resolution = 0.5)
all.rna <- RunUMAP(all.rna, dims = 1:n_lat)

atac3 <- RenameCells(atac3, new.names = paste0(Cells(atac3),"_atac3"))
all.atac <- atac3

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)
all.atac <- FindNeighbors(all.atac, dims = 1:n_lat, reduction = "lsi")
all.atac <- FindClusters(all.atac, resolution = 0.5)
all.atac <- RunUMAP(all.atac, reduction = "lsi", dims = 2:n_lat, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
gene.activities <- GeneActivity(all.atac,features = VariableFeatures(all.rna))
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- FindVariableFeatures(all.atac)
#############################################################
## integration
#############################################################
gene.use <- intersect(VariableFeatures(all.rna), 
                      VariableFeatures(all.atac))


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

df_umap <- as.data.frame(rbind(res$u, res$r))
colnames(df_umap) <- paste0("latent_",1:dim(df_umap)[2])
t2 <- Sys.time()


outdir <- "/path/to/multiomeBench/BMMC_d1/s1d1_paired/bindsc/"
if (!file.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

write.csv(df_umap,paste0(outdir,"coembed_coor.csv"))

write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

