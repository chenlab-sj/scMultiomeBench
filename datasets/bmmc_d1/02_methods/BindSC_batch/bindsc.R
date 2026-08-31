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

n_lat <-30
t1 <- Sys.time()
rna1 <- RenameCells(rna1, new.names = paste0(Cells(rna1),"_rna1")) 
rna2 <- RenameCells(rna2, new.names = paste0(Cells(rna2),"_rna2"))
rna3 <- RenameCells(rna3, new.names = paste0(Cells(rna3),"_rna3")) 

all.rna <-merge(rna1, c(rna2,rna3))
all.rna <-JoinLayers(all.rna)
Idents(all.rna)<- all.rna$orig.ident
all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna)
all.rna[["RNA"]] <- split(all.rna[["RNA"]], f = all.rna$orig.ident)
all.rna  <- IntegrateLayers(object = all.rna , method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                            verbose = FALSE)
all.rna[["RNA"]] <- JoinLayers(all.rna[["RNA"]])
all.rna <- FindNeighbors(all.rna, reduction = "integrated.cca", dims = 1:30)
all.rna <- FindClusters(all.rna, resolution = 0.5)
all.rna<- RunUMAP(all.rna, dims = 1:30, reduction =  "integrated.cca")

atac1 <- RenameCells(atac1 , new.names = paste0(Cells(atac1),"_atac1")) 
atac2 <- RenameCells(atac2, new.names = paste0(Cells(atac2),"_atac2")) 
atac3 <- RenameCells(atac3, new.names = paste0(Cells(atac3),"_atac3")) 
atac1 <- FindTopFeatures(atac1, min.cutoff = 10)
atac1 <- RunTFIDF(atac1)
atac1 <- RunSVD(atac1)
atac2 <- FindTopFeatures(atac2, min.cutoff = 10)
atac2 <- RunTFIDF(atac2)
atac2 <- RunSVD(atac2)
atac3 <- FindTopFeatures(atac3, min.cutoff = 10)
atac3 <- RunTFIDF(atac3)
atac3 <- RunSVD(atac3)

all.atac<- merge(x = atac1,y = c(atac2,atac3))
#all.atac <-JoinLayers(all.atac)
Idents(all.atac)<- all.atac$orig.ident
DefaultAssay(all.atac) <- "ATAC"
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)
all.atac <- FindNeighbors(all.atac, dims = 1:n_lat, reduction = "lsi")
all.atac <- FindClusters(all.atac, resolution = 0.5)
all.atac <- RunUMAP(all.atac, reduction = "lsi", dims = 2:n_lat, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

integration.anchors <- FindIntegrationAnchors(
  object.list = list(atac1, atac2, atac3),
  anchor.features = rownames(all.atac[[1]]$ATAC$counts),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
all.atac <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = all.atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

##################################################################
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# convert to UCSC style
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations

gene.activities <- GeneActivity(all.atac, features = VariableFeatures(all.rna))

all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- FindVariableFeatures(all.atac)

gene.use <- intersect(VariableFeatures(all.rna), 
                      VariableFeatures(all.atac))



DefaultAssay(all.atac) <- "ATAC"
all.atac <- FindNeighbors(all.atac, dims = 1:n_lat, reduction = "integrated_lsi")
all.atac <- FindClusters(all.atac, resolution = 0.5)
#########################################

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
y  <- all.atac@reductions$integrated_lsi@cell.embeddings


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


outdir <- paste0(data.dir,"bindsc2/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}

write.csv(df_umap,paste0(outdir,"coembed_coor.csv"))

write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)


