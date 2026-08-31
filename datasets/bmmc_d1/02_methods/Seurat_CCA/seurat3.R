# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tidyr)
data.dir <- '/path/to/data/BMMC_d1/'
load(paste0(data.dir,"all_seurat.RData" ))
rna1 <- RenameCells(rna1, new.names = paste0(Cells(rna1),"_rna1")) 
rna2 <- RenameCells(rna2, new.names = paste0(Cells(rna2),"_rna2"))
rna3 <- RenameCells(rna3, new.names = paste0(Cells(rna3),"_rna3")) 

all.rna <-merge(rna1, c(rna2,rna3))
all.rna <-JoinLayers(all.rna)
n_lat <-30

all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna)
all.rna <- RunUMAP(all.rna, dims = 1:n_lat)


atac1 <- RenameCells(atac1 , new.names = paste0(Cells(atac1),"_atac1")) 
atac2 <- RenameCells(atac2, new.names = paste0(Cells(atac2),"_atac2"))
atac3 <- RenameCells(atac3, new.names = paste0(Cells(atac3),"_atac3")) 

all.atac<- merge(x = atac1,y = c(atac2,atac3))
#all.atac <-JoinLayers(all.atac)
      
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac)
all.atac <- RunUMAP(all.atac, reduction = "lsi", dims = 2:n_lat, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
gene.activities <- GeneActivity(all.atac,features = VariableFeatures(all.rna))
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- NormalizeData(all.atac)
all.atac <- ScaleData(all.atac, features = rownames(all.atac))

#############################################################
## integration
#############################################################

transfer.anchors <- FindTransferAnchors(reference = all.rna, query = all.atac, features = VariableFeatures(object = all.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
genes.use <- VariableFeatures(all.rna)
refdata <- GetAssayData(all.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = all.atac[["lsi"]],
                           dims = 2:n_lat)
all.atac$RNA <- imputation
coembed <- merge(x = all.rna, y = all.atac)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
t2 <- Sys.time()
df_umap = as.data.frame(coembed@reductions$pca@cell.embeddings[,1:n_lat])
colnames(df_umap) = paste0("latent_",1:ncol(df_umap))


outdir <- paste0(data.dir,"seurat3/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}


write.csv(df_umap,paste0(outdir,"coembed_coor.csv"))

write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)


