# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(DropletUtils)
## seurat all integration
#############################################
set.seed(40)
## input args
bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
n_lat <-30
project <- "pbmc3k_testall"
outdir <- "/path/to/tools/Seuratv3/pbmc3k/rep3/res_pbmc3k_testall/"
if (!file.exists(outdir)) {
  dir.create(outdir)
}


###############################################




t1 <- Sys.time()
inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks

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
#granges(all.atac)
RenameCells(all.atac, add.cell.id='atac') 
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# convert to UCSC style
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations
## normalization 
all.atac <- RunTFIDF(all.atac)
all.atac <- FindTopFeatures(all.atac, min.cutoff = "q0")
all.atac <- RunSVD(all.atac,seed.use = 210)
all.atac <- RunUMAP(all.atac, reduction = "lsi", dims = 2:n_lat, reduction.name = "umap.atac", reduction.key = "atacUMAP_",seed.use = 40)
#gene.activities <- GeneActivity(all.atac)

## create RNA seurat object
all.rna <- CreateSeuratObject(counts = rna_counts,assay = "RNA")
RenameCells(all.rna, add.cell.id='rna') 
# Perform standard analysis of each modality independently RNA analysis
all.rna <- NormalizeData(all.rna)
all.rna <- FindVariableFeatures(all.rna)
all.rna <- ScaleData(all.rna)
all.rna <- RunPCA(all.rna,seed.use = 40)
all.rna <- RunUMAP(all.rna, dims = 1:n_lat,seed.use = 40)

gene.activities <- GeneActivity(all.atac, features = VariableFeatures(all.rna))
# add gene activities as a new assay
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
DefaultAssay(all.atac) <- "ACTIVITY"
all.atac <- NormalizeData(all.atac)
all.atac <- ScaleData(all.atac, features = rownames(all.atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = all.rna, query = all.atac, features = VariableFeatures(object = all.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

#note that we restrict the imputation to variable genes from scRNA-seq, but could impute the full transcriptome if we wanted to
genes.use <- VariableFeatures(all.rna)
refdata <- GetAssayData(all.rna, assay = "RNA", slot = "data")[genes.use, ]

 #refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
 #(output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = all.atac[["lsi"]],
                           dims = 2:n_lat)
all.atac[["RNA"]] <- imputation

coembed <- merge(x = all.rna, y = all.atac)
# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:n_lat)
coembed <- FindNeighbors(coembed, dims = 1:n_lat)

findres<-function(seurat.obj, nclust){
  iteration =0
  resolutions = c(0,100)
  obtained_nclust = -1
  while ((obtained_nclust != nclust) & (iteration < 100)) {
    res = sum(resolutions)/2
    tmp <- FindClusters(seurat.obj, resolution = res)
    obtained_nclust<-length(levels(tmp@active.ident))
    if (obtained_nclust < nclust){
      resolutions[1] = res
    }else{
      resolutions[2]=res
    }
    iteration = iteration + 1
    #print(iteration)
  }
  return(res)
}

res = findres(coembed,14)
coembed<-FindClusters(coembed, resolution = res)

## save cluster result
df_clust<- as.data.frame (coembed@active.ident)

write.csv(df_clust,paste0(outdir,"seurat_nclust.csv"))




t2 <- Sys.time()
print("------ Saving integration result ------")
df_umap = as.data.frame(coembed@reductions$pca@cell.embeddings[,1:n_lat])
colnames(df_umap) = paste0("latent_",1:ncol(df_umap))


write.csv(df_umap,paste0(outdir,"coembed_coor.csv"))
#saveRDS(coembed,paste0(outdir,"coembed.rds"))
write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)


#write10xCounts(x = all.rna@assays['RNA']$RNA@data, path = paste0(outdir,"rna_data.h5"))

