# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).

## https://stuartlab.org/signac/articles/pbmc_multiomic.html
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
## seurat all integration
library(dplyr)
library(DropletUtils)
library(hdf5r)

#############################################
## input args
bcmtx.h5 <- "/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
fpath <-'/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
project <- "pbmc10k_control"
# # 
# train.h5 <- "/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_train.h5"
# fpath <-'/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
# project <- "pbmc10k_test_seurat4"
# test_rna.h5<-'/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5'
# test_atac.h5<-'/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5'
# # # 
#bcmtx.h5 <- '/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5'

t1 <- Sys.time()
## multiomic data
inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks



## create Seurat objet
multi <- CreateSeuratObject(counts = rna_counts)
ATAC_assay <- CreateAssayObject(counts = atac_counts)
multi[['ATAC']]<-ATAC_assay


# extract gene annotations from EnsDb
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# # convert to UCSC style
# seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# genome(annotations) <- "hg38"
# 
# # create atac object
# chrom_assay <- CreateChromatinAssay(
#   counts = atac_counts,
#   sep = c(":", "-"),
#   fragments= fpath,
#   annotation=annotations
# )
# multi[['ATAC']] <- chrom_assay
# #granges(all.atac)

## RNA analysis
DefaultAssay(multi)<-"RNA"
multi <- SCTransform(multi, verbose = FALSE) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

## ATAC analysis
DefaultAssay(multi)<-"ATAC"
multi <- RunTFIDF(multi)
multi <- FindTopFeatures(multi, min.cutoff = 'q0')
multi <- RunSVD(multi)
multi <- RunUMAP(multi, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

## integration of multiomic : WNN
multi <- FindMultiModalNeighbors(multi, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
multi <- RunUMAP(multi, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",return.model=TRUE)
multi <- RunSPCA(multi, assay = 'SCT', graph = 'wsnn') ## supervised projection

# multi@meta.data$ct <- as.character(Idents(multi))
# multi@meta.data$orig.ident <-
####################################################################
## integraion of single omic
# print("project sc-RNA to WNN graph" )
# 
# test.rna <- CreateSeuratObject(counts = test_rna_counts, assay = "RNA")
# test.rna <- SCTransform( test.rna, verbose = FALSE,assay = 'RNA')
# test.rna <- RenameCells(test.rna, add.cell.id="rna", for.merge = FALSE)
# 
# DefaultAssay(test.rna)<-'SCT'
# rna_anchors<-FindTransferAnchors(
#   reference = multi,
#   query = test.rna,
#   normalization.method = 'SCT',
#   reference.assay = 'SCT',
#   query.assay = 'SCT',
#   reference.reduction = 'spca',
#   dims = 1:50
# )
# 
# test.rna<-MapQuery(
#   anchorset = rna_anchors,
#   query = test.rna,
#   reference = multi,
#  # refdata = list(ct = "ct"), ## not right?
#   reference.reduction = 'spca',
#   reduction.model = "wnn.umap"
# )

#########################################
#print("project sc-ATAC to WNN graph" )

# create atac object
# test_chrom_assay <- CreateChromatinAssay(
#   counts = test_atac_counts,
#   sep = c(":", "-"),
#   fragments= fpath
# )
# 
# test.atac <- CreateSeuratObject(
#   counts = test_chrom_assay,
#   assay = "peaks",
# )
# Annotation(test.atac) <- annotations
# 
# DefaultAssay(multi)<-'SCT'
# gene.activities <- GeneActivity(test.atac, features = VariableFeatures(multi))
# test.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# ## normalized test.atac gene activities
# DefaultAssay(test.atac)<-"ACTIVITY"
# test.atac<- SCTransform(test.atac, verbose = FALSE,assay = 'ACTIVITY')
# test.atac<- CreateSeuratObject(counts = test_atac_counts, 
#                                sep = c(":", "-"),
#                                assay = 'ATAC')
# DefaultAssay(test.atac)<-"ATAC"
# test.atac <- RunTFIDF(test.atac)
# test.atac <- FindTopFeatures(test.atac, min.cutoff = "q0")
# test.atac <- RunSVD(test.atac)
# test.atac <- RenameCells(test.atac, add.cell.id="atac", for.merge = FALSE)


DefaultAssay(multi)<-"ATAC"
multi <- RunSLSI(multi, assay = "ATAC", graph = 'wsnn')


# atac_anchors<-FindTransferAnchors(
#   reference = multi, 
#   query = test.atac,
#   reference.assay = 'ATAC',
#   query.assay  = "ATAC",
#   reference.reduction = "slsi",
#   reduction = "lsiproject",
#   dims = 2:50
#   
# )
# 
# test.atac <- MapQuery(
#   anchorset = atac_anchors, 
#   query = test.atac,
#   reference = multi, 
#   #refdata = list(ct = "ct"),
#   reference.reduction = "slsi",
#   reduction.model = "wnn.umap"
# )

t2 <- Sys.time()


print("------ Saving integration result ------")
#df_umap = do.call(rbind,list(multi@reductions$wnn.umap@cell.embeddings,
#                             test.rna@reductions$ref.umap@cell.embeddings,
#                           test.atac@reductions$ref.umap@cell.embeddings))
df_umap_multi <- as.data.frame(multi@reductions$wnn.umap@cell.embeddings)
# df_umap_test.rna<- as.data.frame(test.rna@reductions$ref.umap@cell.embeddings)
# df_umap_test.atac<- as.data.frame(test.atac@reductions$ref.umap@cell.embeddings)

colnames(df_umap_multi) = paste0("latent_",1:ncol(df_umap_multi))
# colnames(df_umap_test.rna) = paste0("latent_",1:ncol(df_umap_test.rna))
# colnames(df_umap_test.atac) = paste0("latent_",1:ncol(df_umap_test.atac))
# df_umap = rbind(df_umap_multi,df_umap_test.rna,df_umap_test.atac)
outdir <- paste0("res_",project,"/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}

write.csv(df_umap_multi,paste0(outdir,"coembed_coor.csv"))
#saveRDS(coembed,paste0(outdir,"coembed.rds"))
write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

#############################################
## imputation of scRNA
# pred.rna <- TransferData(
#   anchorset = atac_anchors,
#   refdata = GetAssayData(multi, assay = "SCT", slot = "data"),
#   weight.reduction = "lsiproject"
# )
# 
# # add predicted values as a new assay
# test.atac[["SCT"]]<-pred.rna


## write to mtx
write10xCounts(x = multi[["SCT"]]@data, path = paste0(outdir,"data_rna.h5"))
#write10xCounts(x = pred.rna, path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_rna_pred.h5")


# ## imputation of scATAC
# pred.atac <- TransferData(
#    anchorset = atac_anchors,
#    refdata = GetAssayData(multi, assay = "ATAC", slot = "data"),
#    weight.reduction = "lsiproject"
#  )

# test.rna[["ilsi"]]<-pred.atac

write10xCounts(x = multi[["ATAC"]]@data, path = paste0(outdir,"data_atac.h5"))
# 
