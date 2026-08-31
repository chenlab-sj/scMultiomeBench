# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(rliger)
library(Seurat)
################################################################
## liger tutorial: http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html
###################################################################
## input files and pram
test.h5<-'/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5'
liger_pre_dir<-'/path/to/tools/liger/pbmc10k/liger_need/'
outdir <-'/path/to/tools/liger/pbmc10k/res_pbmc10k_control/'

rna_test.h5<-'/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5'
atac_test.h5<-'/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5'
liger_pre_dir<-'/path/to/tools/liger/pbmc10k/liger_need/'
outdir <-'/path/to/tools/liger/pbmc10k/res_pbmc10k_test/'
#######################################################################
t1 <- Sys.time()
## import bedmap output to calculate accsisibility counts overlapping by gene for test dataset
genes.bc <- read.table(file = paste0(liger_pre_dir,"atac_genes_bc.bed"), sep = "\t", as.is = c(4,7), header = FALSE)
promoters.bc <- read.table(file = paste0(liger_pre_dir,"atac_promoters_bc.bed"), sep = "\t", as.is = c(4,7), header = FALSE)
# bc <- genes.bc[,7]
# bc_split <- strsplit(bc,";")
# bc_split_vec <- unlist(bc_split)
# bc_unique <- unique(bc_split_vec)
# bc_counts <- table(bc_split_vec)
test.counts<-Read10X_h5(test.h5)
rna.counts<-test.counts$`Gene Expression`
atac.counts<-test.counts$Peaks  
barcodes <- colnames(atac.counts)
gene.counts <- makeFeatureMatrix(genes.bc, barcodes)
promoter.counts <- makeFeatureMatrix(promoters.bc, barcodes)
gene.counts <- gene.counts[order(rownames(gene.counts)),]
promoter.counts <- promoter.counts[order(rownames(promoter.counts)),]
acces.counts <- gene.counts + promoter.counts
colnames(acces.counts)<- paste0(colnames(acces.counts),"_atac")
colnames(rna.counts)<- paste0(colnames(rna.counts),"_rna")
#######################################################################
## start rliger
liger.data <- list(atac = acces.counts, rna = rna.counts)
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
