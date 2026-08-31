# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
library(Seurat)
library(Matrix)
library(rhdf5)

source("/path/to/multiomeBench/common/0_create10x_datafmt.R")
set.seed(420)
## 

## spilt a dataset to two parts to test integration (random)
bcmtx.h5 <- "/path/to/data/pbmc10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"


inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks


sample <- sample(c(TRUE, FALSE), ncol(rna_counts), replace=TRUE, prob=c(0.4,0.6))
train.rna <- rna_counts[,sample]
test.rna  <- rna_counts[,!sample]
train.bc<-colnames(train.rna)
test.bc<-colnames(test.rna)
train.atac<- atac_counts[,train.bc]
test.atac<- atac_counts[,test.bc]
## further spilt test data for rna and atac only
test_spilt<- sample(c(TRUE, FALSE), ncol(test.rna), replace=TRUE, prob=c(0.5,0.5))
test1.rna <- test.rna[,test_spilt]
test2.rna <- test.rna[,!test_spilt]

test1.bc<- colnames(test1.rna)
test2.bc<-colnames(test2.rna)
test1.atac <-test.atac[,test1.bc]
test2.atac <-test.atac[,test2.bc]

## save tran and test bc
x<-NULL
x$train.bc<-train.bc
x$test1.bc<-test1.bc
x$test2.bc<-test2.bc
saveRDS(x,"pbmc10k_data_design_bc0802.rds")



## rename test rna and test atac barcode to ignore the multiomic infor
#colnames(test.rna)<-paste0(colnames(test.rna),"_rna")
#colnames(test.atac)<-paste0(colnames(test.rna),"_atac")
## save as 10xh5 file

train.10x<-sub_h5(barcodes = train.bc,
                  feature_mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv",
                  count_mtx = rbind(train.rna,train.atac),
                  path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_train.h5",
                  test = TRUE)

test.10x<-sub_h5(barcodes = test.bc,
                  feature_mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv",
                  count_mtx = rbind(test.rna,test.atac),
                  path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5",
                  test = TRUE)


test1.rna <-sub_h5(barcodes = colnames(test1.rna),
                  feature_mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv",
                  count_mtx = test1.rna,
                  path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5",
                  module = "Gene Expression",
                  test = TRUE)
test2.rna <-sub_h5(barcodes = colnames(test2.rna),
                   feature_mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv",
                   count_mtx = test2.rna,
                   path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_rna.h5",
                   module = "Gene Expression",
                   test = TRUE)

atac_test1 <-sub_h5(barcodes = colnames(test1.atac),
                   feature_mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv",
                   count_mtx = test1.atac,
                   path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_atac.h5",
                   module = "Peaks",
                   test = TRUE)
atac_test2 <-sub_h5(barcodes = colnames(test2.atac),
                    feature_mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv",
                    count_mtx = test2.atac,
                    path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5",
                    module = "Peaks",
                    test = TRUE)






## save filter_bc_matrix
## train: multiomic
sub_bc_matrix(barcodes = train.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = rbind(train.rna,train.atac),
              output_dir <- "filtered_feature_bc_matrix_train/"
)

sub_bc_matrix(barcodes = test.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = rbind(test.rna,test.atac),
              output_dir <- "filtered_feature_bc_matrix_test/"
)


## train data seperate rna and atac
## train -rna
sub_bc_matrix(barcodes = train.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = train.rna,
              module = "Gene Expression",
              output_dir <- "filtered_feature_bc_matrix_train_rna/"
)


sub_bc_matrix(barcodes = train.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = train.atac,
              module = "Peaks",
              output_dir <- "filtered_feature_bc_matrix_train_atac/"
)





## test -rna
sub_bc_matrix(barcodes = test1.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = test1.rna,
              module = "Gene Expression",
              output_dir <- "filtered_feature_bc_matrix_test1_rna/"
)



sub_bc_matrix(barcodes = test2.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = test2.rna,
              module = "Gene Expression",
              output_dir <- "filtered_feature_bc_matrix_test2_rna/"
)

sub_bc_matrix(barcodes = test1.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = test1.atac,
              module = "Peaks",
              output_dir <- "filtered_feature_bc_matrix_test1_atac/"
)


sub_bc_matrix(barcodes = test2.bc,
              feature.mtx = "/path/to/data/pbmc10k/filtered_feature_bc_matrix/features.tsv", 
              count.mtx = test2.atac,
              module = "Peaks",
              output_dir <- "filtered_feature_bc_matrix_test2_atac/"
)



#################################################
## check h5 file
out_dir <-'/path/to/data/pbmc3k/'
h5_train<-Read10X_h5(paste0(out_dir,'pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_train.h5'))
h5_test<-Read10X_h5(paste0(out_dir,'pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5'))
h5_test1_rna<-Read10X_h5(paste0(out_dir,'pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5'))
h5_test1_atac<-Read10X_h5(paste0(out_dir,'pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_atac.h5'))
h5_test2_rna<-Read10X_h5(paste0(out_dir,'pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_rna.h5'))
h5_test2_atac<-Read10X_h5(paste0(out_dir,'pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5'))
print(dim(h5_train$`Gene Expression`))
print(dim(h5_train$Peaks))
print(dim(h5_test1_rna))
print(dim(h5_test1_atac))
print(dim(h5_test2_rna))
print(dim(h5_test2_atac))
print(dim(h5_test$`Gene Expression`))

print(length(train.bc))
print(length(test1.bc))
print(length(test2.bc))

## to test stability,  do multiple random sample


##########################################
## annotation file prepared
## original celltype label file
library(dplyr)
celltype_label.file <- "/path/to/data/pbmc3k/pbmc3k_celltype_annotation.csv" 
celltype_label <- read.csv(celltype_label.file)
test.bc<-colnames(h5_test$`Gene Expression`)
train.bc<-colnames(h5_train$`Gene Expression`)
testcelltype_label <- celltype_label %>%
  filter(Barcode %in% test.bc)
test_rna_celltype <- 
  testcelltype_label %>%
  mutate(Barcode = paste0(Barcode,"_rna"),
         modality = "test scRNA")

test_atac_celltype <- 
  testcelltype_label %>%
  mutate(Barcode = paste0(Barcode,"_atac"),
         modality = "test scATAC")

train_celltype <-
  celltype_label %>%
  filter(Barcode %in% train.bc)%>%
  mutate(modality ="train multiomics")


pbmc3k_cellannot <-rbind(train_celltype, test_rna_celltype, test_atac_celltype)
write.csv(pbmc3k_cellannot,paste0(out_dir,"pbmc3k_cellannot_1113.csv"), row.names = FALSE)
