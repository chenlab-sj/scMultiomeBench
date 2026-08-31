# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import numpy as np
import scanpy as sc
import scvi

import pandas as pd
import os
import timeit

import anndata
######################################################
## https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/MultiVI_tutorial.html
## input args
train_10x = "/path/to/data/BMMC_d1/s1d2_train" ## train data with combined peak
test1_10x="/path/to/data/BMMC_d1/train_test1_s2d1" ## test rna
test2_10x = "/path/to/data/BMMC_d1/train_test2_s4d1" ## test rna
test3_10x = "/path/to/data/BMMC_d1/train_test3_s1d1" ## test rna
out_dir = "/path/to/multiomeBench/BMMC_d1/s1d1_paired/scVI/"
########################################################
start = timeit.default_timer()

scvi.settings.seed = 420
adata_train = scvi.data.read_10x_multiome(train_10x)
adata_train.var_names_make_unique()

adata_test1 = scvi.data.read_10x_multiome(test1_10x)
adata_test1.var_names_make_unique()
adata_test2 = scvi.data.read_10x_multiome(test2_10x)
adata_test2.var_names_make_unique()
adata_test3 = scvi.data.read_10x_multiome(test3_10x)
adata_test3.var_names_make_unique()

adata_test1_rna =adata_test1[:,adata_test1.var.modality == "Gene Expression"].copy()
adata_test1_atac = adata_test1[:,adata_test1.var.modality == "Peaks"].copy()

adata_test2_rna =adata_test2[:,adata_test2.var.modality == "Gene Expression"].copy()
adata_test2_atac = adata_test2[:,adata_test2.var.modality == "Peaks"].copy()

adata_test3_rna =adata_test3[:,adata_test3.var.modality == "Gene Expression"].copy()
adata_test3_atac = adata_test3[:,adata_test3.var.modality == "Peaks"].copy()

## rename bc and add data batch 
#adata_train.obs.index = adata_train.obs.index +"_train"
adata_train.obs["data"]="train"

adata_test1_rna.obs.index = adata_test1_rna.obs.index + "_rna1"
adata_test1_rna.obs["data"]="test1"
adata_test2_rna.obs.index = adata_test2_rna.obs.index + "_rna2"
adata_test2_rna.obs["data"]="test2"
adata_test3_rna.obs.index = adata_test3_rna.obs.index + "_rna3"
adata_test3_rna.obs["data"]="test3"

adata_test1_atac.obs.index = adata_test1_atac.obs.index + "_atac1"
adata_test1_atac.obs["data"]="test1"
adata_test2_atac.obs.index = adata_test2_atac.obs.index + "_atac2"
adata_test2_atac.obs["data"]="test2"
adata_test3_atac.obs.index = adata_test3_atac.obs.index + "_atac3"
adata_test3_atac.obs["data"]="test3"

adata_rna = adata_test3_rna.copy()
adata_atac = adata_test3_atac.copy()

adata_mvi = scvi.data.organize_multiome_anndatas(adata_train, adata_rna, adata_atac)
## order variable
adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()

print(adata_mvi.shape)
sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
print(adata_mvi.shape)

scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="modality")
#scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="modality",categorical_covariate_keys=["data"])

mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=(adata_mvi.var["modality"] == "Gene Expression").sum(),
    n_regions=(adata_mvi.var["modality"] == "Peaks").sum(),
)
mvi.view_anndata_setup()
mvi.train()
#mvi.save("trained_multivi")
#mvi = scvi.model.MULTIVI.load("trained_multivi", adata=adata_mvi)
adata_mvi.obsm["MultiVI_latent"] = mvi.get_latent_representation()
stop = timeit.default_timer()

## prepare latent dataframe
#adata_mvi.obs = adata_mvi.obs.set_axis([s. split("_", 1)[0] for s in adata_mvi.obs.index], axis='index')
lat_df = pd.DataFrame(adata_mvi.obsm['MultiVI_latent'],index=adata_mvi.obs.index)

## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

## cell embeddings 
latent_csv = os.path.join(out_dir,"latent.csv")
lat_df.to_csv(latent_csv)

### multiVI model
mvi_model = os.path.join(out_dir,"trained_multivi")
mvi.save(mvi_model)

print('Time(s): ', stop - start)  
# record time 
runtime_out = os.path.join(out_dir,"multivi_runtime.txt")
with open(runtime_out, 'w') as file:
    file.write(str(stop-start)) 


## save imputed information
imputed_expression = mvi.get_normalized_expression()
imputed_expression_adata=sc.AnnData(imputed_expression)
imputed_rna_file = os.path.join(out_dir,"imputed_rna.h5ad")
imputed_expression_adata.write_h5ad(imputed_rna_file)


imputed_access = mvi.get_accessibility_estimates()
imputed_access_adata=sc.AnnData(imputed_access)
imputed_atac_file = os.path.join(out_dir,"imputed_atac.h5ad")
imputed_access_adata.write_h5ad(imputed_atac_file)

