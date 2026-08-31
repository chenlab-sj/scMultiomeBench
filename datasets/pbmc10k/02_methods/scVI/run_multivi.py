# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).

import numpy as np
import scanpy as sc
import scvi

import pandas as pd
import os
import timeit
######################################################
## input args
input_10x="/path/to/data/pbmc10k/filtered_feature_bc_matrix"
out_dir = "res_pbmc10k_control"
########################################################
start = timeit.default_timer()
scvi.settings.seed = 420
adata = scvi.data.read_10x_multiome(input_10x)
adata.var_names_make_unique()
########### test with small amount imputation to avoid an error message
n = 10
adata_rna = adata[:n, adata.var.modality == "Gene Expression"].copy()
adata_paired = adata.copy()
adata_atac = adata[n:2 * n, adata.var.modality == "Peaks"].copy()
adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired, adata_rna, adata_atac)
## order variable
adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()

print(adata_mvi.shape)
sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
print(adata_mvi.shape)

scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="modality")
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
adata_mvi.obs = adata_mvi.obs.set_axis([s. split("_", 1)[0] for s in adata_mvi.obs.index], axis='index')
lat_df = pd.DataFrame(adata_mvi.obsm['MultiVI_latent'],index=adata_mvi.obs.index)

## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
### save latent results    
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
