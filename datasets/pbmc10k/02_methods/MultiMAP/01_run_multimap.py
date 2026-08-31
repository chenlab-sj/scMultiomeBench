# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import scanpy as sc
import anndata
import MultiMAP
import pandas as pd
import os

## import data
#rna_h5="/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5"
#atac_peak_h5="/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5"
atac_gene_h5="/path/to/tools/multimap/pbmc10k/multimap_need/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_control_atac_gene.h5"
test_h5="/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5"
out_dir="res_pbmc10k_control"


test_data=sc.read_10x_h5(test_h5, gex_only=False)
rna=test_data[:,test_data.var['feature_types']=='Gene Expression']
atac_peaks=test_data[:,test_data.var['feature_types']=='Peaks']
atac_genes=sc.read_10x_h5(atac_gene_h5, gex_only = False)

rna.var_names_make_unique()
atac_peaks.var_names_make_unique()
atac_genes.var_names_make_unique()

## add suffix to rna data and atac data
rna_index = [obs_name + "_rna" for obs_name in rna.obs.index]
rna.obs.index = rna_index
atac_peaks_index = [obs_name + "_atac" for obs_name in atac_peaks.obs.index]
atac_peaks.obs.index = atac_peaks_index
atac_genes_index = [obs_name + "_atac" for obs_name in atac_genes.obs.index]
atac_genes.obs.index = atac_genes_index
####################################
## start multimap
MultiMAP.TFIDF_LSI(atac_peaks)
atac_genes.obsm['X_lsi'] = atac_peaks.obsm['X_lsi'].copy()
rna_pca = rna.copy()
sc.pp.scale(rna_pca)
sc.pp.pca(rna_pca)
rna.obsm['X_pca'] = rna_pca.obsm['X_pca'].copy()
adata = MultiMAP.Integration([rna, atac_genes], ['X_pca', 'X_lsi'])
adata2 = MultiMAP.Integration([rna, atac_genes], ['X_pca', 'X_lsi'], embedding=False)

lat_df = pd.DataFrame(adata.obsm['X_multimap'],index=adata.obs.index)
lat_df = lat_df.set_axis(["latent_" + s  for s in lat_df.columns.astype("str").tolist()],axis="columns")

## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
### pair latent results    
latent_csv = os.path.join(out_dir,"coembed_coor.csv")
lat_df.to_csv(latent_csv)
out_anndata = os.path.join(out_dir,"connect.h5ad")
adata2.write(out_anndata)