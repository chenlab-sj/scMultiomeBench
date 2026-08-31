# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import os
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import csv
import gzip
import scipy.io

os.environ["CUDA_VISIBLE_DEVICES"] = "3"

np.random.seed(0)

sc.settings.verbosity = 3
sc.logging.print_header()

out_dir='/path/to/tools/portal/Portal/pbmc3k/rep2/'
#test_data_path="/path/to/tools/scJoint/pbmc3k/scJoint_need/"
test_h5="/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
atac_gene_h5="/path/to/tools/multimap/pbmc3k/atac_gene.h5"

test_data=sc.read_10x_h5(test_h5, gex_only=False)
rna=test_data[:,test_data.var['feature_types']=='Gene Expression']
#atac_peaks=test_data[:,test_data.var['feature_types']=='Peaks']
atac_genes=sc.read_10x_h5(atac_gene_h5, gex_only = False)

rna.var_names_make_unique()
#atac_peaks.var_names_make_unique()
atac_genes.var_names_make_unique()

import portal
# Specify the GPU device
#os.environ["CUDA_VISIBLE_DEVICES"] = "6"

# Create a folder for saving results

model = portal.model.Model(training_steps=3000, lambdacos=10.0,seed=0)
model.preprocess(rna, atac_genes)
model.train()
model.eval() 
rna_index = [obs_name + "_rna" for obs_name in rna.obs.index]
atac_genes_index = [obs_name + "_atac" for obs_name in atac_genes.obs.index]
bc= rna_index + atac_genes_index
lat_df = pd.DataFrame(model.latent, index = bc )
lat_df = lat_df.set_axis(["latent_" + s  for s in lat_df.columns.astype("str").tolist()],axis="columns")


## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
### pair latent results    
latent_csv = os.path.join(out_dir,"lat_df.csv")
lat_df.to_csv(latent_csv)
