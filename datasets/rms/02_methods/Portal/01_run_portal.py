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

# Reproducibility replicates: rep1 = base run (SEED=420, REP=""), rep2 (SEED=7), rep3 (SEED=42).
SEED = int(os.environ.get("SEED", "420"))
REP = os.environ.get("REP", "")

np.random.seed(SEED)

sc.settings.verbosity = 3
sc.logging.print_header()

out_dir = '/path/to/data/RMS/Mast607A_TB19_22652/portal/'
test_data_path="/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5"
#test_atac_data_path="/path/to/data/RMS/SJRHB013758_X2_scATAC/SJRHB013758_X2_scATAC_filtered_feature_bc_matrix.h5"
atac_gene_h5="/path/to/tools/multimap/RMS/Mast607A_TB19_22652/atac_gene.h5"

test_data=sc.read_10x_h5(test_data_path, gex_only=False)
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

model = portal.model.Model(training_steps=8000, lambdacos=10.0, seed=SEED)
model.preprocess(rna, atac_genes)
model.train()
model.eval() 
rna_index = [obs_name + "_rna" for obs_name in rna.obs.index]
atac_genes_index = [obs_name + "_atac" for obs_name in atac_genes.obs.index]
bc= rna_index + atac_genes_index
lat_df = pd.DataFrame(model.latent, index = bc )
lat_df = lat_df.set_axis(["latent_" + s  for s in lat_df.columns.astype("str").tolist()],axis="columns")


## save results:
out = os.path.join(out_dir, REP)
os.makedirs(out, exist_ok=True)
### pair latent results
latent_csv = os.path.join(out,"lat_df.csv")
lat_df.to_csv(latent_csv)
