# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import sys, os

import numpy as np
from umap import UMAP
import time
import torch
import matplotlib.pyplot as plt
import pandas as pd  
import scipy.sparse as sp

import scmomat 
import random
#plt.rcParams["font.size"] = 10
import scanpy as sc
from scipy import stats
sys.path.append('/path/to/tools/scMoMaT')
import utils


test_data_path="/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test.h5"
out_dir = '/path/to/tools/scMoMaT/pbmc10k/'

###########################################
#start = timeit.default_timer()
random.seed(420)
test_data=sc.read_10x_h5(test_data_path, gex_only=False)
test_rna=test_data[:,test_data.var['feature_types']=='Gene Expression']
test_rna.var_names_make_unique()
test_atac=test_data[:,test_data.var['feature_types']=='Peaks']
test_atac.obs.index = test_atac.obs.index + '_atac'
test_rna.obs.index = test_rna.obs.index + '_rna'

## find hvg in RNA modality
test_rna.layers["counts"] = test_rna.X.copy()
sc.pp.normalize_total(test_rna, target_sum=1e4)
sc.pp.log1p(test_rna)
sc.pp.highly_variable_genes(test_rna,n_top_genes=7000)
hvg_rna = test_rna.var['highly_variable'].index[test_rna.var['highly_variable']].tolist()


##
A_long = pd.read_csv(os.path.join(out_dir, 'gact.csv'))
A_hvg = A_long[A_long.loc[:, 'gene.name'].isin(hvg_rna)]
A = pd.crosstab(A_hvg.loc[:, 'peak'], A_hvg.loc[:, 'gene.name'])
# filter for peaks that are within 2000bp of TSS or along the gene body, as well as genes that have at least one peak
gene_sel = A.columns.values.squeeze()
peak_sel = A.index.values.squeeze()
counts_rna = test_rna[:,gene_sel].layers['counts'].todense()
counts_rna = utils.preprocess(counts_rna, modality = "RNA", log = False) 

## change from scDART prep
peak_sel = np.array([s.replace('_', ':', 1).replace('_', '-', 1) for s in peak_sel], dtype=object)

counts_atac = test_atac[:,peak_sel].X.todense()
counts_atac = utils.preprocess(counts_atac, modality = "ATAC") 
counts_atac_pseudoRNA = counts_atac @ A
counts_atac_pseudoRNA = np.asmatrix(counts_atac_pseudoRNA.to_numpy())

counts = {"rna":[counts_rna,counts_atac_pseudoRNA], "atac": [None,counts_atac] }
feats_name = {"rna": gene_sel, "atac": peak_sel}
counts["nbatches"] = 2
######################################
## run scMoMaT model
##########################################
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
lamb = 0.001
batchsize = 0.1
# running seed
seed = 0
# number of latent dimensions
K = 30
interval = 1000
T = 10000
lr = 1e-2
#1st stage training, learning cell factors
model = scmomat.scmomat_model(counts = counts, K = K, device = device)
losses = model.train_func(T = T)


print("------ Done ------")
#stop = timeit.default_timer()

# extract cell factors/latent representations
zs = model.extract_cell_factors()
    
cell_bc = np.concatenate([test_rna.obs.index.values.squeeze(),
                          test_atac.obs.index.values.squeeze()])
    
res_df = pd.DataFrame(np.concatenate(zs, axis=0),
                          index=cell_bc)

# set column names as latent_x 
res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")


## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

res_df.to_csv(os.path.join(out_dir,"latent.csv"))

#print('Time(s): ', stop - start)  
# record time 
#runtime_out = os.path.join(out_dir,"runtime.txt")
#with open(runtime_out, 'w') as file:
#    file.write(str(stop-start))
