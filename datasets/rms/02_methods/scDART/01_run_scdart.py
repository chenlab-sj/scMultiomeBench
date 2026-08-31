# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import sys
sys.path.append('../')

import numpy as np
import pandas as pd

import torch
from sklearn.decomposition import PCA

import scDART.utils as utils
import scDART.TI as ti
import scDART

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

import scipy.io
import scanpy as sc
import timeit
import random
import os

############ input parm
test_data_path="/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5"
region2gene = '/path/to/data/RMS/Mast607A_TB19_22652/scDART/region2gene.mtx'

gact_peak=  '/path/to/data/RMS/Mast607A_TB19_22652/scDART/gact_peak.csv'
gact_gene=  '/path/to/data/RMS/Mast607A_TB19_22652/scDART/gact_gene.csv'

# reproducibility reps: SEED + REP env. rep1 = the EXISTING latent.csv (untouched); reps_sub runs
# REP=rep2 SEED=7 and REP=rep3 SEED=42 -> scDART/rep2,rep3/latent.csv (7/42 avoid the original seed 0).
SEED = int(os.environ.get("SEED", "420"))
REP  = os.environ.get("REP", "")
out_dir = os.path.join('/path/to/multiomeBench/RMS/Mast607/script/scDART/', REP)
os.makedirs(out_dir, exist_ok=True)
###########################################
start = timeit.default_timer()
random.seed(SEED)

### train model
# all in one: from scDART example
seeds = [SEED]
latent_dim = 10
learning_rate = 3e-4
n_epochs = 500
use_anchor = False
reg_d = 1
reg_g = 1
reg_mmd = 1
ts = [30, 50, 70]
use_potential = True

test_data=sc.read_10x_h5(test_data_path, gex_only=False)
test_rna=test_data[:,test_data.var['feature_types']=='Gene Expression']
test_atac=test_data[:,test_data.var['feature_types']=='Peaks']

test_rna.var_names_make_unique()
#atac_peaks.var_names_make_unique()
test_atac.var_names_make_unique()

test_atac.obs.index = test_atac.obs.index + '_atac'
test_rna.obs.index = test_rna.obs.index + '_rna'



## find hvg in counts_rna
## find hvg in RNA modality
test_rna.layers["counts"] = test_rna.X.copy()
sc.pp.normalize_total(test_rna, target_sum=1e4)
sc.pp.log1p(test_rna)
sc.pp.highly_variable_genes(test_rna,n_top_genes=1000)
hvg_rna = test_rna.var['highly_variable'].index[test_rna.var['highly_variable']].tolist()

##
peak = pd.read_csv(gact_peak)
peak=peak['x'].tolist()
gene = pd.read_csv(gact_gene)
gene=gene['x'].tolist()
coarse_reg_mtx = scipy.io.mmread(region2gene)
coarse_reg = coarse_reg_mtx.toarray()
gene_mask =np.isin(gene, hvg_rna)
gene_sel = np.array(gene)[gene_mask]
gene_sel = gene_sel.tolist()
coarse_reg_sel = coarse_reg[:,np.isin(gene, hvg_rna)]


counts_rna = test_rna[:,gene_sel].to_df()
counts_atac= test_atac.to_df()
counts_atac= counts_atac[peak]

counts_rna.index = test_rna.obs.index
counts_atac.index = test_atac.obs.index


############
## training model
scDART_op = scDART.scDART(n_epochs = n_epochs, latent_dim = latent_dim, \
        ts = ts, use_anchor = use_anchor, use_potential = use_potential, k = 10, \
        reg_d = 1, reg_g = 1, reg_mmd = 1, l_dist_type = 'kl', seed = seeds[0],\
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu'))

scDART_op = scDART_op.fit(rna_count = counts_rna.values, atac_count = counts_atac.values, reg = coarse_reg_sel, rna_anchor = None, atac_anchor = None)
z_rna, z_atac = scDART_op.transform(rna_count = counts_rna.values, atac_count = counts_atac.values, rna_anchor = None, atac_anchor = None)

print("------ Done ------")
stop = timeit.default_timer()
## export latent results
z_rna_df = pd.DataFrame(z_rna, index = counts_rna.index)
z_atac_df = pd.DataFrame(z_atac, index = counts_atac.index)
df_latent= pd.concat([z_rna_df, z_atac_df])
df_latent = df_latent.set_axis(["latent_" + s  for s in df_latent.columns.astype("str").tolist()],axis="columns")
## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

df_latent.to_csv(os.path.join(out_dir,"latent.csv"))

print('Time(s): ', stop - start)  
# record time 
runtime_out = os.path.join(out_dir,"runtime.txt")
with open(runtime_out, 'w') as file:
    file.write(str(stop-start))
