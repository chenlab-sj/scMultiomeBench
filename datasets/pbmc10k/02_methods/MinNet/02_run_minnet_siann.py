# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import logging
import anndata as ad
import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import sys
import torch
import os
import matplotlib.pyplot as plt
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = '0'
os.environ["CUDA_LAUNCH_BLOCKING"] = '1'

import scanpy as sc

from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from sklearn.metrics import pairwise_distances

from utils.test_util import test_data_multiome, Siamese_Test_multiome
from utils.feature_util import feature_selection_multiome
###############################################################
test_rna='/path/to/tools/MinNet/pbmc10k_control/pbmc10k_10X_GEX.h5ad' 
test_atac='/path/to/tools/MinNet/pbmc10k_control/pbmc10k_10X_ATAC.h5ad'
out_dir='/path/to/tools/MinNet/pbmc10k_test'

###########################################################
input_mod1 = ad.read_h5ad(test_rna)
input_mod2 = ad.read_h5ad(test_atac)

input_mod1, input_mod2 = feature_selection_multiome(input_mod1, input_mod2, path='./utils')

sc.pp.scale(input_mod1, max_value=10)
sc.pp.scale(input_mod2, max_value=10)
input_mod1.X = sp.sparse.csr_matrix(input_mod1.X)
input_mod2.X = sp.sparse.csr_matrix(input_mod2.X)

test_dat = test_data_multiome(mod1=input_mod1,
                              mod2=input_mod2,
                              batch_number=500)
SiaNN = Siamese_Test_multiome(test_dat, num_peaks=6532, num_genes=6532)
ckpt_loader = torch.load('utils/SiaNN.ckpt')
SiaNN.D.load_state_dict(ckpt_loader)

out_x, out_y, class_x, class_y = SiaNN.test()

merge = np.concatenate((out_x, out_y), axis=0)

lat_df = pd.DataFrame(merge,index=np.concatenate((input_mod1.obs_names.to_numpy(), 
                                         input_mod2.obs_names.to_numpy()), axis=0))
lat_df = lat_df.set_axis(["latent_" + s  for s in lat_df.columns.astype("str").tolist()],axis="columns")
latent_csv = os.path.join(out_dir,"coembed_coor.csv")
lat_df.to_csv(latent_csv)