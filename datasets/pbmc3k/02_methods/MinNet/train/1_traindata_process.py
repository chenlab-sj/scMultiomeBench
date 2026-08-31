# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## MiniNet train data process

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp

from matplotlib import colors

import os
import anndata

#####################################################################
## input file
datapath ="/path/to/tools/MinNet/pbmc3k/train"
train_rna_norm = os.path.join(datapath, "train_rna_norm.h5")
train_atac_gene_h5 =os.path.join(datapath, "train_atac_gene.h5")
gene_index_file = os.path.join(datapath, "gene_index.txt")
cell_annot = "/path/to/data/pbmc3k/pbmc3k_cellannot_1113.csv"
#########################################################################

input_mod1 = sc.read_10x_h5(train_rna_norm)
input_mod2_gene = sc.read_10x_h5(train_atac_gene_h5)

## select highly-variable genes and overlapped with ATAC activity score genes
gene_index = open(gene_index_file,'r').readlines()
gene_index = [i.strip() for i in gene_index]
gene_index = np.array(gene_index)

input_mod1 = input_mod1[:,gene_index]
input_mod2_gene = input_mod2_gene[:,gene_index]

##*** MinNet require cell type label in obs for training
labels_annot= pd.read_csv(cell_annot,index_col=0)
train_bc = input_mod1.obs.index
labels_annot = labels_annot[labels_annot.index.isin(train_bc)]
labels_annot = labels_annot.reindex(train_bc)

labels_annot['cell_type']=labels_annot['cell_type'].fillna("other")
input_mod1.obs['cell_type']=labels_annot['cell_type']
input_mod2_gene.obs['cell_type']=labels_annot['cell_type']


#This is not the final data for GEX, but is needed next. It will be overwritten later.
input_mod1.write(os.path.join(datapath, "rna_training.h5"))
#This is the final training data for ATAC, save it in data/
input_mod2_gene.write(os.path.join(datapath, "atac_training.h5"))

#########################################################################
### Define KNN Graph and shortest distance between all cell pairs, need for training
input_train_mod1 = ad.read_h5ad(os.path.join(datapath, "rna_training.h5"))
input_mod1 = input_train_mod1.copy()

sc.pp.highly_variable_genes(input_mod1)
input_mod1 = input_mod1[:,input_mod1.var.highly_variable]
sc.pp.scale(input_mod1, max_value=10)
sc.tl.pca(input_mod1, svd_solver='arpack')
sc.pp.neighbors(input_mod1, n_pcs = 50, n_neighbors = 20)
sc.tl.umap(input_mod1)

input_mod1.obsp['connectivities'][input_mod1.obsp['connectivities'].nonzero()] = 1.01 - input_mod1.obsp['connectivities'][input_mod1.obsp['connectivities'].nonzero()]

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
import time

start_time = time.time()
dist_matrix = dijkstra(csgraph=input_mod1.obsp['connectivities'], directed=False, return_predecessors=False)
print("--- %s seconds ---" % (time.time() - start_time))

input_train_mod1.obsp['distance'] = dist_matrix
input_train_mod1.write(os.path.join(datapath, "rna_training.h5"))