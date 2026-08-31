## interactive on gpu node on hpc
## 
##bsub -P hpcf_interactive -J hpcf_interactive -n 1 -q gpu_short -R "rusage[mem=100001]" -gpu "num=1/host" -Is "bash"
## module load  conda3/202105
## cd  to the test dir
###########################################
import h5py
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

hf = h5py.File('rna.pbmc.h5', 'r')
rna_mtx = hf.get('RNA')
rna_mtx = np.array(rna_mtx).T

rna_meta = pd.read_csv('rna.meta.csv',sep='\t')
rna_gene = pd.read_csv('rna_gene_name.txt',header=None)
rna_gene[0] = rna_gene[0].str.upper()
rna_gene.shape

len(np.unique(rna_gene[0].to_numpy()))

rna_anndat = ad.AnnData(
    X = rna_mtx,
    obs = rna_meta,
)
rna_anndat.var_names = rna_gene[0]

rna_anndat.layers['counts'] = rna_anndat.X.copy()
sc.pp.normalize_total(rna_anndat, target_sum=1e6)
sc.pp.log1p(rna_anndat)

rna_anndat.obs['cell_type'] = rna_anndat.obs['seurat_clusters']

hf = h5py.File('atac.brain.h5', 'r')
atac_mtx = hf.get('ATAC')
atac_mtx = np.array(atac_mtx).T

atac_meta = pd.read_csv('atac.meta.csv',sep='\t')
atac_gene = pd.read_csv('atac_gene_name.txt',header=None)
atac_gene[0] = atac_gene[0].str.upper()

atac_anndat = ad.AnnData(
    X = atac_mtx,
    obs = atac_meta,
)
atac_anndat.var_names = atac_gene[0]

sc.pp.normalize_total(atac_anndat, target_sum=1e6)
sc.pp.log1p(atac_anndat)

rna_anndat.write_h5ad('pbmc10k_10X_GEX.h5ad')
atac_anndat.write_h5ad('pbmc10k_10X_ATAC.h5ad')
