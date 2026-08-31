# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import anndata
import pandas as pd
import os
import scanpy as sc

outdir ="/path/to/tools/scBridge/pbmc10k/"
cell_annot = "/path/to/data/pbmc10k/pbmc10k_cellannot_1113.csv"

rna_h5 = os.path.join(outdir, "rna_comgene.h5")
atac_h5 = os.path.join(outdir, "atac_comgene.h5")

rna = sc.read_10x_h5(rna_h5, gex_only=False)
atac = sc.read_10x_h5(atac_h5, gex_only=False)

rna.obs_names = [name + "_rna" for name in rna.obs_names]
atac.obs_names = [name + "_atac" for name in atac.obs_names]


labels_annot= pd.read_csv(cell_annot,index_col=0)
rna_label=labels_annot.loc[rna.obs_names]
rna.obs['CellType']=rna_label['cell_type'].fillna(value="Unknown")
#atac_peaks_label=labels_annot.loc[atac_peaks.obs_names]
#atac_peaks.obs['CellType']=atac_peaks_label['celltype']

rna.write_h5ad("scBridge_rna.h5ad")
atac.write_h5ad("scBridge_atac.h5ad")