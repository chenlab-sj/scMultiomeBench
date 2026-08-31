# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
from unioncom import UnionCom
import numpy as np
import scanpy as sc
import os
import torch
import pandas as pd
test_data_path="/path/to/data/pbmc3k/"
out_dir = '/path/to/tools/unioncom/pbmc3k/'

###########################################
#start = timeit.default_timer()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


test_rna=sc.read_h5ad(os.path.join(test_data_path,"pbmc3k_testrna_embed.h5ad"))
test_atac=sc.read_h5ad(os.path.join(test_data_path,"pbmc3k_testatac_embed.h5ad"))


embed_rna = test_rna.X
embed_atac = test_atac.X

uc = UnionCom.UnionCom(epoch_pd = 50000)
integrated_data = uc.fit_transform([embed_rna, embed_atac])
z_rna = integrated_data[0]
z_atac = integrated_data[1]

lat_rna = pd.DataFrame(data = z_rna, index = test_rna.obs.index)
lat_atac = pd.DataFrame(data = z_atac, index = test_atac.obs.index)
lat_df = pd.concat([lat_rna, lat_atac], axis=0)
lat_df = lat_df.set_axis(["latent_" + s  for s in lat_df.columns.astype("str").tolist()],axis="columns")
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

lat_df.to_csv(os.path.join(out_dir,"latent.csv"))
