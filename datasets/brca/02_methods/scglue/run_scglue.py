# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import os
import timeit
import random 

rna_h5 = "/path/to/data/HTAN/HT243B1-S1H4/HT243_S1H4_rna.h5ad"
atac_h5 = "/path/to/data/HTAN/HT243B1-S1H4/HT243_S1H4_atac.h5ad"

out_dir = '/path/to/data/HTAN/HT243B1-S1H4/scglue/'
rna = ad.read_h5ad(rna_h5)
rna.var_names_make_unique()
atac = ad.read_h5ad(atac_h5)

## add suffix to rna data and atac data
rna_index = [obs_name + "_rna" for obs_name in rna.obs.index]
rna.obs.index = rna_index
atac_index = [obs_name + "_atac" for obs_name in atac.obs.index]
atac.obs.index = atac_index

rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")

scglue.data.lsi(atac, n_components=100, n_iter=15)

scglue.data.get_gene_annotation(
    rna, gtf="/path/to/tools/scglue/gencode.v32.primary_assembly.annotation.gtf.gz",
    gtf_by="gene_name"
)

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)


######### fix rna issue caused by duplicate var
ann_genes_list = rna.var[~rna.var.loc[:, ["chrom"]].isnull().iloc[:,0]].index.values.tolist()
#Extract annotated genes 
rna = rna[:,ann_genes_list] 
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
########################

## graph construction
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)

## save processed file 
rna.write(os.path.join(out_dir,"rna-pp.h5ad"), compression="gzip")
atac.write(os.path.join(out_dir,"atac-pp.h5ad"), compression="gzip")
nx.write_graphml(guidance, os.path.join(out_dir,"guidance.graphml.gz"))


# configure
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

## train glue module
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"}
)

## save the module
glue.save(os.path.join(out_dir,"glue.dill"))
#glue = scglue.models.load_model("glue.dill")
print("------ Done ------")
stop = timeit.default_timer()

## Check integration diagnostics
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
dx.to_csv(os.path.join(out_dir,"integration_dx.csv"))

## apply for cell/feature embeddings
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

combined = ad.concat([rna, atac])
df_latent = combined.obsm["X_glue"]
df_latent= pd.DataFrame(combined.obsm["X_glue"], index = combined.obs.index)
df_latent = df_latent.set_axis(["latent_" + s  for s in df_latent.columns.astype("str").tolist()],axis="columns")

#print('Time(s): ', stop - start)  
# record time 
#runtime_out = os.path.join(out_dir,"runtime.txt")
#with open(runtime_out, 'w') as file:
#    file.write(str(stop-start))
    
## store cell embeddings
df_latent.to_csv(os.path.join(out_dir,"latent.csv"))
