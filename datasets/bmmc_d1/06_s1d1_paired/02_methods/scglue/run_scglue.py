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

## intake file
test1_h5="/path/to/data/BMMC_d1/test1_s2d1.h5"
test2_h5="/path/to/data/BMMC_d1/test2_s4d1.h5"
test3_h5="/path/to/data/BMMC_d1/test3_s1d1.h5"
out_dir = '/path/to/multiomeBench/BMMC_d1/s1d1_paired/scglue/'

start = timeit.default_timer()
random.seed(420)
## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

## prepocess rna

## load test1 data
test1=sc.read_10x_h5(test1_h5, gex_only=False)
rna1=test1[:,test1.var['feature_types']=='Gene Expression']
atac1=test1[:,test1.var['feature_types']=='Peaks']
rna1.var_names_make_unique()

rna1_index = [obs_name + "_rna1" for obs_name in rna1.obs.index]
rna1.obs.index = rna1_index
atac1_index = [obs_name + "_atac1" for obs_name in atac1.obs.index]
atac1.obs.index = atac1_index

rna1.obs["data"]="test"
atac1.obs["data"]="test"

## load test2 data
test2=sc.read_10x_h5(test2_h5, gex_only=False)
rna2=test2[:,test2.var['feature_types']=='Gene Expression']
atac2=test2[:,test2.var['feature_types']=='Peaks']
rna2.var_names_make_unique()

rna2_index = [obs_name + "_rna2" for obs_name in rna2.obs.index]
rna2.obs.index = rna2_index
atac2_index = [obs_name + "_atac2" for obs_name in atac2.obs.index]
atac2.obs.index = atac2_index

rna2.obs["data"]="test"
atac2.obs["data"]="test"

## load test3 data
test3=sc.read_10x_h5(test3_h5, gex_only=False)
rna3=test3[:,test3.var['feature_types']=='Gene Expression']
atac3=test3[:,test3.var['feature_types']=='Peaks']
rna3.var_names_make_unique()

rna3_index = [obs_name + "_rna3" for obs_name in rna3.obs.index]
rna3.obs.index = rna3_index
atac3_index = [obs_name + "_atac3" for obs_name in atac3.obs.index]
atac3.obs.index = atac3_index

rna3.obs["data"]="test"
atac3.obs["data"]="test"

rna = ad.concat([rna3], axis = 0)
rna.var = rna.var
rna.layers["counts"] =rna.X.copy()
#########################################
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")

sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")

## prepocss atac
atac = ad.concat([atac3], axis = 0)
atac.var = atac3.var
scglue.data.lsi(atac, n_components=100, n_iter=15)

## construct prior regulat graph
##  wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
scglue.data.get_gene_annotation(
    rna, gtf="/path/to/tools/scglue/gencode.v32.primary_assembly.annotation.gtf.gz",
    gtf_by="gene_name"
)
#rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
#atac.var.head()

######### fix rna issue caused by duplicate var
ann_genes_list = rna.var[~rna.var.loc[:, ["chrom"]].isnull().iloc[:,0]].index.values.tolist()
#Extract annotated genes 
rna = rna[:,ann_genes_list] 
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
########################

## graph construction
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
#guidance

scglue.graph.check_graph(guidance, [rna, atac])

## save processed file 
rna.write(os.path.join(out_dir,"rna-pp.h5ad"), compression="gzip")
atac.write(os.path.join(out_dir,"atac-pp.h5ad"), compression="gzip")
nx.write_graphml(guidance, os.path.join(out_dir,"guidance.graphml.gz"))

## read preprocessed file
#rna = ad.read_h5ad("rna-pp.h5ad")
#atac = ad.read_h5ad("atac-pp.h5ad")
#guidance = nx.read_graphml("guidance.graphml.gz")

# configure
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_batch="data",
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_batch="data",
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

print('Time(s): ', stop - start)  
# record time 
runtime_out = os.path.join(out_dir,"runtime.txt")
with open(runtime_out, 'w') as file:
    file.write(str(stop-start))
    
## store cell embeddings
df_latent.to_csv(os.path.join(out_dir,"latent.csv"))

