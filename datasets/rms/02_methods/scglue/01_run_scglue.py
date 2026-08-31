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
#rna_h5="/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test1_rna.h5"
#atac_h5 = "/path/to/data/pbmc10k/pbmc10k_data_design0719/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix_test2_atac.h5"
test_h5="/path/to/data/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5"
out_dir = '/path/to/tools/scglue/RMS/Mast607A_TB19_22652/'

## reproducibility replicates: rep1 = the existing untouched latent (out_dir/latent.csv).
## rep2 (SEED=7) and rep3 (SEED=42) write under out_dir/<REP>/ so reps never clobber each other.
SEED = int(os.environ.get("SEED", "420"))
REP = os.environ.get("REP", "")
out_dir = os.path.join(out_dir, REP)

start = timeit.default_timer()
random.seed(SEED)
## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

## prepocess rna
test_data=sc.read_10x_h5(test_h5, gex_only=False)
rna=test_data[:,test_data.var['feature_types']=='Gene Expression']
atac=test_data[:,test_data.var['feature_types']=='Peaks']
#rna = sc.read_10x_h5(rna_h5)
rna.var_names_make_unique()
#atac = sc.read_10x_h5(atac_h5, gex_only=False)

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

## prepocss atac

scglue.data.lsi(atac, n_components=100, n_iter=15)

## construct prior regulat graph
##  wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
scglue.data.get_gene_annotation(
    rna, gtf="/path/to/tools/scglue/Homo_sapiens.GRCh37.82.gtf.gz",
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
    init_kws={"random_seed": SEED},
    fit_kws={"directory": os.path.join(out_dir, "glue")}
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

