# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import os
import simba as si
import pandas as pd
import scanpy as sc

test_data_path="/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
out_dir = '/path/to/tools/simba/pbmc3k/'

###########################################
#start = timeit.default_timer()
test_data=sc.read_10x_h5(test_data_path, gex_only=False)
test_rna=test_data[:,test_data.var['feature_types']=='Gene Expression']
test_rna.var_names_make_unique()
test_atac=test_data[:,test_data.var['feature_types']=='Peaks']
test_atac.obs.index = test_atac.obs.index + '_atac'
test_rna.obs.index = test_rna.obs.index + '_rna'

## ATAC processing
si.pp.filter_peaks(test_atac,min_n_cells=3)
si.pp.cal_qc_atac(test_atac)

si.pp.pca(test_atac, n_components=50)
si.pp.select_pcs(test_atac,n_pcs=40)
si.pp.select_pcs_features(test_atac)

## RNA process
si.pp.filter_genes(test_rna,min_n_cells=3)
si.pp.cal_qc_rna(test_rna)
si.pp.normalize(test_rna,method='lib_size')
si.pp.log_transform(test_rna)
si.pp.select_variable_genes(test_rna, n_top_genes=4000)
si.tl.discretize(test_rna,n_bins=5) #discretize RNA expression
test_atac.var[['chr', 'start', 'end']] = test_atac.var['gene_ids'].str.extract(r'([^:]+):(\d+)-(\d+)')
test_rna_atac = si.tl.gene_scores(test_atac,genome='hg38',use_gene_weigt=True, use_top_pcs=True)

si.pp.filter_genes(test_rna_atac,min_n_cells=3)
si.pp.cal_qc_rna(test_rna_atac)
si.pp.normalize(test_rna_atac,method='lib_size')
si.pp.log_transform(test_rna_atac)
adata_CrnaCatac = si.tl.infer_edges(test_rna, test_rna_atac, n_components=15, k=15)

# edges can be futhere trimmed if needed. Here we keep all of them
si.tl.trim_edges(adata_CrnaCatac, cutoff=0.5)
## generate graph
si.tl.gen_graph(list_CP=[test_atac],
                list_CG=[test_rna],
                list_CC=[adata_CrnaCatac],
                copy=False,
                use_highly_variable=True,
                use_top_pcs=True,
                dirname='graph0')

print(si.settings.pbg_params)
#si.tl.pbg_train(auto_wd=True, save_wd=True, output='model')

# modify parameters
dict_config = si.settings.pbg_params.copy()
# # dict_config['wd'] = 0.000282
dict_config['workers'] = 1
## start training
si.tl.pbg_train(pbg_params = dict_config, auto_wd=True, save_wd=True, output='model')
# load in graph ('graph0') info
si.load_graph_stats(path='./result_simba/pbg/graph0')
# load in model info for ('graph0')
si.load_pbg_config(path='./result_simba/pbg/graph0/model/')

## post training analysis
dict_adata = si.read_embedding()
adata_C = dict_adata['C']  # embeddings for ATACseq cells
adata_C2 = dict_adata['C2']  # embeddings for RNAseq cells
adata_G = dict_adata['G']  # embeddings for genes
adata_P = dict_adata['P']  # embeddings for peaks

adata_all = si.tl.embed(adata_ref=adata_C2,
                        list_adata_query=[adata_C],
                        use_precomputed=False)

lat_df = pd.DataFrame(adata_all.X, index=adata_all.obs_names, columns=adata_all.var_names)
lat_df = lat_df.set_axis(["latent_" + s  for s in lat_df.columns.astype("str").tolist()],axis="columns")

## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

lat_df.to_csv(os.path.join(out_dir,"latent.csv"))
