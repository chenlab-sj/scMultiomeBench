# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import os
import simba as si
import pandas as pd
import scanpy as sc


def main():
    test_data_path = "/path/to/multiomeBench/pbmc_parse/input_azimuth/"
    rna_h5 = "pbmc_parsed1_rna.h5ad"
    atac_h5 = "pbmc3k_test_atac.h5ad"
    out_dir = "/path/to/multiomeBench/pbmc_parse/simba/"

    test_rna = sc.read_h5ad(os.path.join(test_data_path, rna_h5))
    test_atac = sc.read_h5ad(os.path.join(test_data_path, atac_h5))

    test_rna.var_names_make_unique()
    test_atac.var_names_make_unique()

    test_atac.obs.index = test_atac.obs.index + "_atac"
    test_rna.obs.index = test_rna.obs.index + "_rna"

    ## ATAC processing
    si.pp.filter_peaks(test_atac, min_n_cells=3)
    si.pp.cal_qc_atac(test_atac)
    si.pp.pca(test_atac, n_components=50)
    si.pp.select_pcs(test_atac, n_pcs=40)
    si.pp.select_pcs_features(test_atac)

    ## RNA process
    si.pp.filter_genes(test_rna, min_n_cells=3)
    si.pp.cal_qc_rna(test_rna)
    si.pp.normalize(test_rna, method="lib_size")
    si.pp.log_transform(test_rna)
    si.pp.select_variable_genes(test_rna, n_top_genes=4000)
    si.tl.discretize(test_rna, n_bins=5)  # discretize RNA expression
    test_atac.var["peak"] = test_atac.var.index
    test_atac.var[["chr", "start", "end"]] = test_atac.var["peak"].str.extract(r"^(?P<chr>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    test_rna_atac = si.tl.gene_scores(test_atac, genome="hg38", use_gene_weigt=True, use_top_pcs=True)

    si.pp.filter_genes(test_rna_atac, min_n_cells=3)
    si.pp.cal_qc_rna(test_rna_atac)
    si.pp.normalize(test_rna_atac, method="lib_size")
    si.pp.log_transform(test_rna_atac)
    adata_CrnaCatac = si.tl.infer_edges(test_rna, test_rna_atac, n_components=15, k=15)

    # edges can be further trimmed if needed. Here we keep all of them
    si.tl.trim_edges(adata_CrnaCatac, cutoff=0.5)
    ## generate graph
    si.tl.gen_graph(list_CP=[test_atac],
                    list_CG=[test_rna],
                    list_CC=[adata_CrnaCatac],
                    copy=False,
                    use_highly_variable=True,
                    use_top_pcs=True,
                    dirname="graph0")

    print(si.settings.pbg_params)

    # modify parameters
    dict_config = si.settings.pbg_params.copy()
    dict_config["workers"] = 1
    # NB: do NOT set eval_fraction=0 -- torchbiggraph's _maybe_write_checkpoint assumes the eval
    # ran and references eval_stats_chunk_avg, so disabling eval -> UnboundLocalError. The eval
    # is safe now thanks to the __main__ guard (its spawn worker no longer re-runs the pipeline).
    ## start training
    si.tl.pbg_train(pbg_params=dict_config, auto_wd=True, save_wd=True, output="model")

    ## post training analysis
    dict_adata = si.read_embedding()
    adata_C = dict_adata["C"]    # embeddings for ATACseq cells
    adata_C2 = dict_adata["C2"]  # embeddings for RNAseq cells
    adata_G = dict_adata["G"]    # embeddings for genes
    adata_P = dict_adata["P"]    # embeddings for peaks

    adata_all = si.tl.embed(adata_ref=adata_C2,
                            list_adata_query=[adata_C],
                            use_precomputed=False)

    lat_df = pd.DataFrame(adata_all.X, index=adata_all.obs_names, columns=adata_all.var_names)
    lat_df = lat_df.set_axis(["latent_" + s for s in lat_df.columns.astype("str").tolist()], axis="columns")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    lat_df.to_csv(os.path.join(out_dir, "latent.csv"))


# PyTorch-BigGraph's eval uses a SPAWN subprocess that re-imports this module. Without this
# __main__ guard, the re-import re-runs the whole pipeline -> "_check_not_importing_main"
# RuntimeError (subprocess exited with status 1). The guard makes the re-import a no-op in the
# worker. (set_start_method('fork') did NOT help -- torchbiggraph hardcodes its own spawn context.)
if __name__ == "__main__":
    main()
