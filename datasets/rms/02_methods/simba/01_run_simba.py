# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# TorchBigGraph's PBG eval (si.tl.pbg_train) spawns worker processes via the 'spawn' start method
# (it hard-codes the spawn context, so forcing the global start method to 'fork' does NOT help).
# Each spawned worker re-imports this module; without an `if __name__ == "__main__"` guard it would
# re-run the whole script (-> recursive pbg_train -> "freeze_support()" RuntimeError, subprocess status 1).
# So: keep imports + seeding at module top-level (cheap, safe to re-run on import) and put ALL the
# executable work under the __main__ guard. (fork attempt kept below; harmless if ignored.)
import multiprocessing as _mp
try:
    _mp.set_start_method("fork", force=True)
except (RuntimeError, ValueError):
    pass
import os
import random
import numpy as np
import simba as si
import pandas as pd
import scanpy as sc

# ---- reproducibility replicate controls ----
SEED = int(os.environ.get("SEED", "420"))
REP = os.environ.get("REP", "")
# Seed every RNG SIMBA/TorchBigGraph draw from. NOTE: SIMBA's PBG training
# (si.tl.pbg_train) exposes no seed in its pbg_params, and TorchBigGraph's
# ConfigSchema has no seed field either, so PBG's own training RNG cannot be
# pinned via the API; these global seeds control all upstream/init randomness.
random.seed(SEED)
np.random.seed(SEED)
try:
    import torch
    torch.manual_seed(SEED)
except Exception:
    pass


if __name__ == "__main__":
    # Mast607A TEST multiome = the hg19 cellranger output (8291 cells = matches the original 1st-run
    # latent.csv; the hg38 re-process had a DIFFERENT filtered set of 8215). hg19 peaks are 1:start-end
    # (NO chr prefix) -> add 'chr' below and use genome='hg19'.
    multiome_h5 = "cd RMS/Mast607/data/Mast607A_TB19_22652/Mast607A_TB19_22652/outs/filtered_feature_bc_matrix.h5"
    out_dir = os.path.join('/path/to/multiomeBench/RMS/Mast607/script/simba/', REP)

    # Route SIMBA's working directory (result_simba/, graph0, pbg model) under REP so
    # replicates use SEPARATE workdirs and never clobber each other.
    workdir = os.path.join('./result_simba', REP) if REP else './result_simba'
    si.settings.set_workdir(workdir)

    ###########################################
    #start = timeit.default_timer()
    _x = sc.read_10x_h5(multiome_h5, gex_only=False); _x.var_names_make_unique()
    test_rna  = _x[:, _x.var["feature_types"] == "Gene Expression"].copy()
    test_atac = _x[:, _x.var["feature_types"] == "Peaks"].copy()

    test_rna.var_names_make_unique()
    #atac_peaks.var_names_make_unique()
    test_atac.var_names_make_unique()

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
    si.pp.select_variable_genes(test_rna, n_top_genes=int(os.environ.get("N_HVG", "4000")))   # N_HVG env (default 4000)

    # Drop RNA cells with zero expressed HVG before building the graph. The raw Mast607A multiome
    # (filtered_feature_bc_matrix.h5) keeps sparse cells that gen_graph rejects with
    # "Some nodes contain zero expressed highly_variable features". Must precede discretize +
    # infer_edges so test_rna stays cell-consistent with the C-C edges (adata_CrnaCatac).
    _hvg = test_rna.var['highly_variable'].values
    _hvgsum = np.asarray(test_rna[:, _hvg].X.sum(axis=1)).ravel()
    _keep = _hvgsum > 0
    print(f"[simba] dropping {int((~_keep).sum())} / {test_rna.n_obs} RNA cells with zero expressed HVG", flush=True)
    test_rna = test_rna[_keep].copy()

    si.tl.discretize(test_rna,n_bins=5) #discretize RNA expression

    test_atac.var[['chr', 'start', 'end']] = test_atac.var['gene_ids'].str.extract(r'([^:]+):(\d+)-(\d+)')
    # hg19/outs peaks are 1:start-end (NO chr prefix); simba's gene_scores genome ref is chr-prefixed (the
    # hg38 run matched chr1 peaks), so add 'chr' if missing, then score against hg19.
    if not test_atac.var['chr'].astype(str).str.startswith('chr').any():
        test_atac.var['chr'] = 'chr' + test_atac.var['chr'].astype(str)
    test_rna_atac = si.tl.gene_scores(test_atac,genome='hg19',use_gene_weigt=True, use_top_pcs=True)

    si.pp.filter_genes(test_rna_atac,min_n_cells=3)
    si.pp.cal_qc_rna(test_rna_atac)
    si.pp.normalize(test_rna_atac,method='lib_size')
    si.pp.log_transform(test_rna_atac)

    # infer_edges builds RNA<->ATAC edges from the genes SHARED between test_rna's HVG and test_rna_atac
    # (the ATAC gene-activity). Some cells are empty across that smaller shared set -> infer_edges raises
    # "Some nodes contain zero expressed highly_variable features". Drop those cells from BOTH entities, and
    # keep test_atac consistent with test_rna_atac so the downstream gen_graph C-P/C-C edges still line up.
    _hvg_names = test_rna.var_names[test_rna.var['highly_variable'].values]
    _shared = _hvg_names.intersection(test_rna_atac.var_names)
    _r = np.asarray(test_rna[:, _shared].X.sum(axis=1)).ravel()
    _a = np.asarray(test_rna_atac[:, _shared].X.sum(axis=1)).ravel()
    print(f"[simba] shared HVG/gene-activity features={len(_shared)}; "
          f"RNA 0-shared={int((_r==0).sum())}/{test_rna.n_obs}, "
          f"ATAC 0-shared={int((_a==0).sum())}/{test_rna_atac.n_obs}", flush=True)
    test_rna = test_rna[_r > 0].copy()
    _keep_atac = test_rna_atac.obs_names[_a > 0]
    test_rna_atac = test_rna_atac[_a > 0].copy()
    test_atac = test_atac[test_atac.obs_names.isin(_keep_atac)].copy()

    adata_CrnaCatac = si.tl.infer_edges(test_rna, test_rna_atac, n_components=15, k=15)

    # edges trimmed at TRIM_CUTOFF (env; default 0.5 = manuscript). LOWER it (e.g. TRIM_CUTOFF=0.1) to keep
    # weak cross-modal edges so more ATAC cells stay connected -> survive PBG (the lever for the 1796 collapse).
    _cut = float(os.environ.get("TRIM_CUTOFF", "0.5"))
    si.tl.trim_edges(adata_CrnaCatac, cutoff=_cut)

    # DIAGNOSTIC EXIT (fast, ~5 min): `SIMBA_DIAG=1 [TRIM_CUTOFF=0.1] python 01_run_simba.py` prints the number
    # that PREDICTS the output -- how many cells per modality still have a cross-modal edge after trimming
    # (cells with NO surviving edge are exactly what PBG drops). The dim ~=1796 at cutoff 0.5 is ATAC; watch
    # it climb toward ~8000 as you lower the cutoff. Stops BEFORE gen_graph/PBG.
    if os.environ.get("SIMBA_DIAG"):
        import sys
        _X = adata_CrnaCatac.X
        _obs_conn = int((np.asarray((_X != 0).sum(axis=1)).ravel() > 0).sum())
        _var_conn = int((np.asarray((_X != 0).sum(axis=0)).ravel() > 0).sum())
        print(f"[simba][DIAG] TRIM_CUTOFF={_cut} | adata_CrnaCatac (obs,var)={adata_CrnaCatac.shape} | "
              f"obs-with-edge={_obs_conn}  var-with-edge={_var_conn}  | the dim ~1796 at cutoff 0.5 is ATAC "
              f"(should climb toward ~8000 as cutoff drops). Exiting before gen_graph/PBG.", flush=True)
        sys.exit(0)

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
    #si.load_graph_stats(path='./result_simba/pbg/graph0')
    # load in model info for ('graph0')
    #si.load_pbg_config(path='./result_simba/pbg/graph0/model/')

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
