# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import sys
absolute_path = '/path/to/multiomeBench/common'  
if absolute_path not in sys.path:
    sys.path.append(absolute_path)

import numpy as np
import pandas as pd
import scanpy as sc
import igraph
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances
from sklearn.metrics import silhouette_score
from sklearn.metrics import cohen_kappa_score
import os
import benchmark_fun
from scipy import stats
import re

from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score,precision_score, recall_score
from functools import reduce
from collections import Counter
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix
import seaborn as sns
import scib


## wrap up to function
## wrap up to function
def benchmark_metrics(lat_dfs, pipelines, labels_df, celltype, bc_test, bc_test_rna, bc_test_atac, celltest1, celltest2, celltest3,
                     label_dic,  path):
    filename1 = 'sum_metrics.csv'

    colnames = ['method',
                'asw', 'asw_random', 
                'omics_asw',  'sample_asw',
                'ami',
                'ari', 
                 'ks.statistic']
    filename2 = 'celltype_metrics.csv'
    celltype_dflist = []
    celltype_colnames = ['method','celltype',
                'ks_celltype','ks_sample',
               'ks_inter_celltype', 'celltype_omics_ASW']

    ## other params
    nclust = len(labels_df["cell_type"].unique())
    cell_annot="cell_type"
    data_annot="modality"
    random = 420
    option = "Cosine"
    bc_test1_rna =[bc for bc in bc_test if re.search(r"_rna1", bc)]
    bc_test1_atac =[bc for bc in bc_test if re.search(r"_atac1", bc)]
    bc_test2_rna =[bc for bc in bc_test if re.search(r"_rna2", bc)]
    bc_test2_atac =[bc for bc in bc_test if re.search(r"_atac2", bc)]
    bc_test3_rna =[bc for bc in bc_test if re.search(r"_rna3", bc)]
    bc_test3_atac =[bc for bc in bc_test if re.search(r"_atac3", bc)]
    rows = []

    for i in range(0, len(lat_dfs)):
        lat_df = lat_dfs[i]
        method = pipelines[i]
        print(method)
        bc_test_sel =  list(set(lat_df.index).intersection(bc_test))
        lat_df_test = lat_df.loc[lat_df.index.isin(bc_test_sel)]
        labels_test=labels_df.loc[bc_test_sel]
        adata = sc.AnnData(lat_df_test)
        adata_pp = adata.copy()
        adata.obsm['X_emb'] = lat_df_test

        ## scib: asw/ilsi
        adata.obs = pd.merge(adata.obs,labels_test, how='left', left_index=True, right_index=True)
    
        adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'])
        adata.obs['random_celltype'] = pd.Categorical(adata.obs['random_celltype'])
        asw = scib.me.silhouette(adata, label_key="cell_type", embed="X_emb") ## celltype asw
        asw_random = scib.me.silhouette(adata, label_key="random_celltype", embed="X_emb")
    
        adata.obs['modality'] = pd.Categorical(adata.obs['modality'])
        #adata.obs['random_modality'] = pd.Categorical(adata.obs['random_modality'])
        omics_asw = scib.me.silhouette_batch(adata, batch_key="modality",label_key="cell_type", embed="X_emb",return_all=True) ## celltype asw
        #omics_asw_random = scib.me.silhouette_batch(adata, batch_key="random_modality",label_key="cell_type", embed="X_emb")
        omics_asw_ave = omics_asw[0]
        omics_asw_celltype = omics_asw[1]
        adata.obs['batch'] = pd.Categorical(adata.obs['batch'])
        sample_asw = scib.me.silhouette_batch(adata, batch_key="batch",label_key="cell_type", embed="X_emb")

        #ilisi = scib.me.ilisi_graph(adata, batch_key="modality", type_="embed", use_rep="X_emb")
        #ilisi_random = scib.me.ilisi_graph(adata, batch_key="random_modality", type_="embed", use_rep="X_emb")
    
        ## nmi/ari
        louvain_res = benchmark_fun.find_louvain_res(adata, nclust)
        adata_pp = adata.copy()
        sc.pp.neighbors(adata_pp)
        sc.tl.louvain(adata_pp, resolution = louvain_res, random_state = random)
        labels = adata_pp.obs['louvain']
        labels_test = labels_df.loc[bc_test_sel]
    
        modality_consis=benchmark_fun.cluster_cosistent3(adata_pp.obs,celltest1, celltest2, celltest3, 'louvain')
        ami = modality_consis[0]
        ari = modality_consis[1]
        #nmi= modality_consis[2]

    
        #random_modality_consis = benchmark_fun.cluster_cosistent(adata.obs,cell_test,'random_celltype')
        #random_ami = random_modality_consis[0]
        #random_ari = random_modality_consis[1]
        #random_nmi = random_modality_consis[2]
        ##  distance
        pdist_df=benchmark_fun.lat2pdist(lat_df, bc_test,option)
        samecell_dist = benchmark_fun.get_parallel_dist3(pdist_df,celltest1, celltest2, celltest3)
        pdist_dfintera = pdist_df.loc[pdist_df.index.isin(bc_test_rna), ]
        pdist_dfintera= pdist_dfintera[[col for col in pdist_dfintera.columns if col in bc_test_atac]]
        pdist_df_rna = pdist_df.loc[pdist_df.index.isin(bc_test_rna), ]
        pdist_df_rna= pdist_df_rna[[col for col in pdist_df_rna.columns if col in bc_test_rna]]
        pdist_df_atac = pdist_df.loc[pdist_df.index.isin(bc_test_atac), ]
        pdist_df_atac= pdist_df_atac[[col for col in pdist_df_atac.columns if col in bc_test_atac]]

        print(len(samecell_dist))
        ## same omics distance
        sameomics_dist =benchmark_fun.get_samemodality_dist(pdist_df_rna, pdist_df_atac)
        #print("same omics dist")
        print(len(sameomics_dist))
        #ks=stats.kstest(sameomics_dist,samecell_dist,alternative = 'less') ## no need to random?
        #rows.append([method, asw, asw_random, omics_asw_ave, sample_asw, ami, ari,  ks.statistic])

        # same sample distance
        pdist_df_rna1 = pdist_df.loc[pdist_df.index.isin(bc_test1_rna), ]
        pdist_df_atac1 = pdist_df.loc[pdist_df.index.isin(bc_test1_atac), ]
        pdist_df_rna2 = pdist_df.loc[pdist_df.index.isin(bc_test2_rna), ]
        pdist_df_atac2 = pdist_df.loc[pdist_df.index.isin(bc_test2_atac), ]
        pdist_df_rna3 = pdist_df.loc[pdist_df.index.isin(bc_test3_rna), ]
        pdist_df_atac3 = pdist_df.loc[pdist_df.index.isin(bc_test3_atac), ]
        # pdist_df_rna1_self= pdist_df_rna1[[col for col in pdist_df_rna1.columns if col in bc_test1_rna]]
        # pdist_df_rna2_self= pdist_df_rna2[[col for col in pdist_df_rna2.columns if col in bc_test2_rna]]
        # pdist_df_atac1_self= pdist_df_atac1[[col for col in pdist_df_atac1.columns if col in bc_test1_atac]]
        # pdist_df_atac2_self= pdist_df_atac2[[col for col in pdist_df_atac2.columns if col in bc_test2_atac]]
        # samesample_dist = dist_df2val(pdist_df_rna1_self)+ dist_df2val(pdist_df_rna2_self) + dist_df2val(pdist_df_atac1_self)+ dist_df2val(pdist_df_atac2_self)
        ## cell type distance
        celltype_rows = []
        intra_celltype_sameomics_mean = []
        for j in range(0, len(celltype)):
        #for j in range(0, 1):
            type=celltype[j]
            print(type)
            celltype_atac1 = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_atac1")]
            celltype_rna1 = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_rna1")]
            celltype_atac2 = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_atac2")]
            celltype_rna2 = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_rna2")]
            celltype_atac3 = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_atac3")]
            celltype_rna3 = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_rna3")]
            celltype_atac = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith(("_atac1", "_atac2" ,"_atac3"))]
            #celltype_rna = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith(("rna1", "_rna2" ,"_rna3"))]
            #celltype_rna = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_rna")]
            # type_bc = [key for key in label_dic.keys() if label_dic[key] == type and key in lat_df.index ]
            # lat_df_celltype = lat_df.loc[type_bc]

            ## same sample 

            pdist_dfrna1_celltype = pdist_df_rna1[[col for col in pdist_df_rna1.columns if col in celltype_rna1]]
            pdist_dfatac1_celltype = pdist_df_atac1[[col for col in pdist_df_atac1.columns if col in celltype_atac1]]
            pdist_dfrna2_celltype = pdist_df_rna2[[col for col in pdist_df_rna2.columns if col in celltype_rna2]]
            pdist_dfatac2_celltype = pdist_df_atac2[[col for col in pdist_df_atac2.columns if col in celltype_atac2]]
            pdist_dfrna3_celltype = pdist_df_rna3[[col for col in pdist_df_rna3.columns if col in celltype_rna3]]
            pdist_dfatac3_celltype = pdist_df_atac3[[col for col in pdist_df_atac3.columns if col in celltype_atac3]]
            rna1_celltype_dist = benchmark_fun.get_celltype_dist(pdist_dfrna1_celltype,label_dic)
            atac1_celltype_dist = benchmark_fun.get_celltype_dist(pdist_dfatac1_celltype,label_dic)
            rna2_celltype_dist = benchmark_fun.get_celltype_dist(pdist_dfrna2_celltype,label_dic)
            atac2_celltype_dist = benchmark_fun.get_celltype_dist(pdist_dfatac2_celltype,label_dic)
            rna3_celltype_dist = benchmark_fun.get_celltype_dist(pdist_dfrna3_celltype,label_dic)
            atac3_celltype_dist = benchmark_fun.get_celltype_dist(pdist_dfatac3_celltype,label_dic)  
            wiclust_rna1dist=rna1_celltype_dist[0]
            btclust_rna1dist=rna1_celltype_dist[1]
            wiclust_atac1dist=atac1_celltype_dist[0]
            btclust_atac1dist=atac1_celltype_dist[1] 
            wiclust_rna2dist=rna2_celltype_dist[0]
            btclust_rna2dist=rna2_celltype_dist[1]
            wiclust_atac2dist=atac2_celltype_dist[0]
            btclust_atac2dist=atac2_celltype_dist[1]
            wiclust_rna3dist=rna3_celltype_dist[0]
            btclust_rna3dist=rna3_celltype_dist[1]
            wiclust_atac3dist=atac3_celltype_dist[0]
            btclust_atac3dist=atac3_celltype_dist[1]
            wiclust_sampledist = wiclust_atac1dist + wiclust_atac2dist + wiclust_atac3dist  + wiclust_rna1dist + wiclust_rna2dist + wiclust_rna3dist
            btclust_sampledist = btclust_atac1dist + btclust_atac2dist + btclust_atac3dist  + btclust_rna1dist + btclust_rna2dist + btclust_rna3dist
            ks2 = stats.kstest(btclust_sampledist,wiclust_sampledist,alternative = 'less') ## bio conservation





                        ## same omics, same sample, intra_celltype


            intracelltype_pdist_dfrna1 = pdist_dfrna1_celltype.loc[pdist_dfrna1_celltype.index.isin(celltype_rna1)]
            intra_celltype_rna1_mean = intracelltype_pdist_dfrna1.mean(axis =1)
            intra_celltype_rna1_mean = intra_celltype_rna1_mean.tolist()           
            intracelltype_pdist_dfatac1 = pdist_dfatac1_celltype.loc[pdist_dfatac1_celltype.index.isin(celltype_atac1)]
            intra_celltype_atac1_mean = intracelltype_pdist_dfatac1.mean(axis =1)
            intra_celltype_atac1_mean = intra_celltype_atac1_mean.tolist()

            intra_celltype_sameomics_mean.extend(intra_celltype_rna1_mean)
            intra_celltype_sameomics_mean.extend(intra_celltype_atac1_mean)

            intracelltype_pdist_dfrna2 = pdist_dfrna2_celltype.loc[pdist_dfrna2_celltype.index.isin(celltype_rna2)]
            intra_celltype_rna2_mean = intracelltype_pdist_dfrna2.mean(axis =1)
            intra_celltype_rna2_mean = intra_celltype_rna2_mean.tolist()           
            intracelltype_pdist_dfatac2 = pdist_dfatac2_celltype.loc[pdist_dfatac2_celltype.index.isin(celltype_atac2)]
            intra_celltype_atac2_mean = intracelltype_pdist_dfatac2.mean(axis =1)
            intra_celltype_atac2_mean = intra_celltype_atac2_mean.tolist()

            intra_celltype_sameomics_mean.extend(intra_celltype_rna2_mean)
            intra_celltype_sameomics_mean.extend(intra_celltype_atac2_mean)


            intracelltype_pdist_dfrna3 = pdist_dfrna3_celltype.loc[pdist_dfrna3_celltype.index.isin(celltype_rna3)]
            intra_celltype_rna3_mean = intracelltype_pdist_dfrna3.mean(axis =1)
            intra_celltype_rna3_mean = intra_celltype_rna3_mean.tolist()           
            intracelltype_pdist_dfatac3 = pdist_dfatac3_celltype.loc[pdist_dfatac3_celltype.index.isin(celltype_atac3)]
            intra_celltype_atac3_mean = intracelltype_pdist_dfatac3.mean(axis =1)
            intra_celltype_atac3_mean = intra_celltype_atac3_mean.tolist()

            intra_celltype_sameomics_mean.extend(intra_celltype_rna3_mean)
            intra_celltype_sameomics_mean.extend(intra_celltype_atac3_mean)

            print("same omics inter sample")
            pdist_dfrna_sample_celltype2 = pdist_df_rna2[[col for col in pdist_df_rna2.columns if col in celltype_rna3]]
            pdist_dfatac_sample_celltype2 = pdist_df_atac2[[col for col in pdist_df_atac2.columns if col in celltype_atac3]]

            rna1_intersample = celltype_rna2 + celltype_rna3
            pdist_dfrna_sample_celltype1 = pdist_df_rna1[[col for col in pdist_df_rna1.columns if col in rna1_intersample]]
            atac1_intersample = celltype_atac2 + celltype_atac3
            pdist_dfatac_sample_celltype1 = pdist_df_atac1[[col for col in pdist_df_atac1.columns if col in atac1_intersample]]
            
            rna_btsample_celltype_dist1 = benchmark_fun.get_celltype_dist(pdist_dfrna_sample_celltype1,label_dic)
            rna_btsample_celltype_dist2 = benchmark_fun.get_celltype_dist(pdist_dfrna_sample_celltype2,label_dic)
            atac_btsample_celltype_dist1 = benchmark_fun.get_celltype_dist(pdist_dfatac_sample_celltype1,label_dic)
            atac_btsample_celltype_dist2 = benchmark_fun.get_celltype_dist(pdist_dfatac_sample_celltype2,label_dic)
            wiclust_rnabtsampledist=rna_btsample_celltype_dist1[0] + rna_btsample_celltype_dist2[0]
            #btclust_rnasampledist=rna_celltype_dist[1]
            wiclust_atacbtsampledist=atac_btsample_celltype_dist1[0] + atac_btsample_celltype_dist2[0]
            #btclust_atacsampledist=atac_celltype_dist[1]
            wiclust_btsampledist = wiclust_rnabtsampledist + wiclust_atacbtsampledist
            ks3=stats.kstest(btclust_sampledist,wiclust_btsampledist,alternative = 'less')


            ## scRNA, scATAC alignment
            
            pdist_dfintera_celltype = pdist_dfintera[[col for col in pdist_dfintera.columns if col in celltype_atac]]
            inter_celltype_dist = benchmark_fun.get_celltype_dist(pdist_dfintera_celltype,label_dic)
            wiclust_interdist=inter_celltype_dist[0]
            btclust_interdist=inter_celltype_dist[1]
            ks4=stats.kstest(btclust_interdist,wiclust_interdist,alternative = 'less')


            ## inter sample

        ### get random
            # print(random)
            # random_celltype_atac = [key for key in label_dic.keys() if random_label_dic[key] == type and key.endswith("_atac")]
            # random_celltype_rna = [key for key in label_dic.keys() if random_label_dic[key] == type and key.endswith("_rna")]
            # random_type_bc = [key for key in label_dic.keys() if random_label_dic[key] == type and key in lat_df.index ]
            # random_lat_df_celltype = lat_df.loc[random_type_bc]
            # random_pdist_dfintera_celltype = pdist_dfintera[[col for col in pdist_dfintera.columns if col in random_celltype_atac]]
            # random_pdist_dfrna_celltype = pdist_df_rna[[col for col in pdist_df_rna.columns if col in random_celltype_rna]]
            # random_pdist_dfatac_celltype = pdist_df_atac[[col for col in pdist_df_atac.columns if col in random_celltype_atac]]
            # random_inter_celltype_dist = benchmark_fun.get_celltype_dist(random_pdist_dfintera_celltype,random_label_dic)
            # random_wiclust_interdist=random_inter_celltype_dist[0]
            # random_btclust_interdist=random_inter_celltype_dist[1]
            # random_ks1=stats.kstest(random_btclust_interdist,random_wiclust_interdist,alternative = 'less')
            
            # random_rna_celltype_dist = benchmark_fun.get_celltype_dist(random_pdist_dfrna_celltype,random_label_dic)
            # random_atac_celltype_dist = benchmark_fun.get_celltype_dist(random_pdist_dfatac_celltype,random_label_dic)
            # random_wiclust_rnadist=random_rna_celltype_dist[0]
            # random_btclust_rnadist=random_rna_celltype_dist[1]
            # random_wiclust_atacdist=random_atac_celltype_dist[0]
            # random_btclust_atacdist=random_atac_celltype_dist[1]
            # random_wiclust_modadist = random_wiclust_atacdist + random_wiclust_rnadist
            # random_btclust_modadist = random_btclust_atacdist + random_btclust_rnadist
            # random_ks2=stats.kstest(random_btclust_modadist,random_wiclust_modadist,alternative = 'less')
            celltype_rows.append([method, type, ks2.statistic,ks3.statistic,ks4.statistic   ])
        celltype_stat = pd.DataFrame(celltype_rows)
        celltype_df = pd.merge(celltype_stat, omics_asw_celltype, left_on =1, right_index = True)
        celltype_dflist.append(celltype_df)

        ks=stats.kstest(intra_celltype_sameomics_mean,samecell_dist,alternative = 'less')
        rows.append([method, asw, asw_random, omics_asw_ave, sample_asw, ami, ari,  ks.statistic])
    sum_matrix_df = pd.DataFrame(rows, columns=colnames)
    sum_matrix_df.to_csv(os.path.join(path, filename1),index=False)

    celltype_matrix_df =  pd.concat(celltype_dflist, axis=0, ignore_index=True)
    celltype_matrix_df.columns = celltype_colnames
    celltype_matrix_df.to_csv(os.path.join(path, filename2),index=False)


def main():
    res_dist='/path/to/multiomeBench/common/BMMC_d1/' 
    out_dist= '/path/to/multiomeBench/common/BMMC_d1/benchmark_matrix/' 
    if not os.path.exists(out_dist):
        os.makedirs(out_dist)
    
    cell_annot = os.path.join(res_dist, 'label.csv')
    labels_annot= pd.read_csv(cell_annot,index_col=0)

    bc_test=labels_annot[labels_annot['set']!='train multiomics'].index.values
    bc_train=labels_annot[labels_annot['set']=='train multiomics'].index.values
    bc_sel= labels_annot.index.values
    labels_test = labels_annot.loc[bc_test]
    bc_test_rna =[bc for bc in bc_test if re.search(r"_rna.*$", bc)]
    bc_test_atac =[bc for bc in bc_test if re.search(r"_atac.*$", bc)]
    label_dic = dict(zip(labels_test.index, labels_test['cell_type']))
    random_label_dic = dict(zip(labels_test.index, labels_test['random_celltype']))
    celltype =labels_test['cell_type'].unique().tolist()

    seurat3_lat_file = '/path/to/data/BMMC_d1/seurat3/coembed_coor.csv'
    scBridge_lat_file = "/path/to/data/BMMC_d1/scBridge/latent.csv"
    simba_lat_file = '/path/to/data/BMMC_d1/simba/latent.csv'
    portal_lat_file = '/path/to/data/BMMC_d1/portal/lat_df.csv'
    scdart_lat_file = '/path/to/data/BMMC_d1/scdart/latent.csv'
    scglue_lat_file="/path/to/data/BMMC_d1/scglue/latent.csv"
    scglue2_lat_file="/path/to/data/BMMC_d1/scglue2/latent.csv"
    scglue_paired_lat_file="/path/to/data/BMMC_d1/scglue_paired/scglue_latent.csv"
    scglue2_paired_lat_file="/path/to/data/BMMC_d1/scglue_paired2/scglue_latent.csv"
    scVI_lat_file="/path/to/data/BMMC_d1/scVI/scVI_latent.csv"
    scVI2_lat_file="/path/to/data/BMMC_d1/scVI2/scVI_latent.csv"
    scjoint_lat_file = "/path/to/data/BMMC_d1/scjoint/scJoint_latent.csv"
    scjoint2_lat_file = "/path/to/data/BMMC_d1/scjoint2/scJoint_latent.csv"
    # MinNet_lat_file = '/path/to/data/BMMC_d1/MinNet/coembed_coor.csv'
    # MinNet2_lat_file = '/path/to/data/BMMC_d1/MinNet2/coembed_coor.csv'
    bindsc_lat_file = '/path/to/data/BMMC_d1/bindsc/coembed_coor.csv'
    bindsc2_lat_file = '/path/to/data/BMMC_d1/bindsc2/coembed_coor.csv'
    cobolt_lat_file = "/path/to/data/BMMC_d1/BMMC_d1_train/cobolt/cobolt_latent.csv"

    seurat3_lat_df=pd.read_csv(seurat3_lat_file,index_col=0)
    scBridge_lat_df=pd.read_csv(scBridge_lat_file,index_col=0)
    simba_lat_df=pd.read_csv(simba_lat_file,index_col=0)
    portal_lat_df=pd.read_csv(portal_lat_file,index_col=0)
    scdart_lat_df=pd.read_csv(scdart_lat_file,index_col=0)
    scglue_lat_df=pd.read_csv(scglue_lat_file,index_col=0)
    scglue2_lat_df=pd.read_csv(scglue2_lat_file,index_col=0)
    scglue_paired_lat_df=pd.read_csv(scglue_paired_lat_file,index_col=0)
    scglue2_paired_lat_df=pd.read_csv(scglue2_paired_lat_file,index_col=0)
    scVI_lat_df=pd.read_csv(scVI_lat_file,index_col=0)
    scVI2_lat_df=pd.read_csv(scVI2_lat_file,index_col=0)
    scjoint_lat_df=pd.read_csv(scjoint_lat_file,index_col=0)
    scjoint2_lat_df=pd.read_csv(scjoint2_lat_file,index_col=0)
    # MinNet_lat_df=pd.read_csv(MinNet_lat_file,index_col=0)
    # MinNet2_lat_df=pd.read_csv(MinNet2_lat_file,index_col=0)
    bindsc_lat_df = pd.read_csv(bindsc_lat_file,index_col=0)
    bindsc2_lat_df = pd.read_csv(bindsc2_lat_file,index_col=0)
    cobolt_lat_df = pd.read_csv(cobolt_lat_file,index_col=0)


    # lat_dfs = [seurat3_lat_df, scBridge_lat_df, simba_lat_df,
    #         portal_lat_df, scdart_lat_df,
    #         scglue_lat_df, scglue2_lat_df, scglue_paired_lat_df, scglue2_paired_lat_df,
    #         scVI_lat_df, scVI2_lat_df, scjoint_lat_df, scjoint2_lat_df,
    #         bindsc_lat_df, bindsc2_lat_df,cobolt_lat_df
    #             ]
    # pipelines = ['Seurat(CCA)','scBridge','simba',
    #             'Portal','scDART',
    #             'scglue','scglue(batch)','scglue(multiome)','scglue(multiome,batch)',
    #             'scVI','scVI(batch)','scJoint','scJoint(batch)',
    #             'BindSC',"BindSC(batch)",'Cobolt']
        
    
    lat_dfs = [seurat3_lat_df, scBridge_lat_df, simba_lat_df,
            portal_lat_df, scdart_lat_df
                ]
    pipelines = ['Seurat(CCA)','scBridge','simba',
                'Portal','scDART']
    # conos_rna = conos_lat_df.index[conos_lat_df.index.str.endswith('_rna')]
    # conos_rna = conos_rna.str.replace('_rna', '')
    # conos_atac = conos_lat_df.index[conos_lat_df.index.str.endswith('_atac')]
    # conos_atac = conos_atac.str.replace('_atac', '')
    # conos_common = conos_rna.intersection(conos_atac)
    # conos_common = list(conos_common)
        
    cell_test= [s.split('_')[0] + "_" + s.split('_')[1] for s in bc_test]
    cell_test=list(set(cell_test))
    # cell_test2 = list(set(conos_common) & set(cell_test))
    #print(cell_test2)

        ## get complete benchmark matrix
    celltest1 = []
    celltest2 = []
    celltest3 = []
    for item in bc_test_rna:
        match = re.match(r"^(.*)_rna(\d+)$", item)
        if match:
            prefix = match.group(1)  # Extracted prefix
            suffix = int(match.group(2))  # Extracted number at the end

            # Group based on the extracted suffix
            if suffix == 1:
                celltest1.append(prefix)
            elif suffix == 2:
                celltest2.append(prefix)
            elif suffix == 3:
                celltest3.append(prefix)

    labels_df = labels_annot
    
    benchmark_metrics(lat_dfs, pipelines, labels_df, celltype, bc_test, bc_test_rna, bc_test_atac, celltest1, celltest2, celltest3,
                     label_dic,  out_dist)
if __name__ == "__main__":
    main()
