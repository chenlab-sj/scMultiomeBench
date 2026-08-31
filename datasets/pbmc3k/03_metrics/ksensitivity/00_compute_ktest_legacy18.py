#!/usr/bin/env python
# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# coding: utf-8

# In[1]:


import sys
absolute_path = '/path/to/multiomeBench/common'  
if absolute_path not in sys.path:
    sys.path.append(absolute_path)


# In[2]:


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


# In[13]:


from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score,precision_score, recall_score
from functools import reduce
from collections import Counter
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import confusion_matrix


# In[15]:


from functools import partial


# In[35]:


def knn_ataclabel(k,lat_dfs, pipelines,labels_annot,bc_test_rna, bc_test_atac, res_dist):
    df_preds = []
    df_pred_prob =[]
    for i in range(0, len(lat_dfs)):
        lat_df = lat_dfs[i]
        pipeline = pipelines[i]
        bc_test_rna_valid = [index for index in bc_test_rna if index in lat_df.index]
        bc_test_atac_valid = [index for index in bc_test_atac if index in lat_df.index]
        lat_df_rna =lat_df.loc[bc_test_rna_valid]
        lat_df_atac = lat_df.loc[bc_test_atac_valid]
        labels_df_rna = labels_annot.loc[labels_annot.index.isin(lat_df_rna.index)][["cell_type"]]
        labels_df_atac = labels_annot.loc[labels_annot.index.isin(lat_df_atac.index)][["cell_type"]]
        classifier = KNeighborsClassifier(n_neighbors=k,metric="cosine",weights ='distance')
        classifier.fit(lat_df_rna, labels_df_rna)
        predicted_labels = classifier.predict(lat_df_atac)
        prediction_probabilities = classifier.predict_proba(lat_df_atac)
        hit_prob = np.max(prediction_probabilities, axis=1)
        hit_prob_mean = np.mean(hit_prob)

        df_pred = pd.DataFrame({pipeline: predicted_labels}, index=lat_df_atac.index)
        df_preds.append(df_pred)
        accuracy = accuracy_score(labels_df_atac, predicted_labels)
        #balance_accuracy = balanced_accuracy_score(labels_df_atac, predicted_labels)
        #precision = precision_score(labels_df_atac, predicted_labels, average='weighted')
        #recall = recall_score(labels_df_atac, predicted_labels,average='weighted')
        df_pred_prob.append([pipeline, accuracy, "overall"])
        conf_matrix = confusion_matrix(labels_df_atac, predicted_labels)
        class_wise_accuracy = np.diagonal(conf_matrix) / np.sum(conf_matrix, axis=1)
        type = np.unique(predicted_labels)
        for j in range(0, len(type)):
            df_pred_prob.append([pipeline, class_wise_accuracy[j], type[j]])
    merged_df = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), df_preds)
    ## save KNN results
    knnres_file = "knn_k" + str(k) + "_pred_label.csv"
    knneval_file = "knn_k" + str(k) + "_pred_accu.csv"
    merged_df.to_csv(os.path.join(res_dist, knnres_file))
    df_pred_prob = pd.DataFrame(df_pred_prob, columns=['pipeline',"accuracy",'type'])
    df_pred_prob.to_csv(os.path.join(res_dist, knneval_file))



# In[36]:


def topk_smallest(lst,k):
    return sorted(lst)[:k]


# In[37]:


def col_typedist(col, label_dic):
    wi_celltype_dist=[]
    bt_celltype_dist=[]
    for index, value in col.items():

        col_type = label_dic.get(col.name)
        index_type = label_dic.get(index)


        if col_type == index_type:
            wi_celltype_dist.append(value)
        else:
            bt_celltype_dist.append(value)

    return wi_celltype_dist, bt_celltype_dist

def combine_lists(row):
    combined_list = []
    for sublist in row:
        combined_list.extend(sublist)
    return combined_list


# In[45]:


def get_celltype_disttopk(pdist_df,label_dic,k):
    dist_out = pdist_df.apply(lambda column: col_typedist(column, label_dic))
    #dist_out_topk = dist_out.applymap(topk_smallest)
    dist_out_topk = dist_out.applymap(lambda lst: topk_smallest(lst, k))
    #results= dist_out.apply(lambda row: combine_lists(row), axis=1)
    results= dist_out_topk.apply(lambda row: combine_lists(row), axis=1)
    wiclust_dist= results.loc[0]
    btclust_dist= results.loc[1]
    return([wiclust_dist, btclust_dist])


# In[46]:


def distance_eval(k, res_dist, lat_dfs, pipelines,bc_test, bc_test_rna, bc_test_atac, label_dic, cell_test,celltype):
    n_rows = len(lat_dfs)
    n_cols = len(celltype) 
    path = res_dist
    option = "Cosine"
    filename = option + "disttop" + str(k) + ".png"
    filename2 = option + "disttop" + str(k) + "_stat.csv"
    colnames = ['Pipeline','ks-bt_wi-stat','ks-bt_wi-pval',"type"]
    rows=[]
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(1.5*(n_cols + 1), 1*(n_rows + 1)))

    plt.rc('font', size=9)
    
    for i in range(0, len(lat_dfs)):

        ## calculate distance
        lat_df = lat_dfs[i]
        pipeline = pipelines[i]
        print("work on " + pipeline)
        pdist_df=benchmark_fun.lat2pdist(lat_df, bc_test,option)
        pdist_dfa = pdist_df.loc[pdist_df.index.isin(bc_test_rna), ]
        pdist_dfa= pdist_dfa[[col for col in pdist_dfa.columns if col in bc_test_atac]]
        #parallel_dist = benchmark_fun.get_parallel_dist(pdist_df,cell_test)
        celltype_dist = benchmark_fun.get_celltype_disttopk(pdist_dfa,label_dic,k)
        wiclust_dist=celltype_dist[0]
        btclust_dist=celltype_dist[1]
        ks2=stats.kstest(btclust_dist,wiclust_dist,alternative = 'less')
        rows.append([pipeline, ks2.statistic,ks2.pvalue,"overall"])
        # ax = axs[i, 0]
        # ax.hist(np.array(wiclust_dist), bins = 50, density=True, alpha = 0.5, color = 'tab:orange', label= "within celltype across modalities")
        # ax.hist(np.array(btclust_dist), bins = 50, density=True, alpha = 0.5, color = 'tab:green',label= "between celltype across modalities")
      
        for j in range(0, len(celltype)):
            type = celltype[j]
            print(type)
            celltype_atac = [key for key in label_dic.keys() if label_dic[key] == type and key.endswith("_atac")]
            pdist_dfa_celltype = pdist_dfa[[col for col in pdist_dfa.columns if col in celltype_atac]]
            celltype_dist = benchmark_fun.get_celltype_disttopk(pdist_dfa_celltype,label_dic,k)
            wiclust_dist=celltype_dist[0]
            btclust_dist=celltype_dist[1]
            ks2=stats.kstest(btclust_dist,wiclust_dist,alternative = 'less')
            rows.append([pipeline, ks2.statistic,ks2.pvalue,type])

            ax = axs[i, j]
            if i == 0:
                ax.set_title(type, fontsize = 9)
            if j==0:
                ax.set_ylabel(pipeline, fontsize = 9)

            ax.hist(np.array(wiclust_dist), bins = 50, density=True, alpha = 0.5, color = 'tab:blue', label= "intra-celltype neighbor across modalities")
            ax.hist(np.array(btclust_dist), bins = 50, density=True, alpha = 0.5, color = 'tab:orange',label= "inter-celltype neighbor across modalities")

    
    dist_sum_df = pd.DataFrame(rows, columns=colnames)
    dist_sum_df.to_csv(os.path.join(path, filename2),index=False)
    fig.text(0.5, 0.02, option + " distance ", ha='center', fontsize =10)
    fig.text(0.02, 0.5, 'Density', va='center', rotation='vertical', fontsize = 10)    
    labels = ["intra-celltype neighbor across modalities","inter-celltype neighbor across modalities"]
    fig.legend(labels=labels, loc='upper center', bbox_to_anchor=(0.5, 0.99), fontsize =10, ncol =2)
    
    plt.tight_layout(rect=(0.03,0.02,1,0.96))
    plt.subplots_adjust(hspace=0.3)

    plt.savefig(os.path.join(path, filename),dpi =300)
    plt.close(fig)




# In[ ]:


def main():   
    res_dist='/path/to/multiomeBench/common/pbmc3k' 
    out_dist= '/path/to/multiomeBench/common/pbmc3k/knn_test/' 
    if not os.path.exists(out_dist):
        os.makedirs(out_dist)
    
    cell_annot = os.path.join(res_dist, 'label.csv')
    #cell_annot = '/path/to/data/RMS/SJRHB013758_X1_label_lca0328.csv'
    labels_annot= pd.read_csv(cell_annot,index_col=0)
    bc_test=labels_annot[labels_annot['modality']!="train multiomics"].index.values
    bc_train=labels_annot[labels_annot['modality']=="train multiomics"].index.values
    bc_sel= labels_annot.index.values
    labels_test = labels_annot.loc[bc_test]
    bc_test_rna =[bc for bc in bc_test if bc.endswith("_rna")]
    bc_test_atac =[bc for bc in bc_test if bc.endswith("_atac")]
    label_dic = dict(zip(labels_test.index, labels_test['cell_type']))
    celltype =labels_test['cell_type'].unique().tolist()

    ## load latent file
    seurat3_lat_file = '/path/to/tools/Seuratv3/pbmc3k/res_pbmc3k_testall/lat_df.csv'
    bindsc_lat_file = '/path/to/tools/bindsc/pbmc3k/res_pbmc3k/coembed_coor.csv'
    concos_lat_file = '/path/to/tools/Concos/pbmc3k/res_pbmc3k/coembed_coor.csv'
    liger_lat_file = '/path/to/tools/liger/pbmc3k/res_pbmc3k/coembed_coor.csv'
    scBridge_lat_file =  '/path/to/tools/scBridge/pbmc3k/latent.csv'
    scdart_lat_file = '/path/to/tools/scDART/scDART/pbmc3k/res_pbmc3k/latent.csv'
    multimap_lat_file = '/path/to/tools/multimap/pbmc3k/res_pbmc3k/coembed_coor.csv'
    MinNet_lat_file = '/path/to/tools/MinNet/pbmc3k/coembed_coor.csv'
    scglue_lat_file = "/path/to/tools/scglue/pbmc3k/res_pbmc3k/latent.csv"
    scglue_withpair_lat_file='/path/to/tools/scglue/scglue_withpair/pbmc3k/scglue_latent.csv'
    scJoint_lat_file = "/path/to/tools/scJoint/pbmc3k/output/scJoint_latent.csv"
    cobolt_lat_file = "/path/to/tools/cobolt/pbmc3k/res_pbmc3k/cobolt_latent.csv"
    scVI_lat_file = "/path/to/tools/scvi/pbmc3k/res_pbmc3k/scvi_latent.csv"
    seurat4_lat_file = "/path/to/tools/Seuratv4/pbmc3k/res_pbmc3k_seurat4/seurat4_latent.csv"
    simba_lat_file = '/path/to/tools/simba/pbmc3k/latent.csv'
    scMoMaT_lat_file = '/path/to/tools/scMoMaT/pbmc3k/latent.csv'
    mmdma_lat_file = '/path/to/tools/lsmmdma/pbmc3k/latent.csv'
    unioncom_lat_file = '/path/to/tools/unioncom/pbmc3k/latent.csv'
    portal_lat_file = '/path/to/tools/portal/Portal/pbmc3k/lat_df.csv'

    seurat3_lat_df = pd.read_csv(seurat3_lat_file, index_col = 0)
    bindsc_lat_df = pd.read_csv(bindsc_lat_file, index_col = 0)
    concos_lat_df = pd.read_csv(concos_lat_file, index_col = 0)
    liger_lat_df = pd.read_csv(liger_lat_file, index_col = 0)
    scBridge_lat_df = pd.read_csv(scBridge_lat_file, index_col = 0)
    scdart_lat_df = pd.read_csv(scdart_lat_file, index_col = 0)
    multimap_lat_df = pd.read_csv(multimap_lat_file, index_col = 0)
    MinNet_lat_df = pd.read_csv(MinNet_lat_file, index_col = 0)
    scglue_lat_df = pd.read_csv(scglue_lat_file, index_col = 0)
    simba_lat_df=pd.read_csv(simba_lat_file,index_col=0)
    scMoMaT_lat_df=pd.read_csv(scMoMaT_lat_file,index_col=0)
    mmdma_lat_df=pd.read_csv(mmdma_lat_file,index_col=0)
    unioncom_lat_df=pd.read_csv(unioncom_lat_file,index_col=0)
    portal_lat_df=pd.read_csv(portal_lat_file,index_col=0)
    scglue_withpair_lat_df = pd.read_csv(scglue_withpair_lat_file, index_col = 0)
    scJoint_lat_df = pd.read_csv(scJoint_lat_file, index_col = 0)
    cobolt_lat_df = pd.read_csv(cobolt_lat_file, index_col = 0)
    scVI_lat_df = pd.read_csv(scVI_lat_file, index_col = 0)
    seurat4_lat_df = pd.read_csv(seurat4_lat_file, index_col = 0)

    lat_dfs = [bindsc_lat_df,concos_lat_df, liger_lat_df,multimap_lat_df,seurat3_lat_df,
            scdart_lat_df, scMoMaT_lat_df, simba_lat_df, 
            scglue_lat_df, #mmdma_lat_df, 
            unioncom_lat_df, portal_lat_df,scJoint_lat_df, scBridge_lat_df,
            cobolt_lat_df,scglue_withpair_lat_df,
            scVI_lat_df,seurat4_lat_df,MinNet_lat_df]
    pipelines = ["BindSC", "Conos", "LIGER","MultiMAP","Seurat(CCA)","scDART", "scMoMaT","simba",
                "scglue",#"MMD-MA",
                "Unioncom","Portal",
                "scJoint","scBridge",
                "Cobolt","scglue(multiome)",
                "scVI","Seurat(WNN)","MinNet"]
    
    concos_rna = concos_lat_df.index[concos_lat_df.index.str.endswith('_rna')]
    concos_rna = concos_rna.str.replace('_rna', '')
    concos_atac = concos_lat_df.index[concos_lat_df.index.str.endswith('_atac')]
    concos_atac = concos_atac.str.replace('_atac', '')
    concos_common = concos_rna.intersection(concos_atac)
    concos_common = list(concos_common)

    cell_test= [s.split('_')[0] for s in bc_test]
    cell_test=list(set(cell_test))
    cell_test2 = list(set(concos_common) & set(cell_test))
    print(cell_test2)
    
    labels_df = labels_annot
    cell_annot="cell_type"
    data_annot="modality"
    random = 420

    for k in [5,10,20,40,80]:
        print("k=" + str(k))
        knn_ataclabel(k,lat_dfs, pipelines,labels_annot,bc_test_rna, bc_test_atac, out_dist)
        distance_eval(k, out_dist, lat_dfs, pipelines,bc_test, bc_test_rna, bc_test_atac,label_dic,cell_test2,celltype )


# In[ ]:


if __name__ == "__main__":
    main()

