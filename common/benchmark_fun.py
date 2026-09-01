## benchmark_fun.py
import numpy as np
import pandas as pd
import scanpy as sc
import igraph
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi,adjusted_mutual_info_score as ami
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances
from sklearn.metrics import silhouette_score
from sklearn.metrics import cohen_kappa_score
import os
from scipy.spatial.distance import cdist
import seaborn as sns
###################################################
## functions for cell type clustering evaluation
###################################################
## input: latent in anndata format, number of cluster
## out: resolution value
def find_louvain_res(adata,n_cluster): 
    random = 0
    adata_pp = adata.copy()
    sc.pp.neighbors(adata_pp)
    obtained_clusters = -1
    iteration = 0
    resolutions = [0, 100]
    
    while obtained_clusters != n_cluster and iteration < 200:
        res = sum(resolutions)/2    
        sc.tl.louvain(adata_pp, resolution = res, random_state = random)
        labels = adata_pp.obs['louvain']
        obtained_clusters = len(np.unique(labels))

        if obtained_clusters < n_cluster:
            resolutions[0] = res
        else:
            resolutions[1] = res 
        iteration = iteration + 1

    print("resoultion: ",res," for ",obtained_clusters," clusters after ", iteration, "iteration")
    return res

## input: cluster_res is a dataframe  with predicted cluster result
######### cell_test is a list of 
def cluster_cosistent(cluster_res,cell_test,label):
    ## check whether same cell in same cluster
    ##test bc
    #bc_test=cluster_res.index.values
                
    #count=0
    rna_clust_label=[]
    atac_clust_label=[]
    for cell in cell_test:
        rna_clust = cluster_res.loc[cell+"_rna"][label]
        atac_clust = cluster_res.loc[cell+"_atac"][label]
        rna_clust_label.append(rna_clust)
        atac_clust_label.append(atac_clust)
        #if (rna_clust == atac_clust):
         #   count +=1
            
    
    #modality_consis =count/len(cell_test)
    #kappa = cohen_kappa_score(rna_clust_label, atac_clust_label)
    AMI = ami(rna_clust_label, atac_clust_label)
    ARI = ari(rna_clust_label, atac_clust_label)
    NMI = nmi(rna_clust_label, atac_clust_label)
    
    return([AMI, ARI, NMI])


def cluster_cosistent2(cluster_res,common_bc1, common_bc2, label):

    rna1_clust_label=[]
    atac1_clust_label=[]

    rna2_clust_label=[]
    atac2_clust_label=[]
    for cell in common_bc1:
        rna_clust = cluster_res.loc[cell+"_rna1"][label]
        atac_clust = cluster_res.loc[cell+"_atac1"][label]
        rna1_clust_label.append(rna_clust)
        atac1_clust_label.append(atac_clust)
    for cell in common_bc2:
        rna_clust = cluster_res.loc[cell+"_rna2"][label]
        atac_clust = cluster_res.loc[cell+"_atac2"][label]
        rna2_clust_label.append(rna_clust)
        atac2_clust_label.append(atac_clust)
        
    
    rna_clust_label = rna1_clust_label + rna2_clust_label
    atac_clust_label = atac1_clust_label + atac2_clust_label
    AMI = ami(rna_clust_label, atac_clust_label)
    ARI = ari(rna_clust_label, atac_clust_label)
    NMI = nmi(rna_clust_label, atac_clust_label)
    
    return([AMI, ARI, NMI])


def cluster_cosistent3(cluster_res,common_bc1, common_bc2, common_bc3, label):

    rna1_clust_label=[]
    atac1_clust_label=[]

    rna2_clust_label=[]
    atac2_clust_label=[]

    rna3_clust_label=[]
    atac3_clust_label=[]

    for cell in common_bc1:
        rna_clust = cluster_res.loc[cell+"_rna1"][label]
        atac_clust = cluster_res.loc[cell+"_atac1"][label]
        rna1_clust_label.append(rna_clust)
        atac1_clust_label.append(atac_clust)
    for cell in common_bc2:
        rna_clust = cluster_res.loc[cell+"_rna2"][label]
        atac_clust = cluster_res.loc[cell+"_atac2"][label]
        rna2_clust_label.append(rna_clust)
        atac2_clust_label.append(atac_clust)

    for cell in common_bc3:
        rna_clust = cluster_res.loc[cell+"_rna3"][label]
        atac_clust = cluster_res.loc[cell+"_atac3"][label]
        rna3_clust_label.append(rna_clust)
        atac3_clust_label.append(atac_clust)
    
    rna_clust_label = rna1_clust_label + rna2_clust_label + rna3_clust_label
    atac_clust_label = atac1_clust_label + atac2_clust_label + atac3_clust_label
    AMI = ami(rna_clust_label, atac_clust_label)
    ARI = ari(rna_clust_label, atac_clust_label)
    NMI = nmi(rna_clust_label, atac_clust_label)
    
    return([AMI, ARI, NMI])
    
def get_clusterindex(cluster_res,lat_df):
    #random = 0
    #adata_pp = adata.copy()
    #sc.pp.neighbors(adata_pp)
    #sc.tl.louvain(adata_pp, resolution = res, random_state = random)
    labels = cluster_res['louvain']
    #annot_df_reorder = annot_df.reindex(labels.index)       
    annot = cluster_res['Cell type'].tolist()
    NMI = nmi(annot, labels)
    ARI = ari(annot, labels)
    
    lat_df_test = lat_df.loc[cluster_res.index]
    silhouette_sum= silhouette_score(lat_df_test, annot)
    ASW = silhouette_sum.mean()
    
    return([NMI,ARI,ASW])


def get_clusterindex2(cluster_res,lat_df):
    #random = 0
    #adata_pp = adata.copy()
    #sc.pp.neighbors(adata_pp)
    #sc.tl.louvain(adata_pp, resolution = res, random_state = random)
    labels = cluster_res['louvain']
    #annot_df_reorder = annot_df.reindex(labels.index)       
    annot = cluster_res['Cell type'].tolist()
    NMI = nmi(annot, labels)
    AMI = ami(annot, labels)
    lat_df_test = lat_df.loc[cluster_res.index]
    silhouette_sum= silhouette_score(lat_df_test, annot)
    ASW = silhouette_sum.mean()
    
    cluster_res_rna = cluster_res[cluster_res.index.str.endswith('_rna')]
    cluster_res_atac = cluster_res[cluster_res.index.str.endswith('_atac')]
    
    labels_rna = cluster_res_rna['louvain']
    annot_rna = cluster_res_rna['Cell type'].tolist()
    NMI_rna = nmi(annot_rna, labels_rna)
    ARI_rna = ari(annot_rna, labels_rna)
    lat_df_rna = lat_df.loc[cluster_res_rna.index]
    silhouette_sum_rna= silhouette_score(lat_df_rna, annot_rna)    
    ASW_rna = silhouette_sum_rna.mean()
    
    labels_atac = cluster_res_atac['louvain']
    annot_atac = cluster_res_atac['Cell type'].tolist()
    NMI_atac = nmi(annot_atac, labels_atac)
    ARI_atac = ari(annot_atac, labels_atac)
    lat_df_atac = lat_df.loc[cluster_res_atac.index]
    silhouette_sum_atac= silhouette_score(lat_df_atac, annot_atac)    
    ASW_atac = silhouette_sum_atac.mean()
    
    return([NMI,ARI,ASW, NMI_rna, ARI_rna, ASW_rna, NMI_atac, ARI_atac, ASW_atac])

###################################################
## functions for distance evaluation
###################################################

# ## function to get pair-wise distance from latent datframe of each bc
# ## option: cosine/euclidean
# def lat2pdist(lat_df, option):
#     if (option == "Cosine"):
#         pdist = cosine_distances(lat_df)
#         pdist_df = pd.DataFrame(pdist, index = lat_df.index, columns = lat_df.index)
#     else:
#         pdist = euclidean_distances(lat_df)
#         pdist_df = pd.DataFrame(pdist, index = lat_df.index, columns = lat_df.index)
        
#     return (pdist_df)

# ## same bc distance across modalities
# def get_parallel_dist(pdist_df,bc_test):
#     parallel_dist1 = [] ## rna -> atac
#     parallel_dist2 = [] ## atac -> rna
#     for bc in bc_test:
#         bc_rna= bc + "_rna"
#         bc_atac=bc + "_atac"
#         parallel_dist1.append(pdist_df.loc[bc_rna, bc_atac])
#         #parallel_dist2.append(pdist_df.loc[bc_atac, bc_rna])
#     #if np.array_equal(parallel_dist1, parallel_dist2):
#     #    print("dist of rna -> atac vs atac-> rna are equal")
#     #else:
#      #   print("dist of rna -> atac vs atac-> rna are not equal")
    
#     #parallel_dist=np.concatenate((parallel_dist1, parallel_dist2))
#     return(parallel_dist1)


# def get_celtype_dic(labels_df,bc_test,testBatch=True):   
#     ## filter test bc
#     test_labels_annot = labels_df.loc[bc_test]
#     celltype_grouped = test_labels_annot.groupby('cell_type')
#     celltype_dic = {}
#     for group_name, group_df in celltype_grouped:
#         if len(group_df)>1 :
#             indices = group_df.index.tolist()
#             if testBatch==True:
#                 indices = [s.split('_')[0] + "_" + s.split('_')[1] for s in indices]
#             else:
#                 indices = [s.split('_')[0] for s in indices]
#             indices = list(set(indices))
#             celltype_dic[group_name]=indices
            
#     return (celltype_dic)  

# ## function to covert a dictionary to np.list for histogram
# def combined_list(my_dict):
#     combined_list=[]
#     # Combine all lists into a single large list
#     for sublist in my_dict.values():
#         #print(len(sublist))
#         combined_list.extend(sublist)
#     return(combined_list)

# ## within-cluster distance across modalities
# def get_wiclust_dist(pdist_df,celltype_dic):
#     wiclust_dist_dic = {}
#     for celltype, celltype_bc in celltype_dic.items():
#         bc_rna = [bc + '_rna' for bc in celltype_bc]
#         bc_atac=[bc + '_atac' for bc in celltype_bc]
#         celltype_pdist1 = pdist_df.loc[pdist_df.index.isin(bc_rna),pdist_df.columns.isin(bc_atac) ]
#         #celltype_pdist1 = celltype_pdist1.loc[:, celltype_pdist1.columns.isin(bc_atac)]
        
#         #celltype_pdist2 = pdist_df.loc[pdist_df.index.isin(bc_atac) ,pdist_df.columns.isin(bc_rna)]
#         #celltype_pdist2 = celltype_pdist2.loc[:, celltype_pdist2.columns.isin(bc_rna)]
                
#         #pdist = celltype_pdist.values.tolist()
#         pdist1=celltype_pdist1.stack().tolist()
#         #pdist2=celltype_pdist2.stack().tolist()
#         #pdist=pdist1 + pdist2
#         ## save result to dict
#         wiclust_dist_dic[celltype]=pdist1
      
#     ## convert dic to list
#     wiclust_dist = combined_list(wiclust_dist_dic)
#     return(wiclust_dist)

# ## between-cluster distance across modalities
# def get_btclust_dist(pdist_df, celltype_dic):
#     btclust_dist_dic = {}

#     for celltype, celltype_bc in celltype_dic.items():
#         bc_rna = [bc + '_rna' for bc in celltype_bc]
#         bc_atac = [bc + '_atac' for bc in celltype_bc]
#         celltype_pdist1 = pdist_df.loc[pdist_df.index.isin(bc_rna),pdist_df.columns.str.endswith('_atac') ]
#         #celltype_pdist1 = celltype_pdist1.loc[:, celltype_pdist1.columns.str.endswith('_atac') ]
#         celltype_pdist1 = celltype_pdist1.loc[:, ~celltype_pdist1.columns.isin(bc_atac) ]
        
#         #celltype_pdist2 = pdist_df.loc[pdist_df.index.isin(bc_atac) ]
#         #celltype_pdist2 = celltype_pdist2.loc[:, celltype_pdist2.columns.str.endswith('_rna') ]
#         #celltype_pdist2 = celltype_pdist2.loc[:, ~celltype_pdist2.columns.isin(bc_atac) ]
        
#         pdist1=celltype_pdist1.stack().tolist()
#         #pdist2=celltype_pdist2.stack().tolist()
#         #pdist = pdist1 + pdist2
#         ## save result to dict
#         btclust_dist_dic[celltype]=pdist1
#         #dist_bt.append(celltype_pdist.values)
        
#    ## convert dic to list
#     btclust_dist = combined_list(btclust_dist_dic)
#     return(btclust_dist)


#######################################
## function to get within and betwee celltype distance
def lat2pdist(lat_df, bc_test, option):
    #lat_df = lat_df.loc[bc_test]
    lat_df = lat_df.loc[lat_df.index.isin(bc_test)]
    # lat_df_rna = lat_df[lat_df.index.str.endswith('_rna')]
    # lat_df_atac = lat_df[lat_df.index.str.endswith('_atac')]
    # distances = cdist(lat_df_atac, lat_df_rna, metric= option)
    # pdist_df = pd.DataFrame(distances, index = lat_df_atac.index, columns = lat_df_rna.index)
    distances = cdist(lat_df, lat_df, metric= option)
    pdist_df = pd.DataFrame(distances, index = lat_df.index, columns = lat_df.index)
    return (pdist_df)

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

def get_parallel_dist(pdist_df,common_bc):
    parallel_dist = [] 
    for bc in common_bc:
        bc_rna= bc + "_rna"
        bc_atac=bc + "_atac"
        parallel_dist.append(pdist_df.loc[bc_atac, bc_rna]) 

    return(parallel_dist)
    
def get_parallel_dist2(pdist_df,common_bc1, common_bc2):
    parallel_dist1 = [] 
    for bc in common_bc1:
        bc_rna1= bc + "_rna1"
        bc_atac1=bc + "_atac1"
        parallel_dist1.append(pdist_df.loc[bc_atac1, bc_rna1]) 
    parallel_dist2 = [] 
    for bc in common_bc2:
        bc_rna2= bc + "_rna2"
        bc_atac2=bc + "_atac2"
        parallel_dist2.append(pdist_df.loc[bc_atac2, bc_rna2])
            
    combine = parallel_dist1 + parallel_dist2
    return(combine)

def get_parallel_dist3(pdist_df,common_bc1, common_bc2, common_bc3):
    parallel_dist1 = [] 
    for bc in common_bc1:
        bc_rna1= bc + "_rna1"
        bc_atac1=bc + "_atac1"
        parallel_dist1.append(pdist_df.loc[bc_atac1, bc_rna1]) 
    parallel_dist2 = [] 
    for bc in common_bc2:
        bc_rna2= bc + "_rna2"
        bc_atac2=bc + "_atac2"
        parallel_dist2.append(pdist_df.loc[bc_atac2, bc_rna2])

    parallel_dist3 = [] 
    for bc in common_bc3:
        bc_rna3= bc + "_rna3"
        bc_atac3=bc + "_atac3"
        parallel_dist3.append(pdist_df.loc[bc_atac3, bc_rna3])
    combine = parallel_dist1 + parallel_dist2 + parallel_dist3
    return(combine)

    
def get_celltype_dist(pdist_df,label_dic):
    dist_out = pdist_df.apply(lambda column: col_typedist(column, label_dic))
    results= dist_out.apply(lambda row: combine_lists(row), axis=1)
    wiclust_dist= results.loc[0]
    btclust_dist= results.loc[1]
    return([wiclust_dist, btclust_dist])
    
###########################################################
## KNN test function

## with k value, run KNN to predict ATAC label
def knn_ataclabel(k,lat_dfs, pipelines,labels_annot,bc_test_rna, bc_test_atac, res_dist):
    df_preds = []
    df_pred_prob =[]
    for i in range(0, len(lat_dfs)):
        lat_df = lat_dfs[i]
        pipeline = pipelines[i]
        lat_df_rna =lat_df.loc[lat_df.index.isin(bc_test_rna)]
        lat_df_atac = lat_df.loc[lat_df.index.isin(bc_test_atac)]
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

#########################################################

## KNN test function
#######################################################
## get K nearest neighbor distance
def topk_smallest(lst,k):
    return sorted(lst)[:k]

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

def get_celltype_disttopk(pdist_df,label_dic,k):
    dist_out = pdist_df.apply(lambda column: col_typedist(column, label_dic))
    #dist_out_topk = dist_out.applymap(topk_smallest)
    dist_out_topk = dist_out.applymap(lambda lst: topk_smallest(lst, k))
    #results= dist_out.apply(lambda row: combine_lists(row), axis=1)
    results= dist_out_topk.apply(lambda row: combine_lists(row), axis=1)
    wiclust_dist= results.loc[0]
    btclust_dist= results.loc[1]
    return([wiclust_dist, btclust_dist])


#################################################

## helper function for summary statistics

#################################################
## within modality distance functions

def unique_dist(dist_df):
    unique_dists = []
    for i in range(len(dist_df)):
        for j in range(i + 1, len(dist_df)):
            unique_dists.append(dist_df.iloc[i, j])

    return(unique_dists)

def get_samemodality_dist(dist_dfrna, dist_dfatac):
    rna_dist = unique_dist(dist_dfrna)
    atac_dist = unique_dist(dist_dfatac)
    modality_dist = rna_dist + atac_dist
    return(modality_dist)