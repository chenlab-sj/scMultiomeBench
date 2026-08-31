
import os
import torch
import numpy as np
import scanpy as sc
import seaborn as sns
import scanpy.external as sce
from torch import nn
from model_utils import feature_prototype_similarity, gmm
from sklearn.metrics import accuracy_score, f1_score, silhouette_score
import anndata
import pickle
import pandas as pd

file_path = "tmp_saved_objects.pkl" 
with open(file_path, 'rb') as file:
    loaded_objects = pickle.load(file)

feature_vec=loaded_objects[0]
pred_vec=loaded_objects[1]
reliability_vec=loaded_objects[2]
label_map=loaded_objects[3]
type_num=loaded_objects[4]
source_adata=loaded_objects[5]
target_adata=loaded_objects[6]
args=loaded_objects[7]

adata = sc.AnnData(feature_vec)
adata.obs["Domain"] = np.concatenate(
    (source_adata.obs["Domain"], target_adata.obs["Domain"]), axis=0
)
sc.tl.pca(adata)
sce.pp.harmony_integrate(adata, "Domain", theta=0.0, verbose=False)
feature_vec = adata.obsm["X_pca_harmony"]
source_adata.obsm["Embedding"] = feature_vec[: len(source_adata.obs["Domain"])]
target_adata.obsm["Embedding"] = feature_vec[len(source_adata.obs["Domain"]) :]
predictions = np.empty(len(target_adata.obs["Domain"]), dtype=np.dtype("U30"))

for k in range(type_num):
    predictions[pred_vec == k] = label_map[k]

if args.novel_type:
    predictions[pred_vec == -1] = "Novel (Most Unreliable)"

target_adata.obs["Prediction"] = predictions
target_adata.obs["Reliability"] = reliability_vec

source_save_path = args.data_path + args.source_data[:-5] + "-integrated.h5ad"
target_save_path = args.data_path + args.target_data[:-5] + "-integrated.h5ad"

source_adata.write(source_save_path)
target_adata.write(target_save_path)

## write cell embeddings

source_lat_df = pd.DataFrame(source_adata.obsm["Embedding"],source_adata.obs.index)
target_lat_df = pd.DataFrame(target_adata.obsm["Embedding"],target_adata.obs.index)
combined_lat_df = pd.concat([source_lat_df,target_lat_df],axis=0)

lat_save_path=args.data_path + "latent.csv"
combined_lat_df.to_csv(lat_save_path)
label_pred_df = pd.DataFrame(target_adata.obs)
label_pred_save_path=args.data_path + "label_pred.csv"
label_pred_df.to_csv(label_pred_save_path)





