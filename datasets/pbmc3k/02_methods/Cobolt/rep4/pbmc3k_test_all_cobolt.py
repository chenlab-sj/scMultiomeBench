# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt
import os

from scipy import io
import pandas as pd
import pickle
import timeit
import random

#### user input

train_data_dir = "/path/to/data/pbmc3k/pbmc3k_filtered_feature_bc_matrix_train"
test_data_dir = "/path/to/data/pbmc3k/pbmc3k_filtered_feature_bc_matrix_test"
out_dir = "res_pbmc3k"
dataset_name = "pbmc3k_train"
test_rna_dataset_name='test_rna'
test_atac_dataset_name='test_atac'


## load train data
start = timeit.default_timer()
random.seed(100)
count = io.mmread(os.path.join(train_data_dir, "matrix.mtx")).T.tocsr().astype(float)
feature = pd.read_csv(
    os.path.join(train_data_dir, "features.tsv"),sep="\t", header = None)
barcode = pd.read_csv(
    os.path.join(train_data_dir, "barcodes.tsv"),
    header=None, usecols=[0]
)[0].values.astype('str')

### rna
gene_feature=feature[feature[2]=='Gene Expression']
gene_name = gene_feature[1].values.astype('str')
gene_count = count[:,gene_feature.index]
rna = SingleData("GeneExpr", dataset_name, gene_name, gene_count, barcode)

### atac
peak_feature=feature[feature[2]=='Peaks']
peak_name = peak_feature[1].values.astype('str')
peak_count = count[:,peak_feature.index]
atac = SingleData("Peak", dataset_name, peak_name, peak_count, barcode)


###########################################
## load test data and rename barcode to force unmatched barcode
test_count = io.mmread(os.path.join(test_data_dir, "matrix.mtx")).T.tocsr().astype(float)
test_feature = pd.read_csv(
    os.path.join(test_data_dir, "features.tsv"),sep="\t", header = None)
test_barcode = pd.read_csv(
    os.path.join(test_data_dir, "barcodes.tsv"),
    header=None, usecols=[0]
)[0].values.astype('str')

test_barcode_rna=[item + "_rna" for item in test_barcode]
test_barcode_atac=[item + "_atac" for item in test_barcode]

### test rna
test_rna_feature=test_feature[test_feature[2]=='Gene Expression']
test_gene_name = test_rna_feature[1].values.astype('str')
test_rna_count = test_count[:,test_rna_feature.index]
test_rna = SingleData("GeneExpr", test_rna_dataset_name,test_gene_name, test_rna_count, test_barcode_rna)

### test atac
test_atac_feature=test_feature[test_feature[2]=='Peaks']
test_peak_name = test_atac_feature[1].values.astype('str')
test_atac_count = test_count[:,test_atac_feature.index]
test_atac = SingleData("Peak", test_atac_dataset_name,test_peak_name, test_atac_count, test_barcode_atac)

multi = MultiomicDataset.from_singledata(rna,atac,test_rna,test_atac)
print(multi)

## training
model = Cobolt(dataset=multi,lr=0.001, n_latent=16)
model.train(num_epochs=100)

model.calc_all_latent()
latent = model.get_all_latent()
print("------ Done ------")
stop = timeit.default_timer()

## prepare latent dataframe
lat_df = pd.DataFrame(latent[0],index=latent[1])
lat_df = lat_df.set_axis(["latent_" + s  for s in lat_df.columns.astype("str").tolist()],axis="columns")
# Remove "pbmc3x~" from the index
#prefix=dataset_name + "~"
#lat_df.index =lat_df.index.str.replace(prefix, '')

########### RNA and ATAC raw latent
latent_raw= model.get_all_latent(correction=False)
raw_lat_df= pd.DataFrame(latent_raw[0], index=latent_raw[1])
raw_lat_df = raw_lat_df.set_axis(["latent_" + s  for s in raw_lat_df.columns.astype("str").tolist()],axis="columns")
#raw_lat_df['dataset']=latent_raw[2]

## only get paired data
#lat_df_join =raw_lat_df[raw_lat_df['dataset']=='joint']
#lat_df_join = lat_df_join.drop('dataset', axis=1)

#lat_df_gene =raw_lat_df[raw_lat_df['dataset']=='GeneExpr']
#lat_df_gene = lat_df_gene.drop('dataset', axis=1)

#lat_df_peak =raw_lat_df[raw_lat_df['dataset']=='Peak']
#lat_df_peak = lat_df_peak.drop('dataset', axis=1)

## save results:
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
###  latent results    
latent_csv = os.path.join(out_dir,"latent.csv")
lat_df.to_csv(latent_csv)

####  raw latent result
raw_latent_csv = os.path.join(out_dir,"raw_latent.csv")
raw_lat_df.to_csv(raw_latent_csv)

#gene_latent_csv = os.path.join(out_dir,"gene_latent.csv")
#lat_df_gene.to_csv(gene_latent_csv)

#peak_latent_csv = os.path.join(out_dir,"peak_latent.csv")
#lat_df_peak.to_csv(peak_latent_csv)

## cobolt model
cobolt_model = os.path.join(out_dir,"trained_cobolt.pkl")

with open(cobolt_model, "wb") as file:
    pickle.dump((model), file)

print('Time(s): ', stop - start)  
# record time 
runtime_out = os.path.join(out_dir,"runtime.txt")
with open(runtime_out, 'w') as file:
    file.write(str(stop-start))  
