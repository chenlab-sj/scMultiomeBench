# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
"""
BRCA Cobolt REP runner (rep4/rep5) — identical to 01_run_cobolt.py except out_dir + seed come from env and
we also write the prefix-stripped cobolt_latent.csv the reproduce config reads. Two ADDITIONAL seeds so
MIDAS (5 reps) and Cobolt can be compared like-for-like at 5-vs-5. Outputs to the REPO (not {HTAN}) so all
5 reps are mountable for analysis.
  REP=rep4 SEED=100  /  REP=rep5 SEED=200
"""
from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt
import os
from scipy import io
import pandas as pd
import pickle
import timeit
import random

REPO = "/path/to/multiomeBench/BRCA/HT243B1-S1H4"
REP  = os.environ["REP"]
SEED = int(os.environ["SEED"])

# inputs (read-only) stay on the cluster where the data lives
train_data_dir = "/path/to/data/HTAN/HT243B1-S1H4/train_forHT243B1-S1H4/"
test_rna_data_dir = "/path/to/data/HTAN/HT243B1-S1H4/HT243B1-S1H4_rna/"
test_atac_data_dir = "/path/to/data/HTAN/HT243B1-S1H4/HT243B1-S1H4_commonpeak/"
out_dir = os.path.join(f"{REPO}/cobolt", REP)          # -> repo/cobolt/rep4 or rep5
dataset_name = "train"
test_rna_dataset_name = 'test_rna'
test_atac_dataset_name = 'test_atac'

start = timeit.default_timer()
random.seed(SEED)

## load train rna data
count = io.mmread(os.path.join(train_data_dir, "matrix.mtx")).T.tocsr().astype(float)
feature = pd.read_csv(os.path.join(train_data_dir, "features.tsv"), sep="\t", header=None)
barcode = pd.read_csv(os.path.join(train_data_dir, "barcodes.tsv"), header=None, usecols=[0])[0].values.astype('str')

gene_feature = feature[feature[2] == 'Gene Expression']
gene_name = gene_feature[1].values.astype('str')
gene_count = count[:, gene_feature.index]
rna = SingleData("GeneExpr", dataset_name, gene_name, gene_count, barcode)

peak_feature = feature[feature[2] == 'Peaks']
peak_name = peak_feature[1].values.astype('str')
peak_count = count[:, peak_feature.index]
atac = SingleData("Peak", dataset_name, peak_name, peak_count, barcode)

## test data (barcodes renamed to force unmatched -> unpaired)
test_rna_count = io.mmread(os.path.join(test_rna_data_dir, "matrix.mtx")).T.tocsr().astype(float)
test_rna_feature = pd.read_csv(os.path.join(test_rna_data_dir, "features.tsv"), sep="\t", header=None)
test_rna_barcode = pd.read_csv(os.path.join(test_rna_data_dir, "barcodes.tsv"), header=None, usecols=[0])[0].values.astype('str')
test_barcode_rna = [item + "_rna" for item in test_rna_barcode]
test_gene_name = test_rna_feature[1].values.astype('str')
test_rna = SingleData("GeneExpr", test_rna_dataset_name, test_gene_name, test_rna_count, test_barcode_rna)

test_atac_count = io.mmread(os.path.join(test_atac_data_dir, "matrix.mtx")).T.tocsr().astype(float)
test_atac_feature = pd.read_csv(os.path.join(test_atac_data_dir, "features.tsv"), sep="\t", header=None)
test_atac_barcode = pd.read_csv(os.path.join(test_atac_data_dir, "barcodes.tsv"), header=None, usecols=[0])[0].values.astype('str')
test_barcode_atac = [item + "_atac" for item in test_atac_barcode]
test_peak_name = test_atac_feature[1].values.astype('str')
test_atac = SingleData("Peak", test_atac_dataset_name, test_peak_name, test_atac_count, test_barcode_atac)

multi = MultiomicDataset.from_singledata(rna, atac, test_rna, test_atac)
print(multi)

## training (same hyperparams as 01_run_cobolt.py)
model = Cobolt(dataset=multi, lr=0.0001, n_latent=16)
model.train(num_epochs=300)
model.calc_all_latent()
latent = model.get_all_latent()
stop = timeit.default_timer()

lat_df = pd.DataFrame(latent[0], index=latent[1])
lat_df = lat_df.set_axis(["latent_" + s for s in lat_df.columns.astype("str").tolist()], axis="columns")
latent_raw = model.get_all_latent(correction=False)
raw_lat_df = pd.DataFrame(latent_raw[0], index=latent_raw[1])
raw_lat_df = raw_lat_df.set_axis(["latent_" + s for s in raw_lat_df.columns.astype("str").tolist()], axis="columns")

os.makedirs(out_dir, exist_ok=True)
lat_df.to_csv(os.path.join(out_dir, "latent.csv"))
raw_lat_df.to_csv(os.path.join(out_dir, "raw_latent.csv"))
# prefix-stripped version the reproduce config reads (train~/test_*~ removed -> <bc>-1_rna / <bc>-1_atac)
cobolt_lat_df = lat_df.copy()
cobolt_lat_df.index = cobolt_lat_df.index.str.replace(r"^.*~", "", regex=True)
cobolt_lat_df.to_csv(os.path.join(out_dir, "cobolt_latent.csv"))
with open(os.path.join(out_dir, "trained_cobolt.pkl"), "wb") as file:
    pickle.dump((model), file)
with open(os.path.join(out_dir, "runtime.txt"), "w") as file:
    file.write(str(stop - start))
print(f"DONE {REP} (seed {SEED}): wrote latent.csv + cobolt_latent.csv to {out_dir}  time={stop-start:.1f}s")
