# NOTE: paths below are placeholders. See config/config.yaml and the README.
import process_db
import h5py
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import random
random.seed(1)

rna_h5_files = ["/path/to/data/HTAN/HT137B1-S1H7/scjoint/rna_scjoint.h5"] 
rna_label_files = ["/path/to/data/HTAN/HT137B1-S1H7/scjoint/testrna_celltype_scjoint.csv"] # csv file

atac_h5_files = ["/path/to/data/HTAN/HT137B1-S1H7/scjoint/atac_gene_scjoint.h5"]
atac_label_files = []
process_db.data_parsing(rna_h5_files, atac_h5_files)

rna_label = pd.read_csv(rna_label_files[0], index_col = 0)
rna_label
print(rna_label.value_counts(sort = False))
process_db.label_parsing(rna_label_files, atac_label_files)

