# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import process_db
import h5py
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import random
random.seed(1)




rna1_h5_files = ["/path/to/data/BMMC_d1/scjoint2/rna1_scjoint.h5"] 
rna1_label_files = ["/path/to/data/BMMC_d1/scjoint2/test1rna_celltype_scjoint.csv"] # csv file
atac1_h5_files = ["/path/to/data/BMMC_d1/scjoint2/atac1_gene_scjoint.h5"]
atac1_label_files = []

rna2_h5_files = ["/path/to/data/BMMC_d1/scjoint2/rna2_scjoint.h5"] 
rna2_label_files = ["/path/to/data/BMMC_d1/scjoint2/test2rna_celltype_scjoint.csv"] # csv file
atac2_h5_files = ["/path/to/data/BMMC_d1/scjoint2/atac2_gene_scjoint.h5"]
atac2_label_files = []

rna3_h5_files = ["/path/to/data/BMMC_d1/scjoint2/rna3_scjoint.h5"] 
rna3_label_files = ["/path/to/data/BMMC_d1/scjoint2/test3rna_celltype_scjoint.csv"] # csv file
atac3_h5_files = ["/path/to/data/BMMC_d1/scjoint2/atac3_gene_scjoint.h5"]
atac3_label_files = []


process_db.data_parsing(rna1_h5_files, atac1_h5_files)
process_db.data_parsing(rna2_h5_files, atac2_h5_files)
process_db.data_parsing(rna3_h5_files, atac3_h5_files)


rna1_label = pd.read_csv(rna1_label_files[0], index_col = 0)
#rna_label
print(rna1_label.value_counts(sort = False))
rna2_label = pd.read_csv(rna2_label_files[0], index_col = 0)
#rna_label
print(rna2_label.value_counts(sort = False))
rna3_label = pd.read_csv(rna3_label_files[0], index_col = 0)
#rna_label
print(rna3_label.value_counts(sort = False))
process_db.label_parsing(rna1_label_files, atac1_label_files)
process_db.label_parsing(rna2_label_files, atac2_label_files)
process_db.label_parsing(rna3_label_files, atac3_label_files)

