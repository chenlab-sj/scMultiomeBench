# NOTE: paths below are placeholders. See config/config.yaml and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
import os
import time
import math
import torch
import seaborn as sns

import numpy as np
import pandas as pd
import scipy as sp
import anndata as ad

from utils.train_util import Siamese_Trainer
from utils.data_util import dataset, test_data
import matplotlib.pyplot as plt

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = '0'
os.environ["CUDA_LAUNCH_BLOCKING"] = '1'

datapath ="/path/to/tools/MinNet/pbmc3k/train"
mod1_h5 = os.path.join(datapath, "rna_training.h5")
mod2_h5 = os.path.join(datapath, "atac_training.h5")
dat = dataset(mod1=mod1_h5, 
              mod2=mod2_h5,
              batch_number=750, C=3.0, validation=True)



train = Siamese_Trainer(dat, num_epochs=50, batch_number=750, model_path='./ckpts/',
                        num_peaks=3543, num_genes=3543, margin = 10.0,
                        lamb=0.5, lr=1e-4, validation=True)

log_file = 'ckpts/MinNet_margin10.0_lamb_0.5.txt'
with open(log_file, 'w') as fw:
    #fw.write(str(args)+'\n')
    train.train(fw, intervals=1, initial=True)