#!/bin/bash
#BSUB -P cobolt
#BSUB -J pbmc3k_cobolt_rep4
#BSUB -q gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o pbmc3k_test_all_coboltgpu.log
#BSUB -e pbmc3k_test_all_coboltgpu.err
# Cobolt pbmc3k reproducibility rep4 (random.seed=100; rep1-3 used 420/0/40). Self-contained: absolute
# inputs, relative out_dir=res_pbmc3k -> writes res_pbmc3k/latent.csv here, then strips the Cobolt
# <dataset>~ barcode prefix to make res_pbmc3k/cobolt_latent.csv (the file reproduce_extrareps.py
# reads). Submit from THIS dir; Cobolt reps are independent so rep4/rep5 can run in parallel.
module load conda3/202210
conda activate cobolt
python pbmc3k_test_all_cobolt.py
python -c "import pandas as pd; d=pd.read_csv('res_pbmc3k/latent.csv',index_col=0); d.index=[str(b).split('~')[-1] for b in d.index]; d.to_csv('res_pbmc3k/cobolt_latent.csv')"
conda deactivate
echo "DONE: res_pbmc3k/cobolt_latent.csv"
