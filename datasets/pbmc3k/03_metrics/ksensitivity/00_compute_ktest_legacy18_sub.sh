#!/bin/bash
#BSUB -P benchmark
#BSUB -J knn_test
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o knn_test.log
#BSUB -e knn_test.err


module load conda3/202303 

conda activate benchmark_env

python 00_compute_ktest_legacy18.py

conda deactivate
