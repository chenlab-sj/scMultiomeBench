
#BSUB -P simba
#BSUB -J pbmc3k
#BSUB -q large_mem 
#BSUB -R rusage[mem=200001]
#BSUB -n 4
#BSUB -o simba.log
#BSUB -e simba.err


module load conda3/202402
source activate env_simba
python 01_run_simba.py
conda deactivate

