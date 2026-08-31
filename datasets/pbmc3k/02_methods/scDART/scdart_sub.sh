
#BSUB -P scDART
#BSUB -J pbmc3k
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o scDART.log
#BSUB -e scDART.err


module load conda3/202402
conda activate scdart
python 01_run_scdart.py
conda deactivate

