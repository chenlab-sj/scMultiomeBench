#BSUB -P unioncom
#BSUB -J pbmc3k
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o unioncom.log
#BSUB -e unioncom.err


module load conda3/202402
conda activate unioncom
python run_unioncom.py
conda deactivate
