
#BSUB -P scMomaT
#BSUB -J pbmc10k
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o scMoMaT.log
#BSUB -e scMoMaT.err


module load conda3/202402
conda activate scmomat
python run_scmomat.py
conda deactivate
