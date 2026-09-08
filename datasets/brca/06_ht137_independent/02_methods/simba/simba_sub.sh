#BSUB -P simba
#BSUB -J HT163_S1H6
#BSUB -q standard
#BSUB -R rusage[mem=800001]
#BSUB -n 4
#BSUB -o simba.log
#BSUB -e simba.err


module load conda3/202311
source activate simba_env 
python run_simba.py
conda deactivate
