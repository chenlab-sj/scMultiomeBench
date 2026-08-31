#BSUB -P simba
#BSUB -J HT163_S1H6
#BSUB -q standard
#BSUB -R rusage[mem=800001]
#BSUB -n 4
#BSUB -o simba.log
#BSUB -e simba.err


module load conda3/202311
source activate simba_env 
cd "${LS_SUBCWD:-$PWD}" || exit 1
python 01_run_simba.py
rc=$?
conda deactivate
# Propagate the real status: ending on `conda deactivate` made LSF report
# "Successfully completed" even when 01_run_simba.py had crashed.
[ $rc -eq 0 ] && echo "DONE: latent.csv"
exit $rc
