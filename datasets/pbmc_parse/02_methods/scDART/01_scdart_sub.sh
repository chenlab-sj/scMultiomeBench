# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P scDART
#BSUB -J HT163_S1H6
#BSUB -q dgx
#BSUB -R rusage[mem=50000]
#BSUB -gpu "num=1"
#BSUB -o scDART.log
#BSUB -e scDART.err

cp 01_run_scdart.py ${TOOLS_ROOT}/scDART/scDART/Examples
cd ${TOOLS_ROOT}/scDART/scDART/Examples
module load conda3/202402
conda activate scdart
python 01_run_scdart.py
conda deactivate
