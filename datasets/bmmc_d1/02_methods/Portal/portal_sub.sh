# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P portal
#BSUB -J HT163-S1H6
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o portal.log
#BSUB -e portal.err

cp 01_run_portal.py ${TOOLS_ROOT}/portal/Portal
cd ${TOOLS_ROOT}/portal/Portal
module load conda3/202402
conda activate portal
python 01_run_portal.py
conda deactivate

