# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"
#BSUB -P scDART
#BSUB -J pbmc10k
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o scDART.log
#BSUB -e scDART.err

cp run_scdart.py ${TOOLS_ROOT}/scDART/scDART/Examples
cd ${TOOLS_ROOT}/scDART/scDART/Examples
module load conda3/202402
conda activate scdart
python run_scdart.py
conda deactivate
