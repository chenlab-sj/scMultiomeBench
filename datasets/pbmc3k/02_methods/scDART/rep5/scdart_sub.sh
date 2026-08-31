# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P scDART
#BSUB -J pbmc3k_scdart_rep5
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o scDART.log
#BSUB -e scDART.err

# scDART pbmc3k reproducibility rep5 (seed random.seed=330 / seeds[330]). Submit from THIS dir:
#   cd scripts/scDART/rep5 && bsub < run_scdart_sub.sh
# 01_run_scdart.py (here) has out_dir hardcoded to scripts/scDART/rep5/ -> latent.csv lands beside it.
cp 01_run_scdart.py ${TOOLS_ROOT}/scDART/scDART/Examples
cd ${TOOLS_ROOT}/scDART/scDART/Examples
module load conda3/202402
conda activate scdart
python 01_run_scdart.py
conda deactivate
