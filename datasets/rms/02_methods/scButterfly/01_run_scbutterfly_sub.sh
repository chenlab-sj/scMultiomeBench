#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_run_scbutterfly
#BSUB -q dgx
#BSUB -R rusage[mem=150000]
#BSUB -gpu "num=1"
#BSUB -o run_scbutterfly.log
#BSUB -e run_scbutterfly.err

# scButterfly rep1 on RMS Mast607A (GPU). gpu_interactive for wall-time headroom (avoid the skip-retrain
# pitfall: a tiny time= means it loaded a stale model -> rm scbutterfly/model and re-run).
module load conda3/202311
source activate scbutterfly-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/RMS/Mast607/script/scbutterfly
python 01_run_scbutterfly.py
conda deactivate
echo "DONE: latent.csv"
