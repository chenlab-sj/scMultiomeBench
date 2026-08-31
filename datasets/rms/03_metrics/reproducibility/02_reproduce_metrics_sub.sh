#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_fig4_reproduce
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o fig4_reproduce.log
#BSUB -e fig4_reproduce.err

# RMS Mast607A Fig4A/B reproducibility: KNN-label std (Fig4A) + pairwise-NMI (Fig4B) over 14 methods x 3 reps.
# benchmark_env (sklearn + scanpy + benchmark_fun). The louvain clusterings are the slow part.
# 01_prep_latents.py (re)stages every method's 3 reps -> ../staged/ (normalized to <bc>-1_rna / <bc>-1_atac).
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/RMS/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
export REPRO_OUT="."
export LABEL=label.csv
export MIN_CELLS=100
export FIX_RES=1                  # Fig4B: one louvain resolution per method (from rep1) reused across reps
                                 # -> NMI reflects embedding stability, not per-rep resolution-search noise

python ../prep_latents.py        # (re)stage all methods; scVI/scDART/scJoint/Cobolt/MIDAS now 5 seeds, scButterfly forced
python 02_reproduce_metrics.py
python plot_fig4.py

conda deactivate
echo "DONE: reproducibility.csv + fig4a_reproducibility.pdf + fig4b_nmi.pdf"
