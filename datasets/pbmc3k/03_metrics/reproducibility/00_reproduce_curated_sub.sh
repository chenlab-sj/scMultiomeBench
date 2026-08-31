#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J pbmc3k_fig4_curated
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o fig4_curated.log
#BSUB -e fig4_curated.err
# Redone pbmc3k Fig4 with the CURATED 3-rep selection:
#   scDART = rep3/4/5 (least-reproducible 3, incl the rep4 collapse)
#   Cobolt = rep2/3/4 (least-reproducible 3; note Cobolt is still ~0.98 on Fig4A -- stable in label transfer)
#   scBridge = rep1/3/4, scJoint = rep1/4/5, MIDAS = rep1/4/5 (most-reproducible 3)
#   all other methods = original rep1-3
# Writes to fig4/curated/ so the canonical figures are untouched. Set FIX_RES=1 to also fix the louvain
# resolution per method for Fig4B (reduces the 4A/4B discordance).
module load conda3/202303
conda activate benchmark_env
ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc/pbmc3k/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
export REPRO_OUT="curated"
export LABEL=label.csv
export FIX_RES=1                   # fix louvain resolution per method (from rep1) -> cleaner Fig4B, less 4A/4B discordance
mkdir -p curated
python 00_reproduce_curated.py       # -> curated/{sd_df,nmi_df,reproducibility,...}.csv
python plot_fig4.py               # -> curated/{fig4a_reproducibility,fig4b_nmi}.pdf
conda deactivate
echo "DONE: curated/fig4a_reproducibility.pdf + curated/fig4b_nmi.pdf"
