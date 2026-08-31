#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J bmmc_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=200000]
#BSUB -n 1
#BSUB -o bmmc_metrics.log
#BSUB -e bmmc_metrics.err

# BMMC_d1 stage 2: batch-aware metrics (sample_asw + ks_sample + 3-batch distances) for the 18 staged
# methods. 01_prep_latents.py stages the latents; 02_run_bmmc_metrics.py imports the original benchmark_metrics()
# and scores them -> sum_metrics.csv + celltype_metrics.csv. benchmark_env (scib + benchmark_fun); the
# 3-batch pairwise-distance step is the slow/heavy part (hence large_mem).
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BMMC_d1/benchmark"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}

python 01_prep_latents.py            # stage 18 methods -> staged/<method>/latent.csv
python 02_run_bmmc_metrics.py        # batch-aware scoring -> sum_metrics.csv + celltype_metrics.csv

conda deactivate
echo "DONE: sum_metrics.csv + celltype_metrics.csv (then plot: stage 3)"
