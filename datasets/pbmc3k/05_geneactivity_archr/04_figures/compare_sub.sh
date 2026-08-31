#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J ga_archr_compare
#BSUB -q standard
#BSUB -R rusage[mem=8000]
#BSUB -n 1
#BSUB -o ga_archr_compare.log
#BSUB -e ga_archr_compare.err

# R2.1 final panel: Signac (fig2b) vs ArchR composite. Headline Spearman/scatter over the 5 clean-swap
# methods (Seurat(CCA), bindSC, MaxFuse, scBridge, Portal); scJoint flagged + excluded. Reads
# ../fig2b/sum_metrics_clean.csv + metrics/sum_metrics_clean.csv -> compare_geneactivity.csv + .png.
# Run AFTER the plot job (ga_archr_plot) writes metrics/sum_metrics_clean.csv:
#     bsub -w "done(ga_archr_plot)" < compare_sub.sh
# (compare_geneactivity.py is Python -> benchmark_env, NOT the plot's seurat4.)
module load conda3/202303
conda activate benchmark_env
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr
python compare_geneactivity.py
conda deactivate
echo "DONE: compare_geneactivity.csv + compare_geneactivity.png"
