#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_fig5b_pileups
#BSUB -q large_mem
#BSUB -n 2
#BSUB -R rusage[mem=100000]
#BSUB -o fig5b_pileups.log
#BSUB -e fig5b_pileups.err

# Fig5B: MYOD1 + FOXO1 region pileup tracks for the 3 new methods (MaxFuse/MIDAS/scButterfly), same recipe
# as the published ataclabel2peak.R. seurat4 env (Signac/Seurat). Writes new_methods_{myod1,foxo1}value.csv.
module load conda3/202210
conda activate seurat4

ROOT=${PROJECT_ROOT}
cd "${ROOT}/RMS/benchmark/fig5b"

Rscript 04_compute_new_pileups.R

conda deactivate
echo "DONE: new_methods_myod1value.csv + new_methods_foxo1value.csv"
