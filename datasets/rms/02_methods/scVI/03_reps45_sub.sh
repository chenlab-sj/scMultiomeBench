#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scvi_reps2
#BSUB -q large_mem
#BSUB -R rusage[mem=150000]
#BSUB -n 4
#BSUB -o run_scvi_reps2.log
#BSUB -e run_scvi_reps2.err

# scVI EXTRA seeds (rep4=13, rep5=99) -- the seed-fragility spot-check (same seeds as scDART rep4/rep5),
# to see whether scVI also has a collapse seed that 3 reps missed. -> res_Mast607A/rep4,rep5/latent.csv.
# rep1/rep2/rep3 untouched. CPU (no GPU) -> large_mem. AFTER this, run 04_fix_scvi_barcodes.py to rebuild
# the <bc>-1_rna/_atac scvi_latent.csv for rep4/rep5 (it now covers rep2..rep5).
module load conda3/202210
module load gcc/12.2.0
conda activate scvi-env
cd ${PROJECT_ROOT}/RMS/Mast607/script/scVI
for rs in "rep4 13" "rep5 99"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scVI $REP (seed $SEED) ==="; python 01_run_scvi.py
done
conda deactivate
echo "DONE: res_Mast607A/rep4,rep5/latent.csv  (then: python 04_fix_scvi_barcodes.py)"
