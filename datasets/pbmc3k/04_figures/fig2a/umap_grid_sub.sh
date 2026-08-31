#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J fig2a_umap_grid
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o umap_grid.log
#BSUB -e umap_grid.err

# Fig2A: 5-method UMAP grid (BindSC removed; rows ordered by Fig2b rank ->
# scglue(multiome), scVI, Seurat(CCA), scJoint, Conos). make_umap_grid.py reads the per-method
# latents + label + KNN predictions from the paths hardcoded near the top of that script
# (BASE=${TOOLS_ROOT}); edit that LATENTS/LABEL/KNN block if any path has moved.
# Needs scanpy (neighbors + UMAP per method) -> benchmark_env.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc/pbmc3k/benchmark/fig2a"

python make_umap_grid.py .        # -> fig2a/umap.png

conda deactivate
echo "DONE: fig2a/umap.png"
