#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J pbmc3k_scjoint_s3_rep4
#BSUB -q dgx
#BSUB -R rusage[mem=300000]
#BSUB -gpu "num=1"
#BSUB -o scjointgpu_s3.log
#BSUB -e scjointgpu_s3.err

# scJoint pbmc3k reproducibility rep4 (02_config.py seed=100; rep1-3 used 1/40/0). Swaps THIS rep's
# config into the LSA scJoint dir, trains, assembles the latent, and copies ./output (incl.
# scJoint_latent.csv) back here (the repo archive reproduce_extrareps.py reads). Submit from this dir:
#   cd scripts/scJoint/rep4 && bsub < 02_scjoint_s3_sub.sh
# Do NOT run rep4 and rep5 concurrently -- they share the LSA root 02_config.py and ./output.
module load conda3/202210
source activate scjoint
module load gcc/13.1.0

REPO=${PROJECT_ROOT}/pbmc/pbmc3k/scripts
outdir=$(pwd)
cp 02_config.py ${TOOLS_ROOT}/scJoint/
cd ${TOOLS_ROOT}/scJoint
rm -rf output && mkdir -p output
python main.py                                                    # -> output/*_embeddings.txt
SCJOINT_ROOT=${TOOLS_ROOT}/scJoint python "$REPO/scJoint/assemble_scjoint_latent.py"  # -> output/scJoint_latent.csv
cp -r output/ "$outdir"                                           # -> scripts/scJoint/rep4/output/
conda deactivate
