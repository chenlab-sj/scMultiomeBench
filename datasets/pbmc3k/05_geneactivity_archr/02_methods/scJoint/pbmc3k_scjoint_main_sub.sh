#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J ga_scjoint_main
#BSUB -q dgx
#BSUB -R rusage[mem=300000]
#BSUB -gpu "num=1"
#BSUB -o scjointgpu_main_archr.log
#BSUB -e scjointgpu_main_archr.err

# scJoint training (main.py) for the ArchR gene-activity variant, then assemble the latent in the
# SAME job. main.py reads 02_config.py (cp'd into the package) and writes embeddings to the package's
# shared output dir; we snapshot that into ./output (the shared dir is overwritten by any scJoint
# run) and immediately assemble -> output/scJoint_latent.csv. Run AFTER pbmc3k_scjoint_process.sh.
module load conda3/202210
source activate scjoint
module load gcc/13.1.0

GA=${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/scJoint
cd "$GA"
cp 02_config.py ${TOOLS_ROOT}/scJoint/
cd ${TOOLS_ROOT}/scJoint
python main.py

# snapshot the embeddings into the variant dir (race-safe vs other scJoint runs), then assemble
rm -rf "$GA/output"
cp -r output "$GA/output"
cd "$GA"
python scjoint_assemble.py

conda deactivate
echo "DONE: scJoint/output/scJoint_latent.csv"
