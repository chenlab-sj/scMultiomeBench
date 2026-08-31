# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P portal
#BSUB -J portal_archr_pbmc3k
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o portal_archr.log
#BSUB -e portal_archr.err

GA=${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/portal
module load conda3/202402
conda activate portal
# 1) build the ArchR gene-activity atac_gene h5 (in the variant dir; absolute paths inside).
cd "${GA}"
python make_atac_gene_archr.py
# 2) run portal FROM the portal package dir (01_run_portal.py does `import portal`). 01_run_portal.py uses
#    absolute atac_gene_h5 + out_dir, so outputs still land in the variant dir.
cp "${GA}/run_portal.py" ${TOOLS_ROOT}/portal/Portal/
cd ${TOOLS_ROOT}/portal/Portal
python 01_run_portal.py
conda deactivate
