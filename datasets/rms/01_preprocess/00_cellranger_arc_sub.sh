# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P cellranger
#BSUB -J Mast607A_TB19_22652
#BSUB -q rhel8_standard
#BSUB -R rusage[mem=200000]
#BSUB -n 8
#BSUB -o Mast607A_TB19_22652.log
#BSUB -e Mast607A_TB19_22652.err

module load cellranger-arc/2.0.0 

cellranger-arc count --jobmode=lsf \
   --id=Mast607A_TB19_22652 \
   --reference=${REF_ROOT}/hg19/GRCh37 \
   --libraries=${DATA_ROOT}/RMS/Mast607A_TB19_22652/Mast607A.csv
