#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P velocity
#BSUB -J rms_velocyto
#BSUB -q large_mem
#BSUB -n 8
#BSUB -R rusage[mem=150000]
#BSUB -R span[hosts=1]
#BSUB -o velocyto.log
#BSUB -e velocyto.err

# RMS Mast607A RNA-velocity inputs: velocyto over the GEX BAM -> spliced/unspliced .loom (the long
# confirmation job for the ATAC-precedes-RNA analysis; the fast DPT result already stands on its own).
# Matches the lab's proven velocyto recipe (conda3/202402 + velocity env; unload the default sjcb samtools,
# load samtools/1.17).
#
# *** REFERENCE: this RMS data is hg19 / GRCh37 (BAM chr1 = 249,250,621) -> use the hg19 gene GTF, NOT GRCh38. ***
# No hg19 repeat mask is available (only hg38/mm10 under ~/velocity), so -m is omitted; add `-m <hg19_rmsk.gtf>`
# if you build/obtain one. (-n bumped to 8 vs the GEMM template because this GEX BAM is ~63 GB.)
module load conda3/202402
source activate velocity
module unload sjcb/samtools/1.3.1 2>/dev/null
module load samtools/1.17

OUTS=${DATA_ROOT}/RMS/Mast607A_TB19_22652/Mast607A_TB19_22652/outs
OUTPATH=${PROJECT_ROOT}/RMS/benchmark/atac_precede_rna/velocyto
GTF=${REF_ROOT}/refdata-cellranger-hg19-1.2.0/genes/genes.gtf   # hg19, matches the BAM

mkdir -p "$OUTPATH"
cp "$OUTS/filtered_feature_bc_matrix/barcodes.tsv" "$OUTPATH/barcodes.tsv"
# symlink the 63 GB BAM into OUTPATH so velocyto writes cellsorted_*.bam here (avoids a 63 GB copy; keeps the source dir clean)
ln -sf "$OUTS/gex_possorted_bam.bam" "$OUTPATH/gex_possorted_bam.bam"

velocyto run -b "$OUTPATH/barcodes.tsv" -o "$OUTPATH" -e Mast607A -@ 8 \
  "$OUTPATH/gex_possorted_bam.bam" "$GTF"

conda deactivate
echo "DONE velocyto -> $OUTPATH/Mast607A.loom"
