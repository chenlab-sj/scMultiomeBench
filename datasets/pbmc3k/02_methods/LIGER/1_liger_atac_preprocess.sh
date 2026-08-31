#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P rliger
#BSUB -J pbmc10k_rliger_atacpre
#BSUB -q standard
#BSUB -R rusage[mem=200000]
#BSUB -n 1
#BSUB -o 1_liger_atac.log
#BSUB -e 1_liger_atac.err

module load bedops/2.4.2
## preprocessing for liger: transorm scATAC-seq data to gene level counts
## ref: http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html\

## gene body and promoter index are prepared using '${TOOLS_ROOT}/liger/pbmc10k/0_liger_hg20_promoter_bedmake.R'
## gene body and prmoter has been sorted 
#sort -k 1,1 -k2,2n -k3,3n hg20_genes.bed > hg20_genes.sort.bed
#sort -k 1,1 -k2,2n -k3,3n hg20_promoters.bed > hg20_promoters.sort.bed

atac_frag='${DATA_ROOT}/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv'
hg20_gene_index='${TOOLS_ROOT}/liger/liger_need/hg20_genes.sort.bed'
hg20_promoter_index='${TOOLS_ROOT}/liger/liger_need/hg20_promoters.sort.bed'
out_dir='${TOOLS_ROOT}/liger/pbmc10k/liger_need'

sort -k1,1 -k2,2n -k3,3n $atac_frag > $out_dir/atac_fragments.sort.bed

bedmap --ec --delim "\t" --echo --echo-map-id $hg20_promoter_index $out_dir/atac_fragments.sort.bed > $out_dir/atac_promoters_bc.bed
bedmap --ec --delim "\t" --echo --echo-map-id $hg20_gene_index $out_dir/atac_fragments.sort.bed > $out_dir/atac_genes_bc.bed