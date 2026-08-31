#!/bin/bash
# Submit the ArchR-gene-activity variant runs (R2.1). Run on the cluster login node: bash run_archr_variants.sh
# Each method is identical to its Signac baseline EXCEPT gene activity = the precomputed ArchR GeneScoreMatrix
# (via archr_GeneActivity.R). Latents land in <method>/ ; then score them in metrics/ and diff vs fig2b (Signac).
#
# NOTE two methods need a manual touch (printed at the end):
#   - Portal : run the small h5 builder first (needs scanpy/scipy/h5py).
#   - scJoint: stage-0 prints the ArchR common-gene count -> set 02_config.py input_size before training (like cross-donor).
cd "$(dirname "$0")" || exit 1
id()     { awk -F'[<>]' '/^Job </{print $2; exit}'; }
sub()    { ( cd "$1" && bsub < "$2" ) | id; }
subdep() { ( cd "$1" && bsub -w "done($3)" < "$2" ) | id; }

echo "== single-stage R =="
printf "  Seurat_CCA -> %s\n" "$(sub Seurat_CCA pbmc3k_testall_seurat3_sub.sh)"
printf "  bindsc     -> %s\n" "$(sub bindsc pbmc3k_bindsc_sub.sh)"

echo "== maxfuse: prep -> run =="
mp=$(sub maxfuse 00_prep_maxfuse_input_sub.sh);                 echo "  prep -> $mp"
printf "  run  -> %s (after %s)\n" "$(subdep maxfuse 01_run_maxfuse_sub.sh "$mp")" "$mp"

echo "== scBridge: s0 -> s1 -> s2 =="
b0=$(sub scBridge 0_scBridge_sub.sh);          echo "  s0 -> $b0"
b1=$(subdep scBridge 01_run_scbridge.sh "$b0");  echo "  s1 -> $b1 (after $b0)"
printf "  s2 -> %s (after %s)\n" "$(subdep scBridge 2_save_scBridge_sub.sh "$b1")" "$b1"

echo ""
echo "== MANUAL (2 methods) =="
echo "  Portal : cd portal && conda activate portal && python make_atac_gene_archr.py   # build atac_gene_archr.h5"
echo "           then: bsub < run_portal_sub.sh"
echo "  scJoint: ( cd scJoint && bsub < 0_make_scjoint_h5_sub.sh )   # stage 0 prints the ArchR common-gene count"
echo "           -> set self.input_size in scJoint/config.py to that number"
echo "           -> bsub < scJoint/pbmc3k_scjoint_process.sh ; then scJoint main.py ; then assemble scJoint_latent.csv"
echo ""
echo "latents: Seurat_CCA/res_pbmc3k_testall/coembed_coor.csv  bindsc/res_pbmc3k/coembed_coor.csv  maxfuse/latent.csv"
echo "         scBridge/latent.csv  portal/lat_df.csv  scJoint/output/scJoint_latent.csv"
