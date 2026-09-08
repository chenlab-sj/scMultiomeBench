#!/bin/bash
# Submit the 3 NEW HT137 method runs (MaxFuse, MIDAS, scButterfly) + the scVI re-runs (seeds 0/40),
# to bring HT137B1-S1H7 up to the same 14-method set as HT243 (figS4a). Run on the cluster login node.
#   bash run_methods.sh
# These are heavy GPU/large_mem jobs. All write latent.csv into the method's own folder; the
# figS4b metrics pipeline reads those latents afterward. Main runs only (seed 420) per your choice;
# scVI gets the extra seeds you asked for.
cd "$(dirname "$0")" || exit 1
id() { awk -F'[<>]' '/^Job </{print $2; exit}'; }

# MaxFuse: Signac gene-activity prep (large_mem) -> Fusor (CPU). run must wait for prep.
MP=$(bsub < maxfuse/prep_sub.sh | id);            echo "maxfuse prep    -> $MP"
MR=$(bsub -w "done($MP)" < maxfuse/maxfuse_sub.sh | id); echo "maxfuse run     -> $MR  (after $MP)"

# MIDAS (dgx GPU) + scButterfly (gpu_interactive) -- independent
MI=$(bsub < midas/midas_sub.sh | id);             echo "midas           -> $MI"
SB=$(bsub < scbutterfly/scbutterfly_sub.sh | id); echo "scbutterfly     -> $SB"

# scVI re-runs: seed 0 (rep2) + seed 40 (rep3). The main seed-420 run already exists.
S2=$(bsub < scVI/rep2/scvi_sub.sh | id);          echo "scVI seed0 rep2 -> $S2"
S3=$(bsub < scVI/rep3/scvi_sub.sh | id);          echo "scVI seed40 rep3-> $S3"

echo ""
echo "submitted. watch:  bjobs -w"
echo "outputs (latent.csv each): maxfuse/  midas/  scbutterfly/  scVI/rep2/  scVI/rep3/"
echo "then build/run the metrics pipeline:  cd ../benchmark/figS4b && bash run_all.sh"
