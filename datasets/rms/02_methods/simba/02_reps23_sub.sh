#BSUB -P simba
#BSUB -J HT163_S1H6
#BSUB -q large_mem
#BSUB -R rusage[mem=800001]
#BSUB -n 4
#BSUB -o run_simba_reps.log
#BSUB -e run_simba_reps.err


module load conda3/202311
source activate simba_env

for rs in "rep2 7" "rep3 42"; do
    set -- $rs
    export REP=$1 SEED=$2 TRIM_CUTOFF=0   # 0 = no edge trimming -> keep ALL cells, matching the original's
                                          # balanced 8291/8291 (zero loss). 0.1 keeps most but can drop a few -> asymmetry.
    echo "=== simba $REP (seed $SEED, cutoff $TRIM_CUTOFF) ==="
    python 01_run_simba.py
done

conda deactivate
