#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J scjoint_patch
#BSUB -q standard
#BSUB -R rusage[mem=4000]
#BSUB -n 1
#BSUB -o scjoint_patch.log
#BSUB -e scjoint_patch.err

# Auto-patch 02_config.py's shape constants from the files s1/s2 produced. Run AFTER s2, BEFORE s3.
# number_of_class = # cell-type classes (from label_to_idx.txt); input_size = # common genes (from s1's input_size.txt).
D=${PROJECT_ROOT}/BMMC_d1/s1d1_paired/scjoint
# NOT `cd "$(dirname "$0")"`: with `bsub < script` LSF runs a SPOOLED COPY of this file under
# ~/.lsbatch on the exec host, so $0 is that copy and the cd landed outside the method dir. Both
# seds then failed with "can't read 02_config.py" while the job still reported Successfully completed,
# leaving input_size stale (18366) and crashing s3. LS_SUBCWD is the dir bsub was run from.
cd "${LS_SUBCWD:-$PWD}" || exit 1     # scripts scjoint dir (02_config.py lives here; s3 cp's THIS 02_config.py)
[ -f 02_config.py ] || { echo "ERROR: 02_config.py not found in $PWD"; exit 1; }

NCLASS=$(wc -l < "$D/label_to_idx.txt" 2>/dev/null || sort -u "$D/testrna_celltype_scjoint.txt" | wc -l)
INSIZE=$(cat "$D/input_size.txt")
[ -n "$NCLASS" ] && [ -n "$INSIZE" ] || { echo "ERROR: cannot read NCLASS/INSIZE from $D"; exit 1; }

# Patch only the FIRST occurrence: 02_config.py carries other datasets' branches further down
# (input_size = 18603, 17668) that must not be touched. Matching [0-9]+ keeps this idempotent,
# so re-running after a partial patch still converges.
sed -i "0,/self\.number_of_class = [0-9]\+/s//self.number_of_class = ${NCLASS}/" 02_config.py
sed -i "0,/self\.input_size = [0-9]\+/s//self.input_size = ${INSIZE}/" 02_config.py

# Verify and FAIL the job if the patch did not take. The old version echoed success unconditionally,
# which is why a completely no-op patch looked fine in the log.
GOT=$(grep -m1 -o "self\.input_size = [0-9]\+" 02_config.py | grep -o "[0-9]\+")
[ "$GOT" = "$INSIZE" ] || { echo "ERROR: 02_config.py input_size is '$GOT', expected '$INSIZE'"; exit 1; }
echo "patched 02_config.py: number_of_class=${NCLASS} input_size=${INSIZE}"
