#!/bin/bash
# =============================================================================
# setup_mira_env.sh  --  mira-multiome 2.1.1 on St. Jude LSF (Linux x86_64)
#
# WHY THIS IS PINNED THE WAY IT IS:
#   mira-multiome 2.1.1 has HARD UPPER bounds (verified from PyPI requires_dist):
#       torch<=2.0.0  networkx<=2.5  scipy<=1.10.2  scikit-learn<1.4  sqlalchemy<2
#   So you CANNOT pair it with latest scanpy/torch. The previous failure
#   (torch 2.7.1 / networkx 3.4.2 / scanpy 1.11.5) came from installing the
#   latest of everything. This builds a consistent Python 3.8 / 2022-era stack.
#
# RUN IT WITH bash (NOT sh) on a GPU node so torch.cuda is testable:
#   bsub -P benchmark -q gpu -gpu "num=1" -Is /bin/bash      # interactive GPU node
#   bash setup_mira_env.sh
#   bash verify_mira_env.sh                                   # confirm success
# =============================================================================
set -euo pipefail

ENV_NAME="mira-env"   # matches mira_sub.sh ('source activate mira-env')

# Isolate from ~/.local user-site. The original account had many benchmark-method
# packages there (cobolt, scdart, portal-sc, scvi-tools, tensorflow, ...) installed
# with `pip --user`; without this they leak into AND shadow every py3.8 env
# (this is what hid mira.topics).
export PYTHONNOUSERSITE=1

module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"

# --- 0. wipe the previous BROKEN env so we start clean ----------------------
conda env remove -y -n "${ENV_NAME}" 2>/dev/null || true

# --- 1. bare Python 3.8 env; let pip own the scientific stack ----------------
conda create -y -n "${ENV_NAME}" python=3.8 pip
conda activate "${ENV_NAME}"
python -m pip install --upgrade "pip==24.0" "setuptools<70" "wheel"

# --- 2. constraints file: hard-freeze the load-bearing ceilings --------------
# passed via -c on EVERY install so pip can NEVER bump torch/networkx/numpy.
CFILE="$(mktemp /tmp/mira_constraints.XXXXXX.txt)"
cat > "${CFILE}" <<'EOF'
torch==2.0.0
networkx==2.5
numpy==1.24.4
scipy==1.10.1
scikit-learn==1.3.2
pyro-ppl==1.8.6
matplotlib==3.7.5
numba==0.58.1
llvmlite==0.41.1
umap-learn==0.5.7
scanpy==1.9.8
sqlalchemy==1.4.54
EOF

# --- 3. CUDA torch FIRST (A100 driver 595 -> cu118 is correct) ---------------
# cu118 index serves only +cu118 wheels, so ==2.0.0 -> torch-2.0.0+cu118-cp38.
python -m pip install "torch==2.0.0" --index-url https://download.pytorch.org/whl/cu118

# --- 4. pre-pin the other ceilings so later resolves treat them as satisfied -
python -m pip install -c "${CFILE}" "numpy==1.24.4" "networkx==2.5" "scipy==1.10.1"

# --- 5. MIRA itself (pulls pyro-ppl, anndata, lisa2, mira-moods, optuna, ...) -
python -m pip install -c "${CFILE}" "mira-multiome==2.1.1"

# --- 6. nail pyro-ppl: inside mira's <2 range, 1.8.5 needs torch>=2.0.1 (bad);
#         1.8.6 needs only torch>=1.11.0 -> safe with torch 2.0.0 -------------
python -m pip install -c "${CFILE}" "pyro-ppl==1.8.6"

# --- 7. scanpy analysis stack pinned to the py3.8 + networkx<=2.5 safe set ----
python -m pip install -c "${CFILE}" \
    "scikit-learn==1.3.2" "matplotlib==3.7.5" "numba==0.58.1" "llvmlite==0.41.1" \
    "umap-learn==0.5.7" "scanpy==1.9.8"

# --- 8. keep sqlalchemy under mira's <2 cap (optuna 2.x needs it) -------------
python -m pip install -c "${CFILE}" "sqlalchemy==1.4.54"

# --- 9. consistency gate: MUST print 'No broken requirements found.' ----------
echo "================ pip check ================"
python -m pip check

rm -f "${CFILE}"
echo "DONE: env '${ENV_NAME}' built. Now run:  bash verify_mira_env.sh"
