#!/bin/bash
# =============================================================================
# setup_maxfuse_env.sh -- MaxFuse (Nat Biotech 2024) env on St. Jude LSF.
# Repo: https://github.com/shuxiaoc/maxfuse   import name: `maxfuse`
#
# MaxFuse 0.0.2 (Feb 2023) is CPU-only (sklearn/scipy/scanpy) -- NO torch/GPU.
# Its deps are loose lower-bounds only, so unlike MIRA there are no hard pin
# conflicts. We still build on python 3.8 and cap numpy < 1.24, because the
# Feb-2023 release predates numpy 1.24's removal of np.int/np.float/np.bool
# aliases, which silently break many packages of that era.
#
# RUN WITH bash (not sh), isolation matters (see PYTHONNOUSERSITE below):
#   bash setup_maxfuse_env.sh
#   bash verify_maxfuse_env.sh
# =============================================================================
set -euo pipefail

ENV_NAME="maxfuse-env"

# Isolate from ~/.local user-site (the leak that bit us with MIRA: other methods'
# packages there otherwise pollute/shadow the env). Must be set DURING the build
# so pip installs every dep into the env rather than borrowing from ~/.local.
export PYTHONNOUSERSITE=1

module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"

conda env remove -y -n "${ENV_NAME}" 2>/dev/null || true
conda create  -y -n "${ENV_NAME}" python=3.8 pip
conda activate "${ENV_NAME}"
python -m pip install --upgrade "pip==24.0" "setuptools<70" "wheel"

# Constraints: a coherent py3.8 stack (same versions we verified for MIRA) +
# numpy capped below 1.24. Applied via -c on every install.
CFILE="$(mktemp /tmp/maxfuse_constraints.XXXXXX.txt)"
cat > "${CFILE}" <<'EOF'
numpy==1.23.5
scipy==1.10.1
scikit-learn==1.3.2
scanpy==1.9.8
pandas==1.5.3
EOF

python -m pip install -c "${CFILE}" "numpy==1.23.5"
python -m pip install -c "${CFILE}" "maxfuse==0.0.2"   # pulls scanpy/igraph/leidenalg/matplotlib

echo "================ pip check ================"
python -m pip check

rm -f "${CFILE}"
echo "DONE: env '${ENV_NAME}' built. Now run:  bash verify_maxfuse_env.sh"
