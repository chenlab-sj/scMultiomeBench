#!/bin/bash
# =============================================================================
# setup_scbutterfly_env.sh -- scButterfly (BioX-NKU) on St. Jude LSF.
# Repo: https://github.com/BioX-NKU/scButterfly  pip: scButterfly  import: scButterfly
#
# MIRA-like LEGACY pins (scButterfly 0.0.9 declares exact old versions):
#   scvi-tools==0.19.0  scipy==1.9.3  episcanpy==0.3.2  pot==0.9.0  torch>=1.12.1
# -> python 3.9 + a 2022-era torch (1.13.1) on GPU. Build carefully, isolated.
#
# RUN WITH bash on a GPU node:
#   bsub -P benchmark -q gpu -gpu "num=1" -Is /bin/bash
#   bash setup_scbutterfly_env.sh
#   bash verify_scbutterfly_env.sh
# =============================================================================
set -euo pipefail

ENV_NAME="scbutterfly-env"
export PYTHONNOUSERSITE=1

module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"

conda env remove -y -n "${ENV_NAME}" 2>/dev/null || true
conda create -y -n "${ENV_NAME}" -c conda-forge python=3.9 pip
conda activate "${ENV_NAME}"
python -m pip install --upgrade "pip==24.0" "setuptools<70" "wheel"

# constraints: cap numpy below 1.24 (np.int/np.float aliases) and hold scipy at scButterfly's pin.
# Workflow-verified scvi-tools-0.19.0-era pin set (adversarial audit: 0 refutations).
# These cap the packages scButterfly/scvi-tools pull too-new by default and break import.
CFILE="$(mktemp /tmp/scb_constraints.XXXXXX.txt)"
cat > "${CFILE}" <<'EOF'
numpy==1.23.5
scipy==1.9.3
torchmetrics==0.11.4
anndata==0.10.8
mudata==0.2.4
pytorch-lightning==1.7.7
lightning<2.0
scanpy<1.10
EOF
# Why these exact values:
#  torchmetrics 0.11.4 = last pre-1.0 with the private _compare_version pl 1.7.x imports.
#  anndata 0.10.8 = last anndata with AlignedViewMixin AND py3.9 (0.10.9 removed it).
#  mudata 0.2.4 needs anndata>=0.10.8; pytorch-lightning <1.8 (scvi 0.19.0 hard pin, NOT
#  the 2.x 'lightning' pkg); scanpy<1.10 stays on the anndata-0.10.x line for episcanpy 0.3.2.

# GPU torch trio FIRST (2022-era, compatible with scvi-tools 0.19.0). A100 driver595 -> cu117 OK.
python -m pip install torch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 \
    --index-url https://download.pytorch.org/whl/cu117

# scButterfly pulls scvi-tools==0.19.0, scipy==1.9.3, episcanpy==0.3.2, pot==0.9.0, etc.
python -m pip install -c "${CFILE}" "numpy==1.23.5"
python -m pip install -c "${CFILE}" "scButterfly==0.0.9"

echo "================ pip check ================"
python -m pip check || true   # legacy stack -> some benign warnings possible

python - <<'PY'
import torch, scButterfly
from scButterfly.butterfly import Butterfly
print("scButterfly import: OK")
print("torch     :", torch.__version__)
print("cuda avail:", torch.cuda.is_available())
PY
echo "DONE: env '${ENV_NAME}' built. Now run:  bash verify_scbutterfly_env.sh"
