#!/bin/bash
# =============================================================================
# setup_midas_env.sh -- MIDAS (scmidas, labomics, Nat Biotech 2024) on St. Jude LSF.
# Repo: https://github.com/labomics/midas   pip: scmidas   import: `import scmidas`
#
# Modern GPU stack (like Multigrate): python>=3.10, torch 2.5-2.10, lightning,
# scanpy<2, anndata<0.13, mudata<0.4. Generous upper bounds -> no pin hell.
#
# RUN WITH bash on a GPU node (so torch.cuda verifies):
#   bsub -P benchmark -q gpu -gpu "num=1" -Is /bin/bash
#   bash setup_midas_env.sh
#   bash verify_midas_env.sh
# =============================================================================
set -euo pipefail

ENV_NAME="midas-env"
export PYTHONNOUSERSITE=1          # ~/.local leak fix (as with the other methods)

module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"

conda env remove -y -n "${ENV_NAME}" 2>/dev/null || true
conda create -y -n "${ENV_NAME}" -c conda-forge python=3.12 pip
conda activate "${ENV_NAME}"
python -m pip install --upgrade pip setuptools wheel

# GPU torch trio FIRST, matched versions from one CUDA index (A100/driver595 -> cu121).
# scmidas needs torch/torchvision/torchaudio in 2.5-2.10 / 0.20-0.25 ranges; the cu121
# index serves the matched 2.5.1 / 0.20.1 / 2.5.1 trio.
python -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

python -m pip install scmidas      # pulls lightning, scanpy, anndata, mudata, etc.

echo "================ pip check ================"
python -m pip check || true

python - <<'PY'
import torch, scmidas, mudata, scanpy, anndata
print("scmidas :", scmidas.__version__)
print("torch   :", torch.__version__)
print("mudata  :", mudata.__version__)
print("anndata :", anndata.__version__)
print("cuda    :", torch.cuda.is_available())
PY
echo "DONE: env '${ENV_NAME}' built. Now run:  bash verify_midas_env.sh"
