#!/bin/bash
# =============================================================================
# setup_multigrate_env.sh -- Multigrate (Theis lab) env on St. Jude LSF.
# Repo: https://github.com/theislab/multigrate   import: `import multigrate as mtg`
#
# OPPOSITE of MIRA/MaxFuse: Multigrate 1.0.1 needs a MODERN stack --
# python>=3.12, scvi-tools>=1.4.3, anndata>=0.12, and a GPU (it's a scvi-tools
# VAE). Deps are loose lower-bounds only, so no pin conflicts -- a clean fresh env.
#
# RUN WITH bash on a GPU node (so torch.cuda is testable):
#   bsub -P benchmark -q gpu -gpu "num=1" -Is /bin/bash
#   bash setup_multigrate_env.sh
#   bash verify_multigrate_env.sh
# =============================================================================
set -euo pipefail

ENV_NAME="multigrate-env"
export PYTHONNOUSERSITE=1          # the ~/.local leak fix (as with MIRA/MaxFuse)

module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"

conda env remove -y -n "${ENV_NAME}" 2>/dev/null || true
# py3.12 from conda-forge (the module's base conda is old; -c conda-forge has 3.12)
conda create -y -n "${ENV_NAME}" -c conda-forge python=3.12 pip
conda activate "${ENV_NAME}"
python -m pip install --upgrade pip setuptools wheel

# GPU torch FIRST so scvi-tools sees it satisfied. A100/driver 595 -> cu121 is safe.
# (PyPI torch wheels are CUDA-enabled too, but pin the build explicitly to be sure.)
python -m pip install torch --index-url https://download.pytorch.org/whl/cu121

# Multigrate pulls scvi-tools>=1.4.3, anndata>=0.12, scanpy, jax, etc.
python -m pip install multigrate
python -m pip install muon          # tutorial uses muon/mudata for modality handling

echo "================ pip check ================"
python -m pip check || true   # scvi-tools dep graph is large; warnings here are usually benign

python - <<'PY'
import torch, multigrate, scvi, anndata
print("multigrate :", multigrate.__version__)
print("scvi-tools :", scvi.__version__)
print("anndata    :", anndata.__version__)
print("torch      :", torch.__version__)
print("cuda avail :", torch.cuda.is_available())
PY
echo "DONE: env '${ENV_NAME}' built. Now run:  bash verify_multigrate_env.sh"
