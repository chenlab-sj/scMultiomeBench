#!/bin/bash
# Confirms midas-env works. Run on a GPU node. "ALL CHECKS PASSED" == success.
set -euo pipefail

export PYTHONNOUSERSITE=1
module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate midas-env

echo "===== imports / API / CUDA ====="
python - <<'PY'
import sys, torch, mudata, anndata, scanpy
import scmidas
print("python  :", sys.version.split()[0])
print("scmidas :", scmidas.__version__, "from", scmidas.__file__)
print("torch   :", torch.__version__)

assert "/envs/midas-env/" in scmidas.__file__, scmidas.__file__
# resolve the MIDAS class (top-level or under .model) and check the API used by run_midas.py
MIDAS = getattr(scmidas, "MIDAS", None) or getattr(getattr(scmidas, "model", object), "MIDAS", None)
assert MIDAS is not None, "scmidas.MIDAS not found (check scmidas.model.MIDAS)"
for m in ("setup_mudata", "train", "get_latent_representation"):
    assert hasattr(MIDAS, m), f"MIDAS.{m} missing"
print("scmidas.MIDAS + setup_mudata/train/get_latent_representation: OK")

ok = torch.cuda.is_available()
print("torch.cuda.is_available():", ok)
if ok:
    print("cuda device:", torch.cuda.get_device_name(0))
assert ok, "CUDA not visible -- are you on a gpu/dgx node? (MIDAS is a GPU method)"
print("\nALL CHECKS PASSED")
PY
