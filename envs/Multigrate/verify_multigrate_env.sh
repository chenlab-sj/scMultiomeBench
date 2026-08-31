#!/bin/bash
# Confirms multigrate-env works. Run on a GPU node. "ALL CHECKS PASSED" == success.
set -euo pipefail

export PYTHONNOUSERSITE=1
module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multigrate-env

echo "===== imports / API / CUDA ====="
python - <<'PY'
import sys
import torch, anndata, scanpy, scvi
import multigrate as mtg
print("python     :", sys.version.split()[0])
print("multigrate :", mtg.__version__, "from", mtg.__file__)
print("scvi-tools :", scvi.__version__)
print("torch      :", torch.__version__)

assert "/envs/multigrate-env/" in mtg.__file__, mtg.__file__
# API surface used by run_multigrate.py
assert hasattr(mtg.data, "organize_multimodal_anndatas"), "organize_multimodal_anndatas missing"
assert hasattr(mtg.model, "MultiVAE"), "MultiVAE missing"
for m in ("setup_anndata", "train", "get_model_output"):
    assert hasattr(mtg.model.MultiVAE, m), f"MultiVAE.{m} missing"
print("mtg.data.organize_multimodal_anndatas + mtg.model.MultiVAE API: OK")

ok = torch.cuda.is_available()
print("torch.cuda.is_available():", ok)
if ok:
    print("cuda device:", torch.cuda.get_device_name(0))
assert ok, "CUDA not visible -- are you on a gpu/dgx node? (Multigrate is a GPU VAE)"
print("\nALL CHECKS PASSED")
PY
