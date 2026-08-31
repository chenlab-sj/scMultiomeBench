#!/bin/bash
# Confirms scbutterfly-env works. Run on a GPU node. "ALL CHECKS PASSED" == success.
set -euo pipefail

export PYTHONNOUSERSITE=1
module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scbutterfly-env

echo "===== imports / API / CUDA ====="
python - <<'PY'
import sys, torch, scanpy, anndata
import scButterfly
from scButterfly.butterfly import Butterfly
print("python     :", sys.version.split()[0])
print("scButterfly:", getattr(scButterfly, "__version__", "?"), "from", scButterfly.__file__)
print("torch      :", torch.__version__)
# NOTE: scvi-tools 0.19.0 is installed (scButterfly's hard pin) but its chex/jax stack is
# broken on this cluster. scButterfly imports scvi LAZILY (only for MultiVI_augmentation,
# which we do NOT use -> aug_type=None), so a working scvi is not required. Don't import scvi.

assert "/envs/scbutterfly-env/" in scButterfly.__file__, scButterfly.__file__
b = Butterfly()
for m in ("load_data", "data_preprocessing", "construct_model", "train_model", "test_model"):
    assert hasattr(b, m), f"Butterfly.{m} missing"
print("Butterfly API: OK")

ok = torch.cuda.is_available()
print("torch.cuda.is_available():", ok)
if ok: print("cuda device:", torch.cuda.get_device_name(0))
assert ok, "CUDA not visible -- are you on a gpu/dgx node?"
print("\nALL CHECKS PASSED")
PY
