#!/bin/bash
# Confirms maxfuse-env is genuinely good. Clean exit + "ALL CHECKS PASSED" == success.
#   bash verify_maxfuse_env.sh
set -euo pipefail

export PYTHONNOUSERSITE=1
module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate maxfuse-env

echo "===== 1) pip dependency graph =====";  python -m pip check && echo "pip check: OK"

echo "===== 2) import + Fusor API + tiny smoke ====="
python - <<'PY'
import sys
from packaging.version import Version
assert sys.version_info[:2] == (3, 8), sys.version_info[:2]

import numpy, scipy, sklearn, scanpy
print("python :", sys.version.split()[0])
print("numpy  :", numpy.__version__)
print("scipy  :", scipy.__version__)
print("sklearn:", sklearn.__version__)
print("scanpy :", scanpy.__version__)
assert Version(numpy.__version__) < Version("1.24"), numpy.__version__  # np.int/np.float aliases

import maxfuse as mf
print("maxfuse from:", mf.__file__)
assert "/envs/maxfuse-env/" in mf.__file__, mf.__file__
assert hasattr(mf, "model") and hasattr(mf.model, "Fusor"), "mf.model.Fusor missing"
print("mf.model.Fusor: OK")

# tiny smoke: build a Fusor on dummy aligned-feature data and run the first steps
rng = numpy.random.default_rng(0)
n1, n2, shared, act = 120, 100, 40, 50
s1 = rng.poisson(0.5, size=(n1, shared)).astype("float32")
s2 = rng.poisson(0.5, size=(n2, shared)).astype("float32")   # same #shared features
a1 = rng.poisson(0.5, size=(n1, act)).astype("float32")
a2 = rng.poisson(0.5, size=(n2, act)).astype("float32")
fusor = mf.model.Fusor(shared_arr1=s1, shared_arr2=s2,
                       active_arr1=a1, active_arr2=a2, labels1=None, labels2=None)
fusor.split_into_batches(max_outward_size=5000, matching_ratio=3, metacell_size=2, verbose=False)
print("constructed Fusor + split_into_batches: OK")
print("\nALL CHECKS PASSED")
PY
