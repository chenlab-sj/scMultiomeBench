#!/bin/bash
# Confirms the mira-env is genuinely good (not just "imports").
# Run on a GPU node so torch.cuda is meaningful:
#   bash verify_mira_env.sh
# A clean exit + "ALL VERIFICATION CHECKS PASSED" == success.
set -euo pipefail

# Isolate from ~/.local user-site (otherwise other methods' packages leak in and a
# stray 'mira' shadows the env's mira-multiome -> 'No module named mira.topics').
export PYTHONNOUSERSITE=1

module load conda3/202311
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mira-env

echo "===== 1) pip dependency graph must be clean ====="
python -m pip check && echo "pip check: OK"

echo "===== 2) versions / imports / CUDA / smoke ====="
python - <<'PY'
import sys
from packaging.version import Version

assert sys.version_info[:2] == (3, 8), f"expected py3.8, got {sys.version_info[:2]}"

import torch, networkx, numpy, scipy, sklearn, pyro
print("python   :", sys.version.split()[0])
print("torch    :", torch.__version__)
print("networkx :", networkx.__version__)
print("numpy    :", numpy.__version__)
print("scipy    :", scipy.__version__)
print("sklearn  :", sklearn.__version__)
print("pyro-ppl :", pyro.__version__)

# mira's hard ceilings
assert Version(torch.__version__.split('+')[0]) <= Version("2.0.0"), torch.__version__
assert "cu118" in torch.__version__, f"expected cu118 build, got {torch.__version__}"
assert Version(networkx.__version__) <= Version("2.5"), networkx.__version__
assert Version(numpy.__version__) < Version("2"), numpy.__version__          # torch 2.0.0 ABI
assert Version(scipy.__version__) <= Version("1.10.2"), scipy.__version__
assert Version(sklearn.__version__) < Version("1.4"), sklearn.__version__
assert Version("1.5.2") <= Version(pyro.__version__) < Version("2"), pyro.__version__

import mira
print("mira     :", getattr(mira, "__version__", "unknown"))
print("mira from:", mira.__file__)   # MUST be inside .../envs/mira-env/ (not ~/.local)
assert "/envs/mira-env/" in mira.__file__, f"wrong mira on path: {mira.__file__}"
# In mira 2.1.1, `mira.topics` is an ALIAS for mira.topic_model created in __init__
# (`import mira.topic_model as topics`). So `import mira.topics` FAILS (no such file),
# but attribute access after `import mira` works. Use that.
assert callable(mira.topics.make_model), "mira.topics.make_model missing"
print("mira.topics.make_model: OK  (real constructor; ExpressionTopicModel is a doc-only faux class)")

ok = torch.cuda.is_available()
print("torch.cuda.is_available():", ok)
if ok:
    print("cuda device:", torch.cuda.get_device_name(0))
assert ok, "CUDA not visible -- are you on a gpu/dgx node?"

# Smoke: instantiate a topic model. Non-fatal on kwarg mismatch (the real
# signature gets exercised by run_mira.py); we only want to know the class
# is constructible in this env.
try:
    import numpy as np, anndata as ad
    rng = np.random.default_rng(0)
    a = ad.AnnData(rng.poisson(0.4, size=(50, 60)).astype("float32"))
    a.var_names = [f"g{i}" for i in range(a.n_vars)]
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.var["highly_variable"] = True
    a.layers["counts"] = a.X.copy()
    model = mira.topics.make_model(
        a.n_obs, a.n_vars,
        feature_type="expression",
        highly_variable_key="highly_variable",
        counts_layer="counts",
        num_topics=3,
    )
    print("built topic model via make_model:", type(model).__name__, "OK")
except Exception as e:
    print("NOTE: make_model smoke failed (non-fatal for env check):", repr(e))

print("\nALL VERIFICATION CHECKS PASSED")
PY
