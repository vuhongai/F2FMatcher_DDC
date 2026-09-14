#!/usr/bin/env bash
# Install torch + vismatch + dinov2 into the vismatch ENV site-packages
# (PYTHONNOUSERSITE=1 so the user ~/.local site is ignored and the env's own
# pip is used). Logs to results/benchmark/vismatch_install.log
set -x
PY=/home/avuhong/anaconda3/envs/vismatch/bin/python
export PYTHONNOUSERSITE=1
$PY -m pip install --no-cache-dir torch torchvision \
    --index-url https://download.pytorch.org/whl/cu130
$PY -m pip install --no-cache-dir vismatch dinov2 opencv-python-headless tqdm
echo "=== verify (env site only) ==="
$PY -c "import torch; from vismatch import get_matcher; import dinov2; print('torch', torch.__version__, 'cuda', torch.cuda.is_available(), '| vismatch+dinov2 OK')"
echo "INSTALL_DONE"
