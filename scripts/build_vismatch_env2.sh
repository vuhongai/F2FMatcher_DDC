#!/usr/bin/env bash
# Rebuild the `vismatch` conda env CLEANLY: ensure the env has its OWN pip so
# packages land in the env site-packages (not the user ~/.local site).
set -x
CONDA=/home/avuhong/anaconda3/bin/conda
ENV=vismatch
PY=/home/avuhong/anaconda3/envs/$ENV/bin/python

$CONDA env remove -n $ENV -y 2>/dev/null || true
$CONDA create -n $ENV python=3.10 pip -y

# Install into the ENV site-packages (the env now has its own pip).
$PY -m pip install --no-cache-dir torch torchvision \
    --index-url https://download.pytorch.org/whl/cu130
$PY -m pip install --no-cache-dir vismatch dinov2 opencv-python-headless tqdm

echo "=== torch check (env site only) ==="
PYTHONNOUSERSITE=1 $PY -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
echo "=== vismatch check (env site only) ==="
PYTHONNOUSERSITE=1 $PY -c "from vismatch import get_matcher; print('vismatch OK')"
echo "ENV_BUILD2_DONE"
