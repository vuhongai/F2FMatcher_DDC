#!/usr/bin/env bash
# Build the dedicated `vismatch` conda env for the DINOv2 + VisMatch baselines.
# Runs detached; logs to results/benchmark/vismatch_env.log
set -x
CONDA=/home/avuhong/anaconda3/bin/conda
PYVER=3.10
ENV=vismatch

$CONDA env remove -n $ENV -y 2>/dev/null || true
$CONDA create -n $ENV python=$PYVER -y

# CUDA torch first (GPU driver is CUDA 13.1; cu126 wheels are broadly compatible).
$CONDA run -n $ENV pip install --no-cache-dir torch torchvision \
    --index-url https://download.pytorch.org/whl/cu126

# VisMatch (unified 50+ matcher API) + DINOv2 + image I/O.
$CONDA run -n $ENV pip install --no-cache-dir vismatch dinov2 opencv-python-headless

echo "=== torch check ==="
$CONDA run -n $ENV python -c "import torch; print('torch', torch.__version__, 'cuda_avail', torch.cuda.is_available())"
echo "=== vismatch check ==="
$CONDA run -n $ENV python -c "import vismatch; print('vismatch OK', vismatch.__version__ if hasattr(vismatch,'__version__') else '')"
echo "ENV_BUILD_DONE"
