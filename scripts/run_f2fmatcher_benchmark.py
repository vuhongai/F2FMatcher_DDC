"""Run F2FMatcher on all 38 GT crop pairs (TA) for the §B benchmark.

Calls f2fmatcher's run_single_pair directly (avoids the CLI --param-file bug that
passes the .py param file to load_config->from_yaml). Uses the EXACT annotation
CellPose setup (CP_AV_WGA_Dia_Qua_TA_AxioScan10X, channels=[0,0], flow_threshold=0.4,
cellprob_threshold=0, diameter=model default) so the predicted labels live in the same
label space as the curated GT (verified: every GT label is reproduced).

Device: cpu (GPUs are occupied by the local LLM). Skips pairs that already have a
paired_labels.pkl (resume-safe).

Usage:  python scripts/run_f2fmatcher_benchmark.py
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, "/DATA/F2FMatcher/src")
from f2fmatcher.config import load_config  # noqa: E402
from f2fmatcher.scripts.run_pipeline import run_single_pair  # noqa: E402
import f2fmatcher.matching.spatial as _spatial  # noqa: E402


def _safe_angle(opposite, side1, side2):
    """NaN-safe triangle angle. The stock f2fmatcher angle() divides by
    2*side1*side2, which is 0 for a degenerate triangle (two coincident
    centroids) -> NaN -> downstream 'Input contains NaN' crash. Guard the
    denominator; a degenerate triangle gets angle 0 (cos=1). This is a
    crash-fix only; non-degenerate triangles are computed identically."""
    denom = 2 * side1 * side2
    with np.errstate(divide="ignore", invalid="ignore"):
        cos_val = (side1 ** 2 + side2 ** 2 - opposite ** 2) / denom
    cos_val = np.where(denom > 0, cos_val, 1.0)
    cos_val = np.clip(cos_val, -1.0, 1.0)
    return np.degrees(np.arccos(cos_val))


# Monkey-patch (triangle_geometry calls angle via module-global lookup).
_spatial.angle = _safe_angle

BASE = "/media/DATABRUT/DB_DDC/serverGPU/Cache_GPU_Ai/fiber_matcher/F2FMatcher"
IMG = f"{BASE}/images"
OUT = str(REPO / "results" / "benchmark" / "f2fmatcher_output")
CP = "CP_AV_WGA_Dia_Qua_TA_AxioScan10X"
N_PROC = 16


def main():
    device = torch.device("cpu")
    config = load_config()
    gt = pd.read_csv(f"{BASE}/datasets/training_data_updated.csv")
    pairs = (gt.groupby(["Image1", "Image2"]).size().reset_index()[["Image1", "Image2"]]
             .values.tolist())
    print(f"{len(pairs)} pairs, device={device}, out={OUT}", flush=True)
    t0 = time.time()
    for i, (img1, img2) in enumerate(pairs, 1):
        t1 = time.time()
        print(f"[{i:2d}/{len(pairs)}] {img1}  vs  {img2}", flush=True)
        try:
            n = run_single_pair(
                img1, img2, IMG, IMG, False, False, 0, 0,
                CP, CP, "fluorescence", "fluorescence", "10X", "10X",
                OUT, False, False, N_PROC, device,
                False, False, config, cellpose_gpu=False,
            )
            print(f"      -> {n} pairs  ({time.time()-t1:.0f}s)", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"      ERROR: {type(e).__name__}: {e}  ({time.time()-t1:.0f}s)", flush=True)
    print(f"\nAll done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
