"""Re-run F2FMatcher on all 38 GT crop pairs with save_step_prediction=True.

Purpose: capture the intermediate matching states at each pipeline step so we can
show the mapping quality (P/R/F1 vs GT) improves as steps are added:
  1_initial_guess                      -> classifier seeds (before geometry)
  2_selected_combs_from_initial_guess  -> seeds that pass geometry validation
  3_local_prediction (list, per iter)  -> after neighbour propagation
  4_unannotated_prediction (list)      -> after affine fill of unannotated ROIs

The ablation output dir reuses the existing CellPose masks + segmentation images
(symlinked) so only the VAE embeddings + matching are recomputed (CPU; GPUs are
occupied by the local LLM). The stock angle() NaN crash is guarded with the same
_safe_angle monkey-patch as the main benchmark run (crash-fix only).

Usage:  python scripts/run_f2fmatcher_ablation.py
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
    denom = 2 * side1 * side2
    with np.errstate(divide="ignore", invalid="ignore"):
        cos_val = (side1 ** 2 + side2 ** 2 - opposite ** 2) / denom
    cos_val = np.where(denom > 0, cos_val, 1.0)
    cos_val = np.clip(cos_val, -1.0, 1.0)
    return np.degrees(np.arccos(cos_val))


_spatial.angle = _safe_angle

BASE = "/media/DATABRUT/DB_DDC/serverGPU/Cache_GPU_Ai/fiber_matcher/F2FMatcher"
IMG = f"{BASE}/images"
OUT = str(REPO / "results" / "benchmark" / "f2fmatcher_output_ablation")
CP = "CP_AV_WGA_Dia_Qua_TA_AxioScan10X"
N_PROC = 16


def main():
    device = torch.device("cpu")
    config = load_config()
    gt = pd.read_csv(f"{BASE}/datasets/training_data_updated.csv")
    pairs = (gt.groupby(["Image1", "Image2"]).size().reset_index()[["Image1", "Image2"]]
             .values.tolist())
    print(f"{len(pairs)} pairs, device={device}, out={OUT}, save_step_prediction=True",
          flush=True)
    t0 = time.time()
    for i, (img1, img2) in enumerate(pairs, 1):
        t1 = time.time()
        print(f"[{i:2d}/{len(pairs)}] {img1}  vs  {img2}", flush=True)
        try:
            n = run_single_pair(
                img1, img2, IMG, IMG, False, False, 0, 0,
                CP, CP, "fluorescence", "fluorescence", "10X", "10X",
                OUT, False, True, N_PROC, device,   # export_images=False, save_step_prediction=True
                False, False, config, cellpose_gpu=False,
            )
            print(f"      -> {n} pairs  ({time.time()-t1:.0f}s)", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"      ERROR: {type(e).__name__}: {e}  ({time.time()-t1:.0f}s)", flush=True)
    print(f"\nAll done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
