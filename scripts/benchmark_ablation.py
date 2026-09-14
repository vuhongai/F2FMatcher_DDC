"""Score the F2FMatcher ablation: P/R/F1 vs GT at each pipeline step.

The ablation run (scripts/run_f2fmatcher_ablation.py) saved step_prediction.pkl
per pair, capturing the intermediate matched-label sets as the pipeline adds
steps. We score each stage against the curated GT (with the edge-ROI
correction: recall reported over ALL GT and over the matchable/non-edge GT) to
show the mapping quality improves as steps are added:

  S1 initial_guess                    : classifier seed pairs (pre-geometry)
  S2 geometry-validated seeds         : seeds passing the triangle-geometry check
  S3 + local propagation              : neighbour propagation (final iter)
  S4 + fill unannotated (FULL)        : affine fill of remaining ROIs

Usage:  python scripts/benchmark_ablation.py
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.measure import regionprops

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, str(REPO))
from f2fmatcher.segmentation.cellpose_seg import filter_ROIs  # noqa: E402

BASE = "/media/DATABRUT/DB_DDC/serverGPU/Cache_GPU_Ai/fiber_matcher/F2FMatcher"
ABLA = REPO / "results" / "benchmark" / "f2fmatcher_output_ablation"
PRED = ABLA / "prediction_output"
MASKS = ABLA / "out_CP_masks"          # symlinked to the main run's masks
GT_CSV = f"{BASE}/datasets/training_data_updated.csv"


def load_roi(img):
    with open(MASKS / f"{img}_CP_masks.pkl", "rb") as f:
        masks = pickle.load(f)[0]
    props = regionprops(masks)
    valid = set(filter_ROIs((masks, props, None)))
    is_edge = {int(r.label): (r.label not in valid) for r in props}
    return is_edge


def prf(tp, n_pred, n_gt):
    p = tp / n_pred if n_pred else 0.0
    r = tp / n_gt if n_gt else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f1


STAGES = ["S1 initial_guess", "S2 geometry seeds", "S3 +propagation", "S4 +fill (FULL)"]


def main():
    gt = pd.read_csv(GT_CSV)
    pairs = (gt.groupby(["Image1", "Image2"]).size().reset_index()[["Image1", "Image2"]]
             .values.tolist())
    # pooled accumulators per stage: [tp, n_pred, n_gt_all, n_gt_match]
    acc = {s: [0, 0, 0, 0] for s in STAGES}
    n_pairs = 0
    for img1, img2 in pairs:
        pdir = PRED / f"{img1}___vs___{img2}"
        sp_path = pdir / "step_prediction.pkl"
        pl_path = pdir / "paired_labels.pkl"
        if not sp_path.exists() or not pl_path.exists():
            continue
        n_pairs += 1
        gtset = set(zip(gt[gt.Image1 == img1].ROI_I1, gt[gt.Image2 == img2].ROI_I2))
        e1, e2 = load_roi(img1), load_roi(img2)
        gt_match = {(a, b) for a, b in gtset if not e1.get(a, True) and not e2.get(b, True)}
        sp = pickle.load(open(sp_path, "rb"))
        full = set(pickle.load(open(pl_path, "rb")))
        stages = {
            "S1 initial_guess": set(sp.get("1_initial_guess", [])),
            "S2 geometry seeds": set(sp.get("2_selected_combs_from_initial_guess", [])),
            "S3 +propagation": set(sp["3_local_prediction"][-1]) if sp.get("3_local_prediction") else set(full),
            "S4 +fill (FULL)": full,
        }
        for s in STAGES:
            pred = stages[s]
            tp = len(pred & gtset)
            acc[s][0] += tp
            acc[s][1] += len(pred)
            acc[s][2] += len(gtset)
            acc[s][3] += len(gt_match)

    print(f"\n=== F2FMatcher ablation: P/R/F1 vs GT at each step ({n_pairs} pairs) ===")
    print(f"{'stage':22s} {'#pred':>6s} {'P':>6s} {'R(all)':>7s} {'R(match)':>9s} {'F1(match)':>10s}")
    for s in STAGES:
        tp, np_, ng_all, ng_mat = acc[s]
        p = tp / np_ if np_ else 0.0
        r_all = tp / ng_all if ng_all else 0.0
        r_mat = tp / ng_mat if ng_mat else 0.0
        f1 = 2 * p * r_mat / (p + r_mat) if (p + r_mat) else 0.0
        print(f"{s:22s} {np_:6d} {p:6.3f} {r_all:7.3f} {r_mat:9.3f} {f1:10.3f}")
    print(f"\n(n_gt_all and n_gt_match are pooled over the {n_pairs} pairs with step predictions)")


if __name__ == "__main__":
    main()
