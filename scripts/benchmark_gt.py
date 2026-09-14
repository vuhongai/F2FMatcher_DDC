"""Score F2FMatcher + geometric baselines against the curated GT (TA crop pairs).

For each of the 38 GT pairs, every method produces a set of (label1, label2)
correspondences in the SAME CellPose label space as the GT (verified: every GT
label is reproduced by the annotation CellPose setup). A predicted pair is a true
positive iff it is in the curated GT. We report precision / recall / F1.

Methods:
  - F2FMatcher        : prediction_output/<i1>___vs___<i2>/paired_labels.pkl
  - no-align kNN      : nearest neighbour in raw centroid coordinates
  - affine+kNN        : global moment-based similarity (scale+rot+trans) + kNN
  - random            : null (same # pairs as F2FMatcher, random label pairing)

Usage:  python scripts/benchmark_gt.py
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, str(REPO))
from f2fmatcher.segmentation.cellpose_seg import filter_ROIs  # noqa: E402
from skimage.measure import regionprops  # noqa: E402

BASE = "/media/DATABRUT/DB_DDC/serverGPU/Cache_GPU_Ai/fiber_matcher/F2FMatcher"
OUT = REPO / "results" / "benchmark" / "f2fmatcher_output"
MASKS = OUT / "out_CP_masks"
PRED = OUT / "prediction_output"
GT_CSV = f"{BASE}/datasets/training_data_updated.csv"
RNG = np.random.default_rng(42)


def load_roi(img):
    """Return (cent_valid, is_edge) for an image from the pipeline masks.

    cent_valid: {label: centroid} for ROIs kept by filter_ROIs (non-edge, area>=100).
    is_edge:    {label: bool} True if the ROI is discarded by filter_ROIs (centroid
                within `half` px of any border, or area < 100). F2FMatcher can only
                match non-edge ROIs, so edge GT pairs are not matchable by design.
    """
    with open(MASKS / f"{img}_CP_masks.pkl", "rb") as f:
        masks = pickle.load(f)[0]
    props = regionprops(masks)
    valid = set(filter_ROIs((masks, props, None)))
    cent = {int(r.label): np.array(r.centroid) for r in props if r.label in valid}
    is_edge = {int(r.label): (r.label not in valid) for r in props}
    return cent, is_edge


def fit_similarity(P1, P2):
    """Correspondence-free global similarity (scale+rotation+translation) by
    aligning the two centroid clouds' second moments. Falls back to
    translation-only on degenerate clouds (NaN/zero-variance)."""
    c1, c2 = P1.mean(0), P2.mean(0)
    C1 = np.cov((P1 - c1).T)
    C2 = np.cov((P2 - c2).T)
    tr1, tr2 = np.trace(C1), np.trace(C2)
    if not (np.isfinite(tr1) and np.isfinite(tr2)) or tr1 < 1e-6 or tr2 < 1e-6:
        return 1.0, np.eye(2), c2 - c1
    s = float(np.sqrt(tr2 / tr1))
    _, V1 = np.linalg.eigh(C1)
    _, V2 = np.linalg.eigh(C2)
    R = V2 @ V1.T
    if not np.all(np.isfinite(R)):
        return 1.0, np.eye(2), c2 - c1
    if np.linalg.det(R) < 0:
        V2 = V2.copy()
        V2[:, -1] *= -1
        R = V2 @ V1.T
    t = c2 - s * (R @ c1)
    if not np.all(np.isfinite(t)):
        return 1.0, np.eye(2), c2 - c1
    return s, R, t


def prf(tp, n_pred, n_gt):
    p = tp / n_pred if n_pred else 0.0
    r = tp / n_gt if n_gt else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f1


def main():
    gt = pd.read_csv(GT_CSV)
    pairs = (gt.groupby(["Image1", "Image2"]).size().reset_index()[["Image1", "Image2"]]
             .values.tolist())
    rows = []
    for i1, i2 in pairs:
        gtset = set(zip(gt[gt.Image1 == i1].ROI_I1, gt[gt.Image2 == i2].ROI_I2))
        c1, e1 = load_roi(i1)
        c2, e2 = load_roi(i2)
        # a GT pair is matchable by F2FMatcher only if BOTH ROIs are non-edge
        gt_matchable = {(a, b) for a, b in gtset
                        if not e1.get(a, True) and not e2.get(b, True)}

        # --- F2FMatcher ---
        ppath = PRED / f"{i1}___vs___{i2}" / "paired_labels.pkl"
        f2f = set() if not ppath.exists() else set(
            pickle.load(open(ppath, "rb")))
        f2f = {(a, b) for a, b in f2f if a in c1 and b in c2}

        # --- no-align kNN ---
        P1 = np.array([c1[l] for l in c1]); l1 = list(c1)
        P2 = np.array([c2[l] for l in c2]); l2 = list(c2)
        tree2 = cKDTree(P2)
        nn_raw = tree2.query(P1, k=1)[1]
        raw = {(l1[i], l2[nn_raw[i]]) for i in range(len(P1))}

        # --- affine+kNN (moment-based similarity) ---
        s, R, t = fit_similarity(P1, P2)
        P1t = s * (P1 @ R.T) + t
        nn_aff = cKDTree(P2).query(P1t, k=1)[1]
        aff = {(l1[i], l2[nn_aff[i]]) for i in range(len(P1))}

        # --- random null (same count as F2FMatcher) ---
        n = len(f2f)
        rand = set()
        if n:
            a = RNG.choice(len(P1), n, replace=False)
            b = RNG.choice(len(P2), n, replace=False)
            rand = {(l1[x], l2[y]) for x, y in zip(a, b)}

        for name, pred in [("F2FMatcher", f2f), ("no-align kNN", raw),
                           ("affine+kNN", aff), ("random", rand)]:
            tp = len(pred & gtset)
            p, r, f1 = prf(tp, len(pred), len(gtset))
            r_match = tp / len(gt_matchable) if gt_matchable else 0.0
            rows.append({"Image1": i1, "Image2": i2, "method": name,
                         "n_gt": len(gtset), "n_gt_matchable": len(gt_matchable),
                         "n_pred": len(pred), "tp": tp,
                         "precision": p, "recall": r, "recall_matchable": r_match,
                         "f1": f1, "f2f_ran": ppath.exists()})
    res = pd.DataFrame(rows)
    res.to_csv(OUT.parent / "benchmark_scores.csv", index=False)

    # aggregate (micro / pooled) per method, F2FMatcher only over pairs it ran
    # Two recall denominators:
    #   R_all        : over ALL curated GT pairs (includes edge ROIs the method
    #                  cannot match by design -> underestimates true recall)
    #   R_matchable  : over GT pairs whose two ROIs are both non-edge (fair)
    print("\n=== per-method pooled P / R(all GT) / R(matchable) / F1 (all 38 pairs) ===")
    for m in ["F2FMatcher", "no-align kNN", "affine+kNN", "random"]:
        sub = res[res.method == m]
        if m == "F2FMatcher":
            sub = sub[sub.f2f_ran]
        tp, np_ = sub.tp.sum(), sub.n_pred.sum()
        ng_all, ng_mat = sub.n_gt.sum(), sub.n_gt_matchable.sum()
        p = tp / np_ if np_ else 0.0
        r_all = tp / ng_all if ng_all else 0.0
        r_mat = tp / ng_mat if ng_mat else 0.0
        f1 = 2 * p * r_mat / (p + r_mat) if (p + r_mat) else 0.0
        print(f"  {m:14s} pairs={len(sub):2d}  TP={tp:5d}  "
              f"P={p:.3f}  R_all={r_all:.3f}  R_match={r_mat:.3f}  F1(match)={f1:.3f}")
    n_edge_pairs = int((res[res.method == "F2FMatcher"].n_gt
                        - res[res.method == "F2FMatcher"].n_gt_matchable).sum())
    tot_gt = int(res[res.method == "F2FMatcher"].n_gt.sum())
    print(f"\nF2FMatcher ran on {int(res[res.method=='F2FMatcher'].f2f_ran.sum())}/38 pairs")
    print(f"GT pairs: {tot_gt} total, {tot_gt - n_edge_pairs} matchable (non-edge), "
          f"{n_edge_pairs} with >=1 edge ROI (unmatchable by F2FMatcher design)")
    print(f"saved {OUT.parent/'benchmark_scores.csv'}")


if __name__ == "__main__":
    main()
