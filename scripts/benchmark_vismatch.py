"""Benchmark DINOv2 + VisMatch (dense/sparse matchers) against the curated GT.

Approach (pixel/patch -> ROI), per the user's spec:
  1. Run a matcher on the image pair -> pixel/patch correspondences.
  2. Assign each source pixel/patch to its ROI (CellPose label at that location).
  3. Aggregate to (ROI_0, ROI_1) correspondences.
  4. Score against the curated GT: precision / recall / F1 + coverage.

Methods (add/remove via --models):
  VisMatch (pip `vismatch`, unified 50+ matcher API, https://github.com/gmberton/vismatch):
    dense:      loftr, roma, tiny-roma
    semi-dense: eloftr, xoftr
    sparse:     superpoint-lightglue, superglue, sift-lightglue
  DINOv2 (torch.hub, patch-level embedding matching, ViT-S/14):
    dinov2

The matchers run on GPU (device='cuda'); ROI assignment + scoring are CPU.
Images are loaded at ORIGINAL size so matcher pixel coords live in the same
space as the CellPose label maps (verified: GT labels == annotation CellPose).

Usage:
  python scripts/benchmark_vismatch.py --models loftr,roma,superpoint-lightglue,dinov2
  python scripts/benchmark_vismatch.py --models dinov2 --device cpu   # quick CPU test
"""
import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

REPO = Path("/DATA/F2FMatcher_DDC")
BASE = "/media/DATABRUT/DB_DDC/serverGPU/Cache_GPU_Ai/fiber_matcher/F2FMatcher"
IMG = f"{BASE}/images"
GT_CSV = f"{BASE}/datasets/training_data_updated.csv"
MASKS = REPO / "results" / "benchmark" / "f2fmatcher_output" / "out_CP_masks"
OUT = REPO / "results" / "benchmark" / "vismatch"
OUT.mkdir(parents=True, exist_ok=True)

# DINOv2 normalization (ImageNet)
DINO_MEAN = np.array([0.485, 0.456, 0.406], dtype=np.float32)
DINO_STD = np.array([0.229, 0.224, 0.225], dtype=np.float32)
DINO_PATCH = 14


def load_mask(img):
    """CellPose label map (H, W) for an image, from the pipeline out_CP_masks."""
    with open(MASKS / f"{img}_CP_masks.pkl", "rb") as f:
        return pickle.load(f)[0]


def kpts_to_roi_pairs(kpts0, kpts1, mask0, mask1):
    """Map matched pixel correspondences to (ROI_0, ROI_1) pairs.

    kpts are (N, 2) in (x=col, y=row) in the ORIGINAL image space. A keypoint on
    background (label 0) is dropped. Returns a set of (roi0, roi1) with roi>0.
    """
    kpts0 = np.asarray(kpts0, dtype=np.float64)
    kpts1 = np.asarray(kpts1, dtype=np.float64)
    n = min(len(kpts0), len(kpts1))
    if n == 0:
        return set()
    kpts0, kpts1 = kpts0[:n], kpts1[:n]
    x0, y0 = kpts0[:, 0].astype(int), kpts0[:, 1].astype(int)
    x1, y1 = kpts1[:, 0].astype(int), kpts1[:, 1].astype(int)
    h0, w0 = mask0.shape
    h1, w1 = mask1.shape
    ok = ((x0 >= 0) & (x0 < w0) & (y0 >= 0) & (y0 < h0) &
          (x1 >= 0) & (x1 < w1) & (y1 >= 0) & (y1 < h1))
    x0, y0, x1, y1 = x0[ok], y0[ok], x1[ok], y1[ok]
    r0 = mask0[y0, x0]
    r1 = mask1[y1, x1]
    keep = (r0 > 0) & (r1 > 0)
    return set(zip(r0[keep].tolist(), r1[keep].tolist()))


def run_vismatch(model_name, img0, img1, mask0, mask1, device):
    from vismatch import get_matcher
    matcher = get_matcher(model_name, device=device)
    a = matcher.load_image(f"{IMG}/{img0}.png")   # original size
    b = matcher.load_image(f"{IMG}/{img1}.png")
    res = matcher(a, b)
    return kpts_to_roi_pairs(res["matched_kpts0"], res["matched_kpts1"], mask0, mask1), len(res["matched_kpts0"])


def _dinov2_patch_feats(model, image, size, device):
    """DINOv2 patch tokens for an image resized to `size` (multiple of 14).

    Returns (feats (N, D), centers (N, 2) in the RESIZED space, scale) where
    scale = original / resized (to map centers back to the original space).
    """
    import torch
    from PIL import Image as _Image
    orig = np.asarray(_Image.open(image).convert("RGB"))
    oh, ow = orig.shape[:2]
    im = _Image.open(image).convert("RGB").resize((size, size))
    x = np.asarray(im, dtype=np.float32) / 255.0
    x = (x - DINO_MEAN) / DINO_STD
    x = torch.from_numpy(x).permute(2, 0, 1).unsqueeze(0).to(device)
    with torch.no_grad():
        f = model.forward_features(x)["x_norm_patchtokens"]   # (1, N, D)
    f = f[0].cpu().numpy()
    n = size // DINO_PATCH
    gy, gx = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    centers = np.stack([gx.ravel() * DINO_PATCH + DINO_PATCH / 2,
                        gy.ravel() * DINO_PATCH + DINO_PATCH / 2], axis=1).astype(np.float64)
    scale = np.array([ow / size, oh / size])   # resized -> original
    return f, centers, scale


def run_dinov2(img0, img1, mask0, mask1, device, size=490, topk=1):
    import torch
    model = torch.hub.load("facebookresearch/dinov2", "dinov2_vits14", pretrained=True)
    model.eval().to(device)
    f0, c0, s0 = _dinov2_patch_feats(model, f"{IMG}/{img0}.png", size, device)
    f1, c1, s1 = _dinov2_patch_feats(model, f"{IMG}/{img1}.png", size, device)
    f0 = f0 / (np.linalg.norm(f0, axis=1, keepdims=True) + 1e-8)
    f1 = f1 / (np.linalg.norm(f1, axis=1, keepdims=True) + 1e-8)
    sim = f0 @ f1.T                                   # (N0, N1)
    pairs = set()
    n0 = len(f0)
    for i in range(n0):
        j = int(np.argmax(sim[i]))                    # best match in img1
        x0, y0 = (c0[i] * s0).astype(int)
        x1, y1 = (c1[j] * s1).astype(int)
        if 0 <= x0 < mask0.shape[1] and 0 <= y0 < mask0.shape[0] and \
           0 <= x1 < mask1.shape[1] and 0 <= y1 < mask1.shape[0]:
            r0, r1 = mask0[y0, x0], mask1[y1, x1]
            if r0 > 0 and r1 > 0:
                pairs.add((int(r0), int(r1)))
    return pairs, n0


def prf(tp, n_pred, n_gt):
    p = tp / n_pred if n_pred else 0.0
    r = tp / n_gt if n_gt else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="loftr,roma,superpoint-lightglue,dinov2")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--dinov2-size", type=int, default=490)
    args = ap.parse_args()
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    device = args.device

    gt = pd.read_csv(GT_CSV)
    pairs = (gt.groupby(["Image1", "Image2"]).size().reset_index()[["Image1", "Image2"]]
             .values.tolist())
    print(f"{len(pairs)} pairs, models={models}, device={device}", flush=True)

    rows = []
    for i, (img0, img1) in enumerate(pairs, 1):
        gtset = set(zip(gt[gt.Image1 == img0].ROI_I1, gt[gt.Image2 == img1].ROI_I2))
        mask0, mask1 = load_mask(img0), load_mask(img1)
        for m in models:
            try:
                if m == "dinov2":
                    pred, n_k = run_dinov2(img0, img1, mask0, mask1, device, args.dinov2_size)
                else:
                    pred, n_k = run_vismatch(m, img0, img1, mask0, mask1, device)
            except Exception as e:  # noqa: BLE001
                print(f"[{i}/{len(pairs)}] {m} ERROR {type(e).__name__}: {e}", flush=True)
                pred, n_k = set(), 0
            tp = len(pred & gtset)
            p, r, f1 = prf(tp, len(pred), len(gtset))
            # coverage: fraction of GT ROIs touched by >=1 correspondence
            cov0 = {a for a, _ in pred}
            cov1 = {b for _, b in pred}
            gt0 = {a for a, _ in gtset}
            gt1 = {b for _, b in gtset}
            rows.append({"Image1": img0, "Image2": img1, "method": m,
                         "n_gt": len(gtset), "n_pred": len(pred), "n_kpts": n_k,
                         "tp": tp, "precision": p, "recall": r, "f1": f1,
                         "cov0": len(cov0 & gt0) / len(gt0) if gt0 else 0.0,
                         "cov1": len(cov1 & gt1) / len(gt1) if gt1 else 0.0})
            print(f"[{i:2d}/{len(pairs)}] {m:22s} P={p:.3f} R={r:.3f} F1={f1:.3f} "
                  f"({len(pred)} pairs, {n_k} kpts)", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "vismatch_scores.csv", index=False)
    print("\n=== per-method pooled P/R/F1 + coverage ===")
    for m in models:
        sub = res[res.method == m]
        tp, np_, ng = sub.tp.sum(), sub.n_pred.sum(), sub.n_gt.sum()
        p, r, f1 = prf(tp, np_, ng)
        c0, c1 = sub.cov0.mean(), sub.cov1.mean()
        print(f"  {m:22s} pairs={len(sub):2d}  P={p:.3f} R={r:.3f} F1={f1:.3f}  "
              f"cov0={c0:.3f} cov1={c1:.3f}")
    print(f"saved {OUT/'vismatch_scores.csv'}")


if __name__ == "__main__":
    main()
