"""Precompute per-fiber morphology for the WT evaluation cohort.

Loads each WT image's CellPose mask (large pkl) ONCE, computes regionprops,
applies the same ROI filter used by the mapping (filter_ROIs), and caches a
compact per-fiber table. The §B evaluation then works purely from this cache
+ paired_labels.pkl (no re-loading of the big masks).

Scale-invariant shape features (robust to the global affine between serial
sections) are stored alongside area/centroid:
    area, perimeter, eccentricity, solidity, extent,
    roundness = 4*pi*area/perimeter^2, aspect_ratio = major/minor axis,
    cx, cy (centroid in pixel coords)

Output: results/QUA/eval/fiber_morphology.pkl
    dict { img_name: DataFrame indexed by label_id }

Usage:  python scripts/precompute_fiber_morphology.py [--force]
"""
import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.measure import regionprops

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, str(REPO))
from config.ddc_config import CP_MASKS_DIR_QUA, PAIR_DIRS_BASE_QUA, threshold_fiber_area  # noqa: E402
from f2fmatcher.segmentation.cellpose_seg import filter_ROIs  # noqa: E402

WT_SAMPLES = ["QUAG01", "QUAG02", "QUAG03", "QUAG04", "QUAG05"]
OUT = REPO / "results" / "QUA" / "eval" / "fiber_morphology.pkl"


def wt_pair_dirs():
    """All WT pair dirs that have a paired_labels.pkl (skip 5X_HE)."""
    dirs = []
    for d in sorted(PAIR_DIRS_BASE_QUA.iterdir()):
        if not d.is_dir() or "___vs___" not in d.name:
            continue
        if not any(s in d.name for s in WT_SAMPLES):
            continue
        if (d / "paired_labels.pkl").exists():
            dirs.append(d)
    return dirs


def img_features(img):
    path = CP_MASKS_DIR_QUA / f"{img}_CP_masks.pkl"
    with open(path, "rb") as f:
        masks = pickle.load(f)[0]
    props = regionprops(masks)
    valid = set(filter_ROIs((masks, props, None)))
    rows = []
    for r in props:
        if r.label not in valid:
            continue
        perim = r.perimeter
        roundness = (4 * np.pi * r.area) / (perim ** 2) if perim > 0 else 0.0
        aspect = (r.major_axis_length / r.minor_axis_length) if r.minor_axis_length > 0 else 0.0
        cy, cx = r.centroid
        rows.append({
            "label": int(r.label),
            "area": float(r.area),
            "perimeter": float(perim),
            "eccentricity": float(r.eccentricity),
            "solidity": float(r.solidity),
            "extent": float(r.extent),
            "roundness": float(roundness),
            "aspect_ratio": float(aspect),
            "cx": float(cx),
            "cy": float(cy),
        })
    del masks, props
    df = pd.DataFrame(rows).set_index("label")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="recompute even if cache exists")
    args = ap.parse_args()

    if OUT.exists() and not args.force:
        print(f"cache exists: {OUT} (use --force to recompute)")
        return

    pair_dirs = wt_pair_dirs()
    imgs = []
    for d in pair_dirs:
        img1, img2 = d.name.split("___vs___")
        for img in (img1, img2):
            if img not in imgs:
                imgs.append(img)
    print(f"WT pairs: {len(pair_dirs)}  |  unique images: {len(imgs)}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cache = {}
    t0 = time.time()
    for i, img in enumerate(imgs, 1):
        t1 = time.time()
        df = img_features(img)
        cache[img] = df
        print(f"[{i:2d}/{len(imgs)}] {img[-40:]:40s} n={len(df):5d}  ({time.time()-t1:5.1f}s)")
        del df
    with open(OUT, "wb") as f:
        pickle.dump(cache, f)
    n_tot = sum(len(d) for d in cache.values())
    print(f"\nSaved {len(cache)} images, {n_tot:,} fibers -> {OUT}  ({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()
