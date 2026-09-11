import os, sys, pickle
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from skimage.measure import regionprops

sys.path.append("/DATA/F2FMatcher_DDC")
from config.ddc_config import *

from f2fmatcher.segmentation.cellpose_seg import filter_ROIs

NEW_PANEL = "10X_WGA_Myh7_Myh2_Myh4"


def count_rois(img, dir_cp_masks):
    """Number of ROIs of an image in the same ROI set used by the mapping."""
    path_masks = Path(dir_cp_masks) / f"{img}_CP_masks.pkl"
    with open(path_masks, "rb") as f:
        masks = pickle.load(f)[0]
    return len(filter_ROIs((masks, regionprops(masks), None)))


def panel_of(img):
    # 22-082_10X_WGA_Myh7_Myh2_Myh4_1-Scene-1-QUAG01 -> 10X_WGA_Myh7_Myh2_Myh4
    parts = img.split("_")
    return "_".join(parts[1:-1])


def main():
    dir_pair_output = Path(PAIR_DIRS_BASE_QUA)
    dir_cp_masks = Path(CP_MASKS_DIR_QUA)

    pair_dirs = sorted(d for d in os.listdir(dir_pair_output) if "___vs___" in d)
    pairs = []
    for d in pair_dirs:
        img1, img2 = d.split("___vs___")
        pairs.append((img1, img2, d))
    print(f"Found {len(pairs)} pair directories in {dir_pair_output}")

    all_imgs = sorted({img for p in pairs for img in (p[0], p[1])})
    print(f"Counting ROIs for {len(all_imgs)} image(s) ...")
    roi_counts = Parallel(n_jobs=n_processes, return_as="generator")(
        delayed(count_rois)(img, dir_cp_masks) for img in all_imgs
    )
    rois = dict(zip(all_imgs, roi_counts))

    rows = []
    for img1, img2, d in pairs:
        n1, n2 = rois[img1], rois[img2]
        path = dir_pair_output / d / "paired_labels.pkl"
        base = {
            "panel": panel_of(img2),
            "is_new": panel_of(img2) == NEW_PANEL,
            "img1": img1, "img2": img2,
            "n_rois1": n1, "n_rois2": n2,
        }
        if not path.exists():
            base.update(n_pairs=np.nan, n_matched1=np.nan, n_matched2=np.nan,
                        coverage1=np.nan, coverage2=np.nan)
            rows.append(base)
            continue
        with open(path, "rb") as f:
            pairs_labels = pickle.load(f)
        l1 = {p[0] for p in pairs_labels}
        l2 = {p[1] for p in pairs_labels}
        base.update(n_pairs=len(pairs_labels), n_matched1=len(l1), n_matched2=len(l2),
                    coverage1=len(l1) / n1, coverage2=len(l2) / n2)
        rows.append(base)

    summary = pd.DataFrame(rows)
    out = Path(MAPPING_OUTPUT_DIR_QUA) / "mapping_summary.csv"
    summary.to_csv(out, index=False)
    print(f"\nSaved per-pair summary to {out}\n")

    ok = summary.dropna(subset=["n_pairs"])
    for (panel, is_new), g in ok.groupby(["panel", "is_new"]):
        cov = (g["coverage1"] + g["coverage2"]) / 2
        tag = "NEW" if is_new else "existing"
        print(f"Panel {panel} ({tag}): {len(g)} pairs")
        print(f"  rois ref: mean {g['n_rois1'].mean():.0f} | rois other: mean {g['n_rois2'].mean():.0f}")
        print(f"  pairs/pair: mean {g['n_pairs'].mean():.0f}, min {g['n_pairs'].min():.0f}, max {g['n_pairs'].max():.0f}")
        print(f"  coverage ref:  mean {g['coverage1'].mean():.1%}, min {g['coverage1'].min():.1%}, max {g['coverage1'].max():.1%}")
        print(f"  coverage other: mean {g['coverage2'].mean():.1%}, min {g['coverage2'].min():.1%}, max {g['coverage2'].max():.1%}")

    for is_new, g in ok.groupby("is_new"):
        cov = (g["coverage1"] + g["coverage2"]) / 2
        tag = "NEW (Myh4)" if is_new else "existing"
        print(f"\n{tag}: {len(g)} pairs, {g['n_pairs'].sum():.0f} matched pairs total, "
              f"mean coverage (both sides) {cov.mean():.1%}")

    missing = summary[summary["n_pairs"].isna()]
    if len(missing):
        print(f"\n{len(missing)} pair dir(s) WITHOUT paired_labels.pkl (incomplete/missing):")
        for panel, g in missing.groupby("panel"):
            print(f"  {panel}: {len(g)}")


if __name__ == "__main__":
    main()
