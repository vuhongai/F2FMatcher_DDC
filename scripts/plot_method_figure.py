"""Real-data method-input figure (companion to fig1b_worked_example).

For one randomly chosen matched fibre pair (Laminin anchor vs NADH target,
WT QUAG01), draws the exact 256x256 model input for both sides:
image crop, target-ROI binary mask, Cellpose flow_x, flow_y, magnitude, angle.

Uses only precomputed artifacts (CP mask/flow pkls + paired_labels.pkl).

Outputs: manuscript/figures/fig1b2_real_input_pair.png (+ .pdf)
"""
import argparse
import os
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from PIL import Image
from skimage.measure import regionprops

warnings.filterwarnings("ignore", category=UserWarning, module="PIL")
Image.MAX_IMAGE_PIXELS = None

PROJECT_ROOT = Path(__file__).resolve().parent.parent
QUA_DIR = PROJECT_ROOT / "results" / "QUA"

IMG1 = "22-082_10X_DAPI_LAM_DYS_COL4_5-Scene-1-QUAG01"   # Laminin (IF), slide 1
IMG2 = "22-082_10X_NADH_13-Scene-1-QUAG01"               # NADH (brightfield)
CROP = 256


def load_side(img_name):
    with open(f"{QUA_DIR}/out_CP_masks/{img_name}_CP_masks.pkl", "rb") as f:
        masks = pickle.load(f)[0]
    with open(f"{QUA_DIR}/out_CP_masks/{img_name}_CP_flows.pkl", "rb") as f:
        flows = pickle.load(f)[0]
    flow_x = flows[1][0]
    flow_y = flows[1][1]
    image = np.array(Image.open(f"{QUA_DIR}/images_segmentation/{img_name}.png"))
    props = {p.label: p for p in regionprops(masks)}
    return {"masks": masks, "flow_x": flow_x, "flow_y": flow_y,
            "image": image, "props": props}


def valid_labels(masks, props):
    """Labels whose 256px window (filter_ROIs rule) fits in the image, area >= 100."""
    h, w = masks.shape
    half = CROP // 2
    out = []
    for label, p in props.items():
        cy, cx = p.centroid
        if (int(cy - half) >= 0 and int(cy - half) + CROP <= h
                and int(cx - half) >= 0 and int(cx - half) + CROP <= w
                and p.area >= 100):
            out.append(label)
    return out


def crop_side(side, label):
    cy, cx = side["props"][label].centroid
    y0, x0 = int(cy - CROP // 2), int(cx - CROP // 2)
    sl = (slice(y0, y0 + CROP), slice(x0, x0 + CROP))
    return {
        "image": side["image"][sl],
        "mask": (side["masks"][sl] == label).astype(np.uint8),
        "flow_x": side["flow_x"][sl],
        "flow_y": side["flow_y"][sl],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--label1", type=int, default=None, help="force anchor ROI label (Laminin)")
    ap.add_argument("--out", type=str, default=str(PROJECT_ROOT / "manuscript" / "figures" / "fig1b2_real_input_pair"))
    args = ap.parse_args()

    pair_dir = QUA_DIR / "prediction_output" / f"{IMG1}___vs___{IMG2}"
    with open(f"{pair_dir}/paired_labels.pkl", "rb") as f:
        pairs = pickle.load(f)
    pairs = np.array(pairs)

    print("Loading Laminin side ...", flush=True)
    s1 = load_side(IMG1)
    print("Loading NADH side ...", flush=True)
    s2 = load_side(IMG2)

    v1, v2 = set(valid_labels(s1["masks"], s1["props"])), set(valid_labels(s2["masks"], s2["props"]))
    ok = pairs[np.isin(pairs[:, 0], list(v1)) & np.isin(pairs[:, 1], list(v2))]
    print(f"matched pairs: {len(pairs)}, with valid 256px crops on both sides: {len(ok)}")

    if args.label1 is not None:
        sel = ok[ok[:, 0] == args.label1]
        if len(sel) == 0:
            sys.exit(f"label {args.label1} not in the matched set")
        l1, l2 = sel[0]
    else:
        rng = np.random.default_rng(args.seed)
        l1, l2 = ok[rng.integers(len(ok))]
    print(f"selected pair: Laminin ROI {l1} -> NADH ROI {l2} (seed={args.label1 if args.label1 is not None else args.seed})")

    c1, c2 = crop_side(s1, int(l1)), crop_side(s2, int(l2))

    fx = np.concatenate([c1["flow_x"].ravel(), c2["flow_x"].ravel()])
    fy = np.concatenate([c1["flow_y"].ravel(), c2["flow_y"].ravel()])
    lim = float(np.percentile(np.abs(np.concatenate([fx, fy])), 99.9))
    mag = np.hypot(c1["flow_x"], c1["flow_y"])
    mag2 = np.hypot(c2["flow_x"], c2["flow_y"])
    mag_max = max(float(mag.max()), float(mag2.max()))
    ang = np.arctan2(c1["flow_y"], c1["flow_x"])
    ang2 = np.arctan2(c2["flow_y"], c2["flow_x"])

    cols = ["256x256 crop", "mask (target ROI)", "flow_x", "flow_y", "magnitude", "angle"]
    rows = [f"Laminin (IF) - anchor ROI {l1}", f"NADH (brightfield) - matched ROI {l2}"]

    fig, axes = plt.subplots(2, 6, figsize=(18, 7.2), dpi=150)
    fig.subplots_adjust(hspace=0.42, wspace=0.18, left=0.05, right=0.985, top=0.88, bottom=0.07)

    for r, (c, row_title) in enumerate([(c1, rows[0]), (c2, rows[1])]):
        ax = axes[r, 0]
        if c["image"].ndim == 2:
            ax.imshow(c["image"], cmap="gray")
        else:
            ax.imshow(c["image"])
        ax.set_title(cols[0], fontsize=11)

        ax = axes[r, 1]
        ax.imshow(c["mask"], cmap="gray")
        ax.set_title(cols[1], fontsize=11)

        ax = axes[r, 2]
        ax.imshow(c["flow_x"], cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax.set_title(cols[2], fontsize=11)

        ax = axes[r, 3]
        ax.imshow(c["flow_y"], cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax.set_title(cols[3], fontsize=11)

        ax = axes[r, 4]
        m = mag if r == 0 else mag2
        ax.imshow(m, cmap="magma", vmin=0, vmax=mag_max)
        ax.set_title(cols[4], fontsize=11)

        ax = axes[r, 5]
        a = ang if r == 0 else ang2
        ax.imshow(a, cmap="twilight", vmin=-np.pi, vmax=np.pi)
        ax.set_title(cols[5], fontsize=11)

        for cc in range(6):
            axes[r, cc].set_xticks([])
            axes[r, cc].set_yticks([])
        axes[r, 0].set_ylabel(row_title, fontsize=12, fontweight="bold", rotation=0, ha="right", va="center", labelpad=80)

    for cc, title in enumerate(cols):
        axes[0, cc].set_title(title, fontsize=11, pad=8)

    arrow = FancyArrowPatch((0.115, 0.56), (0.115, 0.36), arrowstyle="-|>",
                            mutation_scale=22, lw=2.2, color="tab:blue",
                            shrinkA=0, shrinkB=0, figure=fig)
    fig.add_artist(arrow)
    fig.text(0.128, 0.46, "matched pair", color="tab:blue", fontsize=10, rotation=90, va="center")

    fig.suptitle("F2FMatcher input for one matched fibre pair (WT QUAG01): "
                 "256x256 Cellpose flow-field crops, target-ROI mask and flow components",
                 fontsize=13, y=0.965)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"saved {out.with_suffix('.png')} and {out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
