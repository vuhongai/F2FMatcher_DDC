"""Per-step prediction-progression figures (real data).

Reads step_prediction.pkl (produced by run_step_prediction_single.py) and draws
one full-slide figure per image per step, highlighting the predicted fibre
pairs — same style as the original Visualize_step_prediction notebook:
  IHF (fluorescence):  red contours, white fill   (confirmed alpha 0.5 / on-going 0.2)
  brightfield:         black contours, red fill   (confirmed alpha 0.5 / on-going 0.2)

Steps:
  Step0   plain image (no highlights)
  Step1   1_initial_guess                              (on-going)
  Step2   2_selected_combs_from_initial_guess          (confirmed)
  Step3-i 3_local_prediction[i]                        (confirmed)
  Step4   4_unannotated_prediction[-1]                 (confirmed)

Outputs (manuscript/figures/step_prediction_real/):
  {img1}_Step0.png, {img1}_Step1.png, ... , {img2}_Step4.png
  side_by_side.gif  (img1 | img2 per step)
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
from matplotlib.patches import Polygon
from PIL import Image
from joblib import Parallel, delayed
from skimage.measure import regionprops, find_contours
from tqdm import tqdm

warnings.filterwarnings("ignore")
Image.MAX_IMAGE_PIXELS = None

PROJECT_ROOT = Path(__file__).resolve().parent.parent
QUA_DIR = PROJECT_ROOT / "results" / "QUA"

IMG1 = "22-082_10X_DAPI_LAM_DYS_COL4_5-Scene-1-QUAG01"   # Laminin, IHF
IMG2 = "22-082_10X_NADH_13-Scene-1-QUAG01"               # NADH, brightfield
IHF1, IHF2 = True, False

STYLE = {
    "IHF": {
        "confirmed": {"color": "red", "lw": 1.0, "face": "white", "alpha": 0.5},
        "on-going": {"color": "red", "lw": 0.5, "face": "white", "alpha": 0.2},
    },
    "BrightField": {
        "confirmed": {"color": "black", "lw": 1.0, "face": "red", "alpha": 0.5},
        "on-going": {"color": "black", "lw": 0.5, "face": "red", "alpha": 0.2},
    },
}


def roi_contours(args):
    """find_contours restricted to the ROI bounding box (fast), coords offset back.

    NOTE: regionprops bbox is (min_row, min_col, max_row, max_col).
    The crop is padded by 1 px: if the region touches the crop edge,
    find_contours returns OPEN fragments (truncated at the array boundary),
    and matplotlib 'closes' them with a straight chord -> partial fill.
    """
    label, masks, bbox = args
    y0, x0, y1, x1 = bbox
    H, W = masks.shape
    py0, px0 = max(y0 - 1, 0), max(x0 - 1, 0)
    py1, px1 = min(y1 + 1, H), min(x1 + 1, W)
    if (py1 - py0) < 2 or (px1 - px0) < 2:
        return label, []
    sub = (masks[py0:py1, px0:px1] == label)
    cs = find_contours(sub, 0.5)
    return label, [c + np.array([py0, px0]) for c in cs]


def precompute_contours(masks, n_jobs=32):
    props = regionprops(masks)
    bboxes = {p.label: p.bbox for p in props}
    args = [(lbl, masks, bb) for lbl, bb in bboxes.items()]
    out = {}
    for lbl, cs in Parallel(n_jobs=n_jobs)(
            delayed(roi_contours)(a) for a in tqdm(args, desc="contours")):
        out[lbl] = cs
    return out


def draw_step(image, contours, labels, style_key, roi_type, title, out_path, dpi=100):
    h, w = image.shape[:2]
    fig, ax = plt.subplots(figsize=(w / 300, h / 300), dpi=dpi)
    if image.ndim == 2:
        ax.imshow(image, cmap="gray")
    else:
        ax.imshow(image)
    st = STYLE[style_key][roi_type]
    for lbl in labels:
        for c in contours.get(lbl, []):
            xy = np.c_[c[:, 1], c[:, 0]]
            ax.add_patch(Polygon(xy, closed=True, facecolor=st["face"],
                                 edgecolor="none", alpha=st["alpha"]))
            ax.plot(xy[:, 0], xy[:, 1], color=st["color"], linewidth=st["lw"])
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.text(0.01, 0.985, title, transform=ax.transAxes, fontsize=14,
                color="white" if style_key == "IHF" else "black",
                va="top", ha="left", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="black" if style_key == "IHF" else "white",
                          alpha=0.6))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default=str(QUA_DIR))
    ap.add_argument("--img1", default=IMG1)
    ap.add_argument("--img2", default=IMG2)
    ap.add_argument("--out", default=str(PROJECT_ROOT / "manuscript" / "figures" / "step_prediction_real"))
    ap.add_argument("--dpi", type=int, default=100)
    ap.add_argument("--gif-height", type=int, default=768)
    ap.add_argument("--steps", default=None,
                    help="comma-separated suffixes to (re)draw, e.g. 'Step4' (gif is always built from all steps)")
    args = ap.parse_args()
    only = set(args.steps.split(",")) if args.steps else None

    data_root = Path(args.data_root)
    pair_dir = data_root / "prediction_output" / f"{args.img1}___vs___{args.img2}"
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(f"{pair_dir}/step_prediction.pkl", "rb") as f:
        sp = pickle.load(f)
    print("steps:", {k: (len(v) if isinstance(v, list) else 1) for k, v in sp.items()})

    imgs = [args.img1, args.img2]
    # IHF (fluorescence) images carry the DAPI_LAM prefix; NADH/HE/COX are brightfield
    ihs = ["DAPI_LAM" in i for i in imgs]

    images, masks, contours, style_keys = [], [], [], []
    for img, ihf in zip(imgs, ihs):
        image = np.array(Image.open(f"{data_root}/images_segmentation/{img}.png"))
        with open(f"{data_root}/out_CP_masks/{img}_CP_masks.pkl", "rb") as f:
            m = pickle.load(f)[0]
        print(f"precomputing contours for {img} ...", flush=True)
        c = precompute_contours(m)
        n_empty = sum(1 for v in c.values() if not v)
        print(f"  contours: {len(c)} regions, {n_empty} with empty contour "
              f"(expect only 1-px slivers)")
        images.append(image)
        masks.append(m)
        contours.append(c)
        style_keys.append("IHF" if ihf else "BrightField")

    # build the step list: (suffix, roi_type, pairs)
    steps = [("Step0", "on-going", [])]
    steps.append(("Step1", "on-going", sp["1_initial_guess"]))
    steps.append(("Step2", "confirmed", sp["2_selected_combs_from_initial_guess"]))
    for i, pl in enumerate(sp["3_local_prediction"]):
        steps.append((f"Step3-{i}", "confirmed", pl))
    # Step4 = the ACTUAL final mapping (paired_labels.pkl). The last
    # 4_unannotated_prediction snapshot is taken BEFORE the final re-validation
    # pass, so it can contain pairs that were later removed.
    final_path = pair_dir / "paired_labels.pkl"
    if final_path.exists():
        with open(final_path, "rb") as f:
            final_pairs = pickle.load(f)
    else:
        unannot = sp.get("4_unannotated_prediction") or []
        final_pairs = unannot[-1] if unannot else sp["3_local_prediction"][-1]
    steps.append(("Step4", "confirmed", final_pairs))

    step_titles = {
        "Step0": "Step 0 - start (no pairs)",
        "Step1": f"Step 1 - initial guess ({len(sp['1_initial_guess'])} candidate pairs)",
        "Step2": f"Step 2 - triangle-consistent seeds ({len(sp['2_selected_combs_from_initial_guess'])} pairs)",
        "Step4": f"Step 4 - final after affine fill ({len(final_pairs)} pairs)",
    }

    for suffix, roi_type, pairs in steps:
        if only is not None and suffix not in only:
            continue
        n = len(pairs)
        for img_index, img in enumerate(imgs):
            labels = [p[img_index] for p in pairs]
            title = step_titles.get(suffix)
            if suffix.startswith("Step3-"):
                title = f"Step 3.{suffix.split('-')[1]} - local propagation ({n} pairs)"
            out_path = out_dir / f"{img}_{suffix}.png"
            draw_step(images[img_index], contours[img_index], labels,
                      style_keys[img_index], roi_type, title, out_path, dpi=args.dpi)
            print(f"saved {out_path.name} ({n} pairs)", flush=True)

    # side-by-side gif
    print("building side-by-side gif ...", flush=True)
    from PIL import Image as PILImage
    frames = []
    for suffix, _, pairs in steps:
        row = []
        for img_index, img in enumerate(imgs):
            p = out_dir / f"{img}_{suffix}.png"
            im = PILImage.open(p)
            scale = args.gif_height / im.height
            row.append(im.resize((int(im.width * scale), args.gif_height), PILImage.BILINEAR))
        wtot = sum(r.width for r in row) + 8
        canvas = PILImage.new("RGB", (wtot, args.gif_height), (30, 30, 30))
        x = 0
        for r in row:
            canvas.paste(r, (x, 0))
            x += r.width + 8
        frames.append(canvas)
    gif_path = out_dir / "side_by_side.gif"
    frames[0].save(gif_path, save_all=True, append_images=frames[1:],
                   duration=1200, loop=0, optimize=True)
    print(f"saved {gif_path} ({gif_path.stat().st_size/1e6:.1f} MB, {len(frames)} frames)")


if __name__ == "__main__":
    main()
