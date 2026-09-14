"""Whole-slide baseline matching for the §B coverage-vs-correctness comparison (WT).

For each WT whole-slide section-pair (anchor laminin -> panel), run each learned
matcher on the FULL sections and reduce the pixel/patch correspondences to a
**1-to-1 per-fibre assignment** {anchor_label -> panel_label} by majority vote of
keypoint hits (the anchor ROI is assigned to the panel ROI it hits most). We save
that assignment and the wall-clock runtime per pair. The geometric baselines and
F2FMatcher itself are scored separately by `wholeslide_score_and_plot.py` (they
need no GPU / no images), so THIS script only runs the appearance/dense matchers.

Whole slides are large; generic matchers cannot run at native resolution, so each
side is downscaled so its long edge <= --max-size (keypoints are scaled back into
the ORIGINAL CellPose-mask space before the pixel->ROI vote). The chosen downscale
and the runtime are recorded — they are part of the scalability story, not hidden.

Outputs (under results/QUA/eval/wholeslide/):
    {method}/{img1}___vs___{img2}.pkl   {anchor_label:int -> panel_label:int}
    runtime.csv                          per (pair, method): seconds, n_kpts, max_size

Usage (on the GPU box, WT data present):
    python scripts/wholeslide_run_baselines.py --models dinov2,loftr,roma,superpoint-lightglue
    python scripts/wholeslide_run_baselines.py --models dinov2 --device cpu --limit-pairs 2   # smoke test
"""
import argparse
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None  # whole-slide PNGs are large

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, str(REPO))
from config.ddc_config import CP_MASKS_DIR_QUA, PAIR_DIRS_BASE_QUA, MAPPING_OUTPUT_DIR_QUA  # noqa: E402

WT_SAMPLES = ["QUAG01", "QUAG02", "QUAG03", "QUAG04", "QUAG05"]
OUT = REPO / "results" / "QUA" / "eval" / "wholeslide"
DINO_MEAN = np.array([0.485, 0.456, 0.406], np.float32)
DINO_STD = np.array([0.229, 0.224, 0.225], np.float32)
DINO_PATCH = 14


def wt_pairs():
    out = []
    for d in sorted(PAIR_DIRS_BASE_QUA.iterdir()):
        if not d.is_dir() or "___vs___" not in d.name:
            continue
        if not any(s in d.name for s in WT_SAMPLES):
            continue
        if not (d / "paired_labels.pkl").exists():
            continue
        img1, img2 = d.name.split("___vs___")
        out.append((img1, img2, d))
    return out


def load_mask(img):
    with open(CP_MASKS_DIR_QUA / f"{img}_CP_masks.pkl", "rb") as f:
        return pickle.load(f)[0]


def load_image_scaled(img_dir, img, max_size):
    """Return (RGB float array HxWx3 in [0,1] resized, scale=(orig/resized) x,y)."""
    im = Image.open(Path(img_dir) / f"{img}.png").convert("RGB")
    ow, oh = im.size
    s = max(ow, oh) / float(max_size)
    if s > 1:
        im = im.resize((int(round(ow / s)), int(round(oh / s))))
    rw, rh = im.size
    return np.asarray(im, np.float32) / 255.0, np.array([ow / rw, oh / rh])


def vote_assign(kpts0, kpts1, scale0, scale1, mask0, mask1):
    """Pixel correspondences (in resized space) -> 1:1 {anchor_label: panel_label}
    by majority vote. kpts are (N,2) (x,y); scale maps resized->original."""
    k0 = np.asarray(kpts0, float); k1 = np.asarray(kpts1, float)
    n = min(len(k0), len(k1))
    if n == 0:
        return {}
    x0 = (k0[:n, 0] * scale0[0]).astype(int); y0 = (k0[:n, 1] * scale0[1]).astype(int)
    x1 = (k1[:n, 0] * scale1[0]).astype(int); y1 = (k1[:n, 1] * scale1[1]).astype(int)
    h0, w0 = mask0.shape; h1, w1 = mask1.shape
    ok = ((x0 >= 0) & (x0 < w0) & (y0 >= 0) & (y0 < h0) &
          (x1 >= 0) & (x1 < w1) & (y1 >= 0) & (y1 < h1))
    r0 = mask0[y0[ok], x0[ok]]; r1 = mask1[y1[ok], x1[ok]]
    keep = (r0 > 0) & (r1 > 0)
    votes = defaultdict(lambda: defaultdict(int))
    for a, b in zip(r0[keep].tolist(), r1[keep].tolist()):
        votes[a][b] += 1
    return {int(a): int(max(bs, key=bs.get)) for a, bs in votes.items()}


def run_vismatch(model_name, img_dir, img0, img1, mask0, mask1, device, max_size):
    from vismatch import get_matcher
    matcher = get_matcher(model_name, device=device, max_num_keypoints=8192)
    a, s0 = _vismatch_load(matcher, img_dir, img0, max_size)
    b, s1 = _vismatch_load(matcher, img_dir, img1, max_size)
    res = matcher(a, b)
    k0, k1 = res["matched_kpts0"], res["matched_kpts1"]
    return vote_assign(k0, k1, s0, s1, mask0, mask1), int(len(k0))


def _vismatch_load(matcher, img_dir, img, max_size):
    """Load via the matcher's own loader when possible; else fall back to PIL.
    Returns (image_tensor, scale resized->original)."""
    from PIL import Image as _I
    p = str(Path(img_dir) / f"{img}.png")
    ow, oh = _I.open(p).size
    try:
        t = matcher.load_image(p, resize=max_size)  # most vismatch loaders accept resize
        # infer resized long edge from tensor shape (C,H,W)
        rh, rw = int(t.shape[-2]), int(t.shape[-1])
        return t, np.array([ow / rw, oh / rh])
    except TypeError:
        t = matcher.load_image(p)
        rh, rw = int(t.shape[-2]), int(t.shape[-1])
        return t, np.array([ow / rw, oh / rh])


def run_dinov2(img_dir, img0, img1, mask0, mask1, device, max_size):
    import torch
    size = (max_size // DINO_PATCH) * DINO_PATCH  # multiple of patch
    model = getattr(run_dinov2, "_model", None)
    if model is None:
        model = torch.hub.load("facebookresearch/dinov2", "dinov2_vits14", pretrained=True).eval().to(device)
        run_dinov2._model = model

    def feats(img):
        from PIL import Image as _I
        ow, oh = _I.open(Path(img_dir) / f"{img}.png").size
        im = _I.open(Path(img_dir) / f"{img}.png").convert("RGB").resize((size, size))
        x = (np.asarray(im, np.float32) / 255.0 - DINO_MEAN) / DINO_STD
        x = torch.from_numpy(x).permute(2, 0, 1).unsqueeze(0).to(device)
        with torch.no_grad():
            f = model.forward_features(x)["x_norm_patchtokens"][0].cpu().numpy()
        n = size // DINO_PATCH
        gy, gx = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        c = np.stack([gx.ravel() * DINO_PATCH + DINO_PATCH / 2,
                      gy.ravel() * DINO_PATCH + DINO_PATCH / 2], 1).astype(float)
        return f / (np.linalg.norm(f, axis=1, keepdims=True) + 1e-8), c, np.array([ow / size, oh / size])

    f0, c0, s0 = feats(img0); f1, c1, s1 = feats(img1)
    j = (f0 @ f1.T).argmax(1)                      # best panel patch per anchor patch
    return vote_assign(c0, c1[j], s0, s1, mask0, mask1), int(len(f0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="dinov2,loftr,roma,superpoint-lightglue")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--max-size", type=int, default=2048,
                    help="downscale each side so long edge <= this (keypoints scaled back to mask space)")
    ap.add_argument("--img-dir", default=str(MAPPING_OUTPUT_DIR_QUA / "images_segmentation"),
                    help="dir of whole-slide PNGs named {img}.png")
    ap.add_argument("--limit-pairs", type=int, default=0, help="0 = all WT pairs")
    args = ap.parse_args()
    models = [m.strip() for m in args.models.split(",") if m.strip()]

    pairs = wt_pairs()
    if args.limit_pairs:
        pairs = pairs[:args.limit_pairs]
    print(f"{len(pairs)} WT pairs | models={models} | device={args.device} | "
          f"max_size={args.max_size} | img_dir={args.img_dir}", flush=True)

    rows = []
    for i, (img1, img2, _) in enumerate(pairs, 1):
        mask0, mask1 = load_mask(img1), load_mask(img2)
        for m in models:
            (OUT / m).mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            try:
                if m == "dinov2":
                    assign, nk = run_dinov2(args.img_dir, img1, img2, mask0, mask1, args.device, args.max_size)
                else:
                    assign, nk = run_vismatch(m, args.img_dir, img1, img2, mask0, mask1, args.device, args.max_size)
                err = ""
            except Exception as e:  # noqa: BLE001
                assign, nk, err = {}, 0, f"{type(e).__name__}: {e}"
            dt = time.time() - t0
            with open(OUT / m / f"{img1}___vs___{img2}.pkl", "wb") as f:
                pickle.dump(assign, f)
            rows.append({"img1": img1, "img2": img2, "method": m, "seconds": round(dt, 2),
                         "n_kpts": nk, "n_assigned": len(assign), "max_size": args.max_size, "error": err})
            print(f"[{i:2d}/{len(pairs)}] {m:22s} {dt:7.1f}s  assigned={len(assign):6d}  "
                  f"kpts={nk}{('  ERR '+err) if err else ''}", flush=True)
        del mask0, mask1

    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    rt = OUT / "runtime.csv"
    if rt.exists():  # append across runs (different --models)
        df = pd.concat([pd.read_csv(rt), df], ignore_index=True)
        df = df.drop_duplicates(["img1", "img2", "method"], keep="last")
    df.to_csv(rt, index=False)
    print(f"\nsaved assignments under {OUT}/<method>/ and runtime -> {rt}")
    print(df.groupby("method")[["seconds", "n_assigned"]].mean().round(1).to_string())


if __name__ == "__main__":
    main()
