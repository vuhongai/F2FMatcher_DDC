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

Memory note (this box: 2x A30 24 GB, no xformers in the `vismatch` env): DINOv2
and RoMa use EAGER attention, whose (H x N x N) tensor at a 2048-px long edge is
~22 GB — an OOM even on a FREE 24 GB GPU. Both are therefore capped at
METHOD_MAX_SIZE (1024 px, ~1.3 GB attention) regardless of --max-size; the
keypoint/semi-dense methods (LoFTR, SuperPoint-LightGlue) run at the full
--max-size. This cap is part of the scalability story and is logged per method.

Device: --device auto (default) uses CUDA only when >= 12 GB is free on GPU 0
(the local LLM normally occupies both A30s); otherwise it falls back to CPU.
An explicit --device cuda|cpu overrides the probe.

Outputs (under results/QUA/eval/wholeslide/):
    {method}/{img1}___vs___{img2}.pkl   {anchor_label:int -> panel_label:int}
    runtime.csv                          per (pair, method): seconds, n_kpts, max_size

Usage (on the GPU box, WT data present):
    python scripts/wholeslide_run_baselines.py --models dinov2,loftr,roma,superpoint-lightglue
    python scripts/wholeslide_run_baselines.py --models dinov2 --device cpu --limit-pairs 2   # smoke test

Resume: pairs whose assignment pkl already exists and is non-empty are skipped
(override with --force). runtime.csv is rewritten after every (pair, method), so
a killed run keeps its data.
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

# Eager-attention (no xformers) memory ceiling on a 24 GB GPU:
# 2048-px long edge -> ~22 GB attention tensor -> OOM. 1024 px -> ~1.3 GB, safe.
METHOD_MAX_SIZE = {"dinov2": 1024, "roma": 1024}
CUDA_FREE_MIN = 12 * 1024**3  # bytes; enough for LoFTR@2048 + RoMa/DINOv2@1024


def pick_device(requested):
    """'auto' -> cuda if >= CUDA_FREE_MIN bytes free on GPU 0, else cpu."""
    if requested != "auto":
        return requested
    try:
        import torch
        if torch.cuda.is_available():
            free, _ = torch.cuda.mem_get_info(0)
            if free >= CUDA_FREE_MIN:
                return "cuda"
    except Exception:  # noqa: BLE001  (mem_get_info itself can raise when saturated)
        pass
    return "cpu"


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


def build_matcher(model_name, device):
    """Instantiate the (expensive) matcher once; reused across all pairs."""
    if model_name == "dinov2":
        import torch
        model = torch.hub.load("facebookresearch/dinov2", "dinov2_vits14",
                               pretrained=True).eval().to(device)
        return model
    from vismatch import get_matcher
    return get_matcher(model_name, device=device, max_num_keypoints=8192)


def target_size(ow, oh, max_size, mult):
    """Aspect-preserving target (w, h): long edge -> max_size, dims floored to mult.
    (Square resizing would anisotropically stretch sections by up to ~16%,
    distorting the very appearance the matchers rely on.)"""
    s = max(ow, oh) / float(max_size)
    if s <= 1:
        return ow, oh
    return int((ow / s) // mult) * mult, int((oh / s) // mult) * mult


def _vismatch_load(matcher, img_dir, img, max_size):
    """Load via the matcher's own loader (aspect-preserving resize, long edge <= max_size).
    Returns (image_tensor, scale resized->original)."""
    from PIL import Image as _I
    p = str(Path(img_dir) / f"{img}.png")
    ow, oh = _I.open(p).size
    rw, rh = target_size(ow, oh, max_size, 8)
    t = matcher.load_image(p, resize=(rh, rw))  # vismatch resize tuple is (H, W)
    rh, rw = int(t.shape[-2]), int(t.shape[-1])
    return t, np.array([ow / rw, oh / rh])


def run_vismatch(matcher, img_dir, img0, img1, mask0, mask1, max_size):
    a, s0 = _vismatch_load(matcher, img_dir, img0, max_size)
    b, s1 = _vismatch_load(matcher, img_dir, img1, max_size)
    res = matcher(a, b)
    k0, k1 = res["matched_kpts0"], res["matched_kpts1"]
    return vote_assign(k0, k1, s0, s1, mask0, mask1), int(len(k0))


def run_dinov2(model, img_dir, img0, img1, mask0, mask1, device, max_size):
    import torch

    def feats(img):
        from PIL import Image as _I
        ow, oh = _I.open(Path(img_dir) / f"{img}.png").size
        rw, rh = target_size(ow, oh, max_size, DINO_PATCH)  # multiple of patch
        im = _I.open(Path(img_dir) / f"{img}.png").convert("RGB").resize((rw, rh))
        x = (np.asarray(im, np.float32) / 255.0 - DINO_MEAN) / DINO_STD
        x = torch.from_numpy(x).permute(2, 0, 1).unsqueeze(0).to(device)
        with torch.no_grad():
            f = model.forward_features(x)["x_norm_patchtokens"][0].cpu().numpy()
        nw, nh = rw // DINO_PATCH, rh // DINO_PATCH
        gy, gx = np.meshgrid(np.arange(nh), np.arange(nw), indexing="ij")
        c = np.stack([gx.ravel() * DINO_PATCH + DINO_PATCH / 2,
                      gy.ravel() * DINO_PATCH + DINO_PATCH / 2], 1).astype(float)
        return f / (np.linalg.norm(f, axis=1, keepdims=True) + 1e-8), c, np.array([ow / rw, oh / rh])

    f0, c0, s0 = feats(img0); f1, c1, s1 = feats(img1)
    j = (f0 @ f1.T).argmax(1)                      # best panel patch per anchor patch
    return vote_assign(c0, c1[j], s0, s1, mask0, mask1), int(len(f0))


def save_runtime(rows):
    """Merge with any existing runtime.csv (append across --models runs)."""
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    rt = OUT / "runtime.csv"
    if rt.exists():
        df = pd.concat([pd.read_csv(rt), df], ignore_index=True)
        df = df.drop_duplicates(["img1", "img2", "method"], keep="last")
    df.to_csv(rt, index=False)
    return df


def load_existing(path, force):
    """Resume: return the saved assignment if present and non-empty, else None."""
    if force or not path.exists():
        return None
    try:
        with open(path, "rb") as f:
            d = pickle.load(f)
        return d if d else None
    except Exception:  # noqa: BLE001  (corrupt/partial file -> redo)
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="dinov2,loftr,roma,superpoint-lightglue")
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    ap.add_argument("--max-size", type=int, default=2048,
                    help="downscale each side so long edge <= this (keypoints scaled back to mask space); "
                         "dinov2/roma are additionally capped at 1024 (eager-attention memory)")
    ap.add_argument("--img-dir", default=str(MAPPING_OUTPUT_DIR_QUA / "images_segmentation"),
                    help="dir of whole-slide PNGs named {img}.png")
    ap.add_argument("--limit-pairs", type=int, default=0, help="0 = all WT pairs")
    ap.add_argument("--force", action="store_true", help="re-run pairs even if a non-empty pkl exists")
    args = ap.parse_args()
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    device = pick_device(args.device)

    pairs = wt_pairs()
    if args.limit_pairs:
        pairs = pairs[:args.limit_pairs]
    print(f"{len(pairs)} WT pairs | models={models} | device={device} "
          f"(requested {args.device}) | max_size={args.max_size} | img_dir={args.img_dir}",
          flush=True)
    for m in models:
        cap = METHOD_MAX_SIZE.get(m)
        if cap and cap < args.max_size:
            print(f"  note: {m} capped at {cap}px (eager-attention memory ceiling, see docstring)",
                  flush=True)

    # Build each (expensive) matcher once, before the pair loop.
    handles = {}
    for m in models:
        print(f"loading {m} on {device} ...", flush=True)
        handles[m] = build_matcher(m, device)
        print(f"  {m} ready", flush=True)

    rows = []
    n_done = n_skip = 0
    for i, (img1, img2, d) in enumerate(pairs, 1):
        mask0, mask1 = load_mask(img1), load_mask(img2)
        for m in models:
            (OUT / m).mkdir(parents=True, exist_ok=True)
            out_pkl = OUT / m / f"{img1}___vs___{img2}.pkl"
            existing = load_existing(out_pkl, args.force)
            if existing is not None:
                n_skip += 1
                print(f"[{i:2d}/{len(pairs)}] {m:22s} SKIP (exists, {len(existing)} assigned)",
                      flush=True)
                continue
            t0 = time.time()
            try:
                size = min(args.max_size, METHOD_MAX_SIZE.get(m, args.max_size))
                if m == "dinov2":
                    assign, nk = run_dinov2(handles[m], args.img_dir, img1, img2, mask0, mask1, device, size)
                else:
                    assign, nk = run_vismatch(handles[m], args.img_dir, img1, img2, mask0, mask1, size)
                err = ""
            except Exception as e:  # noqa: BLE001
                assign, nk, err = {}, 0, f"{type(e).__name__}: {e}"
            dt = time.time() - t0
            with open(out_pkl, "wb") as f:
                pickle.dump(assign, f)
            rows.append({"img1": img1, "img2": img2, "method": m, "seconds": round(dt, 2),
                         "n_kpts": nk, "n_assigned": len(assign),
                         "max_size": min(args.max_size, METHOD_MAX_SIZE.get(m, args.max_size)),
                         "device": device, "error": err})
            save_runtime(rows)
            n_done += 1
            print(f"[{i:2d}/{len(pairs)}] {m:22s} {dt:7.1f}s  assigned={len(assign):6d}  "
                  f"kpts={nk}{('  ERR '+err) if err else ''}", flush=True)
        del mask0, mask1

    rt = OUT / "runtime.csv"
    df = pd.read_csv(rt) if rt.exists() else pd.DataFrame()
    print(f"\nsaved assignments under {OUT}/<method>/ and runtime -> {rt} "
          f"({n_done} run, {n_skip} skipped)")
    if len(df):
        print(df.groupby("method")[["seconds", "n_assigned"]].mean().round(1).to_string())


if __name__ == "__main__":
    main()
