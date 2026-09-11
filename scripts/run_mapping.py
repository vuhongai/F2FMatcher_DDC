"""F2FMatcher mapping runner with a segmentation/matching pipeline.

Producer threads (n_seg_workers) precompute the PERSISTENT part of each
pair — PNG export + CellPose segmentation (masks/flows pkls, cached) —
while the main thread (consumer) runs the matching stage on the previous
pair: 256px ROI crops (npz) + VAE embeddings (transient, deleted after
each pair) + geometry-aware matching. Segmentation time is thus hidden
behind matching time.

Resume: a pair is skipped entirely if its paired_labels.pkl already exists;
segmentation is skipped per image if its mask pkl already exists.
"""
import argparse
import os
import pickle
import queue
import sys
import threading
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from joblib import Parallel, delayed

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(PROJECT_ROOT))
from config.ddc_config import *

from f2fmatcher.config import load_config
from f2fmatcher.io.czi_reader import find_img_path, import_resize_export_czi, resize_image
from f2fmatcher.segmentation.cellpose_seg import (
    segment_image, filter_ROIs, _save_single_vae_input,
)
from f2fmatcher.vae.embed import generate_embedding_dir
from f2fmatcher.matching.matcher import match_fibers

PARAM_REF = ["fluorescence", "10X", 1.0]


def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", flush=True)


def cell(row, col, default):
    if col in row.index and pd.notna(row[col]):
        return str(row[col]).strip()
    return default


def load_pairs(csv_path):
    df = pd.read_csv(csv_path)
    if df.shape[1] < 6:
        raise ValueError(
            f"Pairing CSV needs at least 6 columns "
            f"(img1, img2, source1, source2, seg_channel1, seg_channel2), "
            f"got {df.shape[1]}: {list(df.columns)}"
        )
    pairs = []
    for _, row in df.iterrows():
        if pd.isna(row.iloc[0]) or pd.isna(row.iloc[1]):
            continue  # skip empty placeholder rows
        pairs.append({
            "img1": str(row.iloc[0]).strip(),
            "img2": str(row.iloc[1]).strip(),
            "source1": str(row.iloc[2]).strip(),
            "source2": str(row.iloc[3]).strip(),
            "ch1": int(row.iloc[4]),
            "ch2": int(row.iloc[5]),
            "cp1": cell(row, "cp_model1", DEFAULT_CP_MODEL),
            "cp2": cell(row, "cp_model2", DEFAULT_CP_MODEL),
            "param1": cell(row, "param1", DEFAULT_PARAM),
            "param2": cell(row, "param2", DEFAULT_PARAM),
            "obj1": cell(row, "obj1", DEFAULT_OBJECTIVE),
            "obj2": cell(row, "obj2", DEFAULT_OBJECTIVE),
        })
    return pairs


def validate_pairs(pairs):
    """Resolve each image path on disk. Returns a list of problems (empty if OK)."""
    problems = []
    for p in pairs:
        for side in ("1", "2"):
            src = Path(p[f"source{side}"])
            if not src.is_dir():
                problems.append(f"source{side} dir not found: {src}")
                continue
            try:
                path = find_img_path(p[f"img{side}"], str(src))
            except FileNotFoundError as e:
                problems.append(str(e))
                continue
            p[f"path{side}"] = path
            p[f"czi{side}"] = path.lower().endswith(".czi")
    return problems


def clean_transient_files(dir_output):
    """Remove leftover ROI crops/embeddings from a previously interrupted run."""
    n = 0
    for sub, pattern in [("npz_256", "*_roi_*.npz"), ("VAE_embed", "*_roi_*.npy")]:
        d = Path(dir_output) / sub
        if d.is_dir():
            for f in d.glob(pattern):
                f.unlink()
                n += 1
    if n:
        log(f"Removed {n} leftover transient file(s) from a previous interrupted run.")


def pick_device():
    """Use a GPU only if it has >= 8 GB free (the local LLM may occupy the GPUs)."""
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            try:
                free, _ = torch.cuda.mem_get_info(i)
            except Exception:
                continue  # probe allocation failed -> GPU saturated
            if free >= 8 * 1024**3:
                return torch.device(f"cuda:{i}")
    return torch.device("cpu")


def pair_dir_of(p):
    return Path(MAPPING_OUTPUT_DIR_QUA) / "prediction_output" / f"{p['img1']}___vs___{p['img2']}"


def export_png(p, side, images_dir):
    img, src, czi, param, obj, ch = (
        p[f"img{side}"], p[f"source{side}"], p[f"czi{side}"],
        p[f"param{side}"], p[f"obj{side}"], p[f"ch{side}"],
    )
    png = Path(images_dir) / f"{img}.png"
    if png.exists():
        return
    param_img = [param, obj, 1.0]
    if czi:
        import_resize_export_czi(find_img_path(img, src), param == "fluorescence",
                                 ch, str(images_dir), param_img, PARAM_REF)
    else:
        resize_image(find_img_path(img, src), str(images_dir), param_img, PARAM_REF)
    log(f"  exported PNG {img}")


def ensure_vae_inputs(img_name, cp_model_name, images_dir, dir_npz, dir_cp_masks,
                      n_processes, gpu, cp_model_path):
    """Segment (or load cached masks/flows) + write 256px ROI npz crops (skip existing)."""
    seg_kwargs = {}
    if cp_model_path is not None:
        seg_kwargs["CP_model_path"] = cp_model_path
    cp_output = segment_image(img_path=f"{images_dir}/{img_name}.png",
                              CP_model_name=cp_model_name,
                              savedir=str(dir_cp_masks), gpu=gpu, **seg_kwargs)
    masks, props, img_rec = cp_output
    labels = filter_ROIs(cp_output, 256)
    missing = [l for l in labels if not os.path.exists(f"{dir_npz}/{img_name}_roi_{l}.npz")]
    if missing:
        log(f"  {img_name}: {len(labels)} ROIs, generating {len(missing)} npz crop(s)")
        centroids_yx = {pr.label: pr.centroid for pr in props if pr.label in labels}
        Parallel(n_jobs=n_processes)(
            delayed(_save_single_vae_input)(lbl, centroids_yx, masks, img_rec, img_name,
                                            str(dir_npz), 256)
            for lbl in missing
        )
    else:
        log(f"  {img_name}: {len(labels)} ROIs, npz crops already present")
    del img_rec  # free the big flow array; matching only needs props
    return cp_output, labels


def segment_pair(p, config, device):
    """Producer stage: PNG export + CellPose segmentation (persistent pkls) for one pair."""
    dir_output = Path(MAPPING_OUTPUT_DIR_QUA)
    images_dir = dir_output / "images_segmentation"
    dir_cp = dir_output / "out_CP_masks"
    for d in (images_dir, dir_cp):
        d.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    log(f"SEGMENT {p['img1']} vs {p['img2']}")
    for side in ("1", "2"):
        export_png(p, side, images_dir)
    cp_model_path = config.get("cellpose.model_path")
    seg_kwargs = {"CP_model_path": cp_model_path} if cp_model_path is not None else {}
    for side in ("1", "2"):
        img, cp_model = p[f"img{side}"], p[f"cp{side}"]
        segment_image(img_path=f"{images_dir}/{img}.png", CP_model_name=cp_model,
                      savedir=str(dir_cp), gpu=device.type == "cuda", **seg_kwargs)
        log(f"  segmented {img}")
    log(f"SEGMENTED {p['img2']}: {time.time()-t0:.0f}s")


def match_pair(p, config, device):
    """Consumer stage: npz crops + VAE embeddings (transient) + matching + save + cleanup."""
    dir_output = Path(MAPPING_OUTPUT_DIR_QUA)
    images_dir = dir_output / "images_segmentation"
    dir_npz = dir_output / "npz_256"
    dir_embed = dir_output / "VAE_embed"
    dir_cp = dir_output / "out_CP_masks"
    for d in (dir_npz, dir_embed):
        d.mkdir(parents=True, exist_ok=True)
    pair_dir = pair_dir_of(p)
    pair_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    log(f"MATCH {p['img1']} vs {p['img2']}")
    cp_model_path = config.get("cellpose.model_path")
    cp1, labels1 = ensure_vae_inputs(p["img1"], p["cp1"], images_dir, dir_npz, dir_cp,
                                     n_processes_seg, device.type == "cuda", cp_model_path)
    cp2, labels2 = ensure_vae_inputs(p["img2"], p["cp2"], images_dir, dir_npz, dir_cp,
                                     n_processes_seg, device.type == "cuda", cp_model_path)
    generate_embedding_dir(str(dir_npz), str(dir_embed), config.get("vae.checkpoint"),
                           device, img_names=[p["img1"], p["img2"]])
    label2index = {
        "img1": {l: i for i, l in enumerate(labels1)},
        "img2": {l: i for i, l in enumerate(labels2)},
    }
    cls_ckpt = config.get("classifier.checkpoint")
    matching_cfg = config.get("matching", {})
    if isinstance(matching_cfg, dict):
        matching_cfg = type("obj", (object,), matching_cfg)

    matched_labels, scores, spatial_dist, cp_outputs = match_fibers(
        p["img1"], labels1, cp1,
        p["img2"], labels2, cp2,
        label2index, str(dir_embed), cls_ckpt, device,
        list_k=matching_cfg.get("list_k", [3, 5, 7]) if hasattr(matching_cfg, "get") else matching_cfg.list_k,
        n_initial_guess=getattr(matching_cfg, "n_initial_guess", 40),
        n_pair_selected=getattr(matching_cfg, "n_pair_selected", 3),
        min_cls_logit_init=getattr(matching_cfg, "min_cls_logit_init", 0.75),
        distance_neighbors_ref=getattr(matching_cfg, "distance_neighbors_ref", 200),
        max_distance_affine=getattr(matching_cfg, "max_distance_affine", 150),
        max_cost_geo_neighbors_sides=getattr(matching_cfg, "max_cost_geo_neighbors_sides", 30),
        max_cost_geo_neighbors_angles=getattr(matching_cfg, "max_cost_geo_neighbors_angles", 0.15),
        min_cls_logit=getattr(matching_cfg, "min_cls_logit", 0.5),
        patience_label=getattr(matching_cfg, "patience_label", 5),
        n_neighbors_validation=getattr(matching_cfg, "n_neighbors_validation", 3),
        n_processes=n_processes,
        n_try_unannotated=getattr(matching_cfg, "n_try_unannotated", 1),
        patience_prediction_neighbors=getattr(matching_cfg, "patience_prediction_neighbors", 5),
        use_multiprocessing_for_local_prediction=getattr(matching_cfg, "use_multiprocessing", True),
        save_step_prediction=save_step_prediction,
        dir_save_prediction_output=str(pair_dir),
    )
    with open(f"{pair_dir}/paired_labels.pkl", "wb") as f:
        pickle.dump(matched_labels, f)

    # cleanup: drop all transient files of this pair (npz, embeddings, scores cache)
    for f in os.listdir(dir_npz):
        if f.split("_roi_")[0] in (p["img1"], p["img2"]):
            os.remove(os.path.join(dir_npz, f))
    for f in os.listdir(dir_embed):
        if f.split("_roi_")[0] in (p["img1"], p["img2"]) \
                or f == f"{p['img1']}_vs_{p['img2']}_scores.npy":
            os.remove(os.path.join(dir_embed, f))

    log(f"MATCHED {p['img2']}: {len(matched_labels)} pairs in {time.time()-t0:.0f}s")
    return len(matched_labels)


def main():
    parser = argparse.ArgumentParser(
        description="Run F2FMatcher fiber mapping for all pairs in the pairing CSV "
                    "(segmentation pipelined with matching)."
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate all pairs and print the run plan, then exit (no GPU needed).")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only the first N pairs.")
    parser.add_argument("--csv", type=str, default=None,
                        help="Pairing CSV path (overrides CSV_FILE in config).")
    args = parser.parse_args()

    csv_path = Path(args.csv) if args.csv else Path(CSV_FILE)
    pairs = load_pairs(csv_path)
    if args.limit is not None:
        pairs = pairs[: args.limit]

    log(f"Loaded {len(pairs)} pair(s) from {csv_path}")
    problems = validate_pairs(pairs)
    if problems:
        log(f"Validation failed with {len(problems)} problem(s):")
        for pr in problems:
            log(f"  - {pr}")
        sys.exit(1)

    if args.dry_run:
        for i, p in enumerate(pairs, 1):
            print(f"[{i}/{len(pairs)}] {p['img1']}  vs  {p['img2']}")
            print(f"    czi1={p['czi1']} czi2={p['czi2']}  ch1={p['ch1']} ch2={p['ch2']}  "
                  f"cp1={p['cp1']} cp2={p['cp2']}  "
                  f"param1={p['param1']} param2={p['param2']}  obj1={p['obj1']} obj2={p['obj2']}")
        print("Dry run OK.")
        return

    config = load_config()
    device = pick_device()
    log(f"Device: {device} | n_seg_workers: {n_seg_workers} | "
        f"n_processes (matching): {n_processes} | output: {MAPPING_OUTPUT_DIR_QUA}")

    clean_transient_files(MAPPING_OUTPUT_DIR_QUA)

    todo = [p for p in pairs if not (pair_dir_of(p) / "paired_labels.pkl").exists()]
    n_done_already = len(pairs) - len(todo)
    if n_done_already:
        log(f"Resuming: {n_done_already} pair(s) already have paired_labels.pkl, skipped.")
    if not todo:
        log("Nothing to do.")
        return

    prep_q = queue.Queue()
    counter = 0
    counter_lock = threading.Lock()

    def producer():
        nonlocal counter
        while True:
            with counter_lock:
                if counter >= len(todo):
                    break
                idx = counter
                counter += 1
            p = todo[idx]
            try:
                segment_pair(p, config, device)
                prep_q.put((idx, None, None))
            except Exception as e:
                log(f"SEGMENT FAILED {p['img1']} vs {p['img2']}: {e}")
                traceback.print_exc()
                prep_q.put((idx, None, e))
        prep_q.put(None)  # sentinel

    threads = [threading.Thread(target=producer, name=f"prep-{i}", daemon=True)
               for i in range(n_seg_workers)]
    for t in threads:
        t.start()

    ok, failed = [], []
    sentinels = 0
    n_finished = 0
    while sentinels < n_seg_workers:
        item = prep_q.get()
        if item is None:
            sentinels += 1
            continue
        idx, _, err = item
        p = todo[idx]
        n_finished += 1
        tag = f"[{n_finished}/{len(todo)}]"
        if err is not None:
            failed.append((p["img1"], p["img2"], f"segment: {err}"))
            continue
        try:
            n = match_pair(p, config, device)
            ok.append((p["img1"], p["img2"], n))
            log(f"{tag} OK {p['img2']}: {n} pairs")
        except Exception as e:
            log(f"{tag} MATCH FAILED {p['img1']} vs {p['img2']}: {e}")
            traceback.print_exc()
            failed.append((p["img1"], p["img2"], f"match: {e}"))

    log(f"\nDone: {len(ok)} OK, {len(failed)} failed, {n_done_already} already done "
        f"({len(pairs)} total).")
    for a, b, why in failed:
        log(f"  failed: {a} vs {b} ({why})")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
