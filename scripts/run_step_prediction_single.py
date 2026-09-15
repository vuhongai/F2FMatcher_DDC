"""Re-run F2FMatcher matching for ONE pair with save_step_prediction=True.

Regenerates the transient npz crops + VAE embeddings, runs match_fibers with
step saving enabled, and writes step_prediction.pkl (+ paired_labels.pkl) into
the pair's prediction_output directory. The existing paired_labels.pkl is
backed up first, and the re-run result is verified against it (as a set).

Usage:
  python scripts/run_step_prediction_single.py [--img1 NAME] [--img2 NAME] [--device cpu|auto]
"""
import argparse
import os
import pickle
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import torch

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "scripts"))
from config.ddc_config import *  # noqa: F403

from run_mapping import (  # noqa: E402
    log, load_pairs, validate_pairs, ensure_vae_inputs,
)
from f2fmatcher.config import load_config  # noqa: E402
from f2fmatcher.vae.embed import generate_embedding_dir  # noqa: E402
from f2fmatcher.matching.matcher import match_fibers  # noqa: E402

DEFAULT_IMG1 = "22-082_10X_DAPI_LAM_DYS_COL4_5-Scene-1-QUAG01"
DEFAULT_IMG2 = "22-082_10X_NADH_13-Scene-1-QUAG01"

TA_IMG1 = "22-082_10X_DAPI_LAM_DYS_COL4_1-Scene-1-TAG01"
TA_IMG2 = "22-082_10X_NADH_9-Scene-1-TAG01"
TA_IMG1_05 = "22-082_10X_DAPI_LAM_DYS_COL4_1-Scene-5-TAG05"
TA_IMG2_05 = "22-082_10X_NADH_9-Scene-5-TAG05"

# Pairs not in the production CSV (e.g. the original 6-panel run used a wider CSV).
# Specs mirror config/ddc_config.py SLIDES entries.
FALLBACK_PAIRS = {
    (DEFAULT_IMG1, DEFAULT_IMG2): {
        "img1": DEFAULT_IMG1, "img2": DEFAULT_IMG2,
        "source1": "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_QUA/IHF_Lam-Dys-Col4",
        "source2": "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_QUA/NADH",
        "ch1": 1, "ch2": 0,
        "cp1": "CP_AV_Laminin_Dia_Qua_TA_AxioScan10X",
        "cp2": "CP_AV_TA_COX-SDH-NADH_AxioScan10X",
        "param1": "fluorescence", "param2": "brightfield",
        "obj1": "10X", "obj2": "10X",
    },
    (TA_IMG1, TA_IMG2): {
        "img1": TA_IMG1, "img2": TA_IMG2,
        "source1": "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_TA/IHF_Lam-Dys-Col4",
        "source2": "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_TA/NADH",
        "ch1": 1, "ch2": 0,
        "cp1": "CP_AV_Laminin_Dia_Qua_TA_AxioScan10X",
        "cp2": "CP_AV_TA_COX-SDH-NADH_AxioScan10X",
        "param1": "fluorescence", "param2": "brightfield",
        "obj1": "10X", "obj2": "10X",
    },
    (TA_IMG1_05, TA_IMG2_05): {
        "img1": TA_IMG1_05, "img2": TA_IMG2_05,
        "source1": "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_TA/IHF_Lam-Dys-Col4",
        "source2": "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_TA/NADH",
        "ch1": 1, "ch2": 0,
        "cp1": "CP_AV_Laminin_Dia_Qua_TA_AxioScan10X",
        "cp2": "CP_AV_TA_COX-SDH-NADH_AxioScan10X",
        "param1": "fluorescence", "param2": "brightfield",
        "obj1": "10X", "obj2": "10X",
    },
}


def pick_device(arg):
    if arg == "cpu":
        return torch.device("cpu")
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            try:
                free, _ = torch.cuda.mem_get_info(i)
            except Exception:
                continue
            if free >= 8 * 1024**3:
                return torch.device(f"cuda:{i}")
    return torch.device("cpu")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--img1", default=DEFAULT_IMG1)
    ap.add_argument("--img2", default=DEFAULT_IMG2)
    ap.add_argument("--device", default="cpu", choices=["cpu", "auto"])
    ap.add_argument("--data-root", default=str(MAPPING_OUTPUT_DIR_QUA),
                    help="results dir containing images_segmentation/, out_CP_masks/, prediction_output/ "
                         "(default: QUA; use results/TA for the TA cohort)")
    args = ap.parse_args()

    p = None
    pairs = load_pairs(CSV_FILE)
    problems = validate_pairs(pairs)
    if problems:
        log(f"Validation failed: {problems}")
        sys.exit(1)
    p = next((x for x in pairs if x["img1"] == args.img1 and x["img2"] == args.img2), None)
    if p is None and (args.img1, args.img2) in FALLBACK_PAIRS:
        p = FALLBACK_PAIRS[(args.img1, args.img2)]
        p["czi1"] = p["czi2"] = False
        log("pair not in CSV -> using FALLBACK_PAIRS spec (PNG + masks must be cached)")
    if p is None:
        sys.exit(f"pair {args.img1} vs {args.img2} not in {CSV_FILE} and no fallback spec")

    config = load_config()
    device = pick_device(args.device)
    log(f"device={device}")

    dir_output = Path(args.data_root)
    images_dir = dir_output / "images_segmentation"
    dir_npz = dir_output / "npz_256"
    dir_embed = dir_output / "VAE_embed"
    dir_cp = dir_output / "out_CP_masks"
    for d in (images_dir, dir_npz, dir_embed, dir_cp):
        d.mkdir(parents=True, exist_ok=True)
    pair_dir = dir_output / "prediction_output" / f"{p['img1']}___vs___{p['img2']}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    cp_model_path = config.get("cellpose.model_path")
    seg_kwargs = {"CP_model_path": cp_model_path} if cp_model_path is not None else {}

    log("SEGMENT/LOAD img1 (masks+flows cached, npz crops regenerated)")
    cp1, labels1 = ensure_vae_inputs(p["img1"], p["cp1"], images_dir, dir_npz, dir_cp,
                                     n_processes_seg, device.type == "cuda", cp_model_path)
    log("SEGMENT/LOAD img2")
    cp2, labels2 = ensure_vae_inputs(p["img2"], p["cp2"], images_dir, dir_npz, dir_cp,
                                     n_processes_seg, device.type == "cuda", cp_model_path)
    log(f"labels: img1={len(labels1)} img2={len(labels2)}  ({time.time()-t0:.0f}s)")

    log("VAE embeddings")
    generate_embedding_dir(str(dir_npz), str(dir_embed), config.get("vae.checkpoint"),
                           device, img_names=[p["img1"], p["img2"]])

    backup = pair_dir / "paired_labels_backup.pkl"
    existing = pair_dir / "paired_labels.pkl"
    if existing.exists():
        shutil.copy2(existing, backup)
        log(f"backed up existing paired_labels.pkl -> {backup.name}")

    label2index = {
        "img1": {l: i for i, l in enumerate(labels1)},
        "img2": {l: i for i, l in enumerate(labels2)},
    }
    cls_ckpt = config.get("classifier.checkpoint")
    matching_cfg = config.get("matching", {})
    if isinstance(matching_cfg, dict):
        matching_cfg = type("obj", (object,), matching_cfg)

    log("MATCH (save_step_prediction=True)")
    t1 = time.time()
    matched_labels, scores, spatial_dist, _ = match_fibers(
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
        save_step_prediction=True,
        dir_save_prediction_output=str(pair_dir),
    )
    with open(f"{pair_dir}/paired_labels.pkl", "wb") as f:
        pickle.dump(matched_labels, f)
    log(f"MATCHED: {len(matched_labels)} pairs in {time.time()-t1:.0f}s (total {time.time()-t0:.0f}s)")

    # cleanup transient files of this pair
    for f in os.listdir(dir_npz):
        if f.split("_roi_")[0] in (p["img1"], p["img2"]):
            os.remove(os.path.join(dir_npz, f))
    for f in os.listdir(dir_embed):
        if f.split("_roi_")[0] in (p["img1"], p["img2"]) \
                or f == f"{p['img1']}_vs_{p['img2']}_scores.npy":
            os.remove(os.path.join(dir_embed, f))

    # verify against backup
    if backup.exists():
        with open(backup, "rb") as f:
            old = pickle.load(f)
        same = set(map(tuple, old)) == set(map(tuple, matched_labels))
        log(f"VERIFICATION vs backup: {'IDENTICAL' if same else 'DIFFERENT'} "
            f"(old={len(old)}, new={len(matched_labels)}, "
            f"only_old={len(set(map(tuple, old)) - set(map(tuple, matched_labels)))}, "
            f"only_new={len(set(map(tuple, matched_labels)) - set(map(tuple, old)))})")

    step_pkl = pair_dir / "step_prediction.pkl"
    log(f"step_prediction.pkl: {step_pkl.exists()} ({step_pkl.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
