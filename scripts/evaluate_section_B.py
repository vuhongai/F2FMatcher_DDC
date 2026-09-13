"""Section B — F2FMatcher method evaluation (WT cohort only).

Everything is derived from PRECOMPUTED artifacts (no re-matching, no re-clustering):
  - results/QUA/eval/fiber_morphology.pkl   per-fiber morphology (scale-invariant + area + centroid)
  - results/QUA/prediction_output/*/paired_labels.pkl   F2FMatcher matches
  - results/QUA/mapping_summary.csv         per-pair coverage

Produces (all under results/QUA/eval/ + visualizations/eval/):
  B1  coverage_WT.csv / coverage_WT_by_panel.csv        per-pair + per-panel coverage (WT)
  B2  morph_consistency.csv / crosspanel_consistency.csv
      self-supervised accuracy:
        (a) matched-vs-random morphological consistency (scale-invariant shape)
        (b) cross-panel spatial consistency (star-topology cycle check)
  B3  baselines.csv / baseline_compare.csv
      F2FMatcher vs global-affine+kNN vs no-align kNN, scored on the same
      morphological-consistency metric
  section_B_results.json                            all headline numbers

Usage:  python scripts/evaluate_section_B.py
"""
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from sklearn.metrics import roc_auc_score

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, str(REPO))
from config.ddc_config import PAIR_DIRS_BASE_QUA  # noqa: E402

EVAL = REPO / "results" / "QUA" / "eval"
FIG = REPO / "visualizations" / "eval"
MORPH = EVAL / "fiber_morphology.pkl"
SUMMARY = REPO / "results" / "QUA" / "mapping_summary.csv"

WT_SAMPLES = ["QUAG01", "QUAG02", "QUAG03", "QUAG04", "QUAG05"]
SI_FEATURES = ["eccentricity", "solidity", "extent", "roundness", "aspect_ratio"]
RNG = np.random.default_rng(42)


# --------------------------------------------------------------------------- #
# loading + numpy fast-structures
# --------------------------------------------------------------------------- #
def wt_pairs():
    """[(sample, panel, img1, img2, pair_dir)] for WT pairs with labels."""
    out = []
    for d in sorted(PAIR_DIRS_BASE_QUA.iterdir()):
        if not d.is_dir() or "___vs___" not in d.name:
            continue
        if not any(s in d.name for s in WT_SAMPLES):
            continue
        if not (d / "paired_labels.pkl").exists():
            continue
        img1, img2 = d.name.split("___vs___")
        sample = next(s for s in WT_SAMPLES if s in d.name)
        panel = "_".join(img2.split("_")[1:-1])
        out.append((sample, panel, img1, img2, d))
    return out


def load_labels(pair_dir):
    with open(pair_dir / "paired_labels.pkl", "rb") as f:
        return pickle.load(f)


def build_arrays(morph):
    """img -> dict(labels, si (n x 5 standardized), cent (n x 2), pos {label: i})."""
    allf = np.vstack([df[SI_FEATURES].to_numpy(dtype=float) for df in morph.values()])
    mu, sd = allf.mean(0), allf.std(0)
    arr = {}
    for img, df in morph.items():
        labels = df.index.to_numpy()
        si = (df[SI_FEATURES].to_numpy(dtype=float) - mu) / sd
        cent = df[["cx", "cy"]].to_numpy(dtype=float)
        arr[img] = {"labels": labels, "si": si, "cent": cent,
                    "pos": {int(l): i for i, l in enumerate(labels)}}
    return arr, mu, sd


# --------------------------------------------------------------------------- #
# B1 — coverage
# --------------------------------------------------------------------------- #
def b1_coverage():
    df = pd.read_csv(SUMMARY)
    df["sample"] = df["img1"].str.extract(r"(QUAG\d+)")
    wt = df[df["sample"].isin(WT_SAMPLES)].copy()
    ok = wt.dropna(subset=["n_pairs"]).copy()
    ok["coverage_mean"] = (ok["coverage1"] + ok["coverage2"]) / 2

    per_pair = ok[["sample", "panel", "n_rois1", "n_rois2", "n_pairs",
                   "coverage1", "coverage2", "coverage_mean"]].reset_index(drop=True)
    per_pair.to_csv(EVAL / "coverage_WT.csv", index=False)

    by_panel = (ok.groupby("panel")
                   .agg(n_pairs=("n_pairs", "size"),
                        mean_pairs=("n_pairs", "mean"),
                        cov_ref=("coverage1", "mean"),
                        cov_other=("coverage2", "mean"),
                        cov_mean=("coverage_mean", "mean"))
                   .round(4).reset_index())
    by_panel.to_csv(EVAL / "coverage_WT_by_panel.csv", index=False)

    overall = {
        "n_samples": int(ok["sample"].nunique()),
        "n_pairs": int(len(ok)),
        "total_matched_pairs": float(ok["n_pairs"].sum()),
        "mean_coverage_ref": float(ok["coverage1"].mean()),
        "mean_coverage_other": float(ok["coverage2"].mean()),
        "mean_coverage_both": float(ok["coverage_mean"].mean()),
    }
    print("\n[B1] coverage (WT)")
    print(by_panel.to_string(index=False))
    print(f"    overall: {overall['n_pairs']} pairs, "
          f"cov ref {overall['mean_coverage_ref']:.1%}, "
          f"other {overall['mean_coverage_other']:.1%}, "
          f"both {overall['mean_coverage_both']:.1%}")
    return per_pair, by_panel, overall


# --------------------------------------------------------------------------- #
# B2a — morphological consistency (matched vs random)
# --------------------------------------------------------------------------- #
def _si_dist(si1, si2):
    d = si1 - si2
    return np.sqrt((d * d).sum(axis=1))


def b2_morph_consistency(arr, n_random=20000):
    rows = []
    for sample, panel, img1, img2, d in wt_pairs():
        a1, a2 = arr[img1], arr[img2]
        labels = load_labels(d)
        p1 = a1["pos"]; p2 = a2["pos"]
        i1 = np.array([p1[l1] for l1, l2 in labels if l1 in p1 and l2 in p2])
        i2 = np.array([p2[l2] for l1, l2 in labels if l1 in p1 and l2 in p2])
        dm = _si_dist(a1["si"][i1], a2["si"][i2])
        # random null: uniform pairs from the two pools
        n = min(n_random, len(a1["labels"]) * len(a2["labels"]))
        r1 = RNG.integers(0, len(a1["labels"]), n)
        r2 = RNG.integers(0, len(a2["labels"]), n)
        dr = _si_dist(a1["si"][r1], a2["si"][r2])
        y = np.concatenate([np.ones(len(dm)), np.zeros(len(dr))])
        auc = roc_auc_score(y, -np.concatenate([dm, dr]))
        rows.append({
            "sample": sample, "panel": panel, "n_matched": int(len(dm)),
            "matched_mean": float(dm.mean()), "matched_med": float(np.median(dm)),
            "random_mean": float(dr.mean()), "random_med": float(np.median(dr)),
            "auc": float(auc), "lift": float((dr.mean() - dm.mean()) / dr.mean()),
        })
    res = pd.DataFrame(rows)
    res.to_csv(EVAL / "morph_consistency.csv", index=False)
    pooled = {
        "n_pairs": int(len(res)),
        "n_matched_total": int(res["n_matched"].sum()),
        "matched_mean": float(res["matched_mean"].mean()),
        "matched_median": float(res["matched_med"].median()),
        "random_mean": float(res["random_mean"].mean()),
        "auc_mean": float(res["auc"].mean()),
        "auc_min": float(res["auc"].min()),
        "lift_mean": float(res["lift"].mean()),
    }
    print("\n[B2a] morphological consistency (scale-invariant shape), matched vs random")
    print(res[["sample", "panel", "n_matched", "matched_mean", "random_mean", "auc", "lift"]]
          .round(4).to_string(index=False))
    print(f"    pooled: AUC {pooled['auc_mean']:.3f} (min {pooled['auc_min']:.3f}), "
          f"matched {pooled['matched_mean']:.3f} vs random {pooled['random_mean']:.3f} "
          f"({pooled['lift_mean']:.1%} reduction)")
    return res, pooled


# --------------------------------------------------------------------------- #
# B2b — cross-panel spatial consistency (star-topology cycle check)
# --------------------------------------------------------------------------- #
def fit_affine(src, dst, n_rounds=3, z_thr=4.0):
    """Robust affine (linear+translation) src->dst from PAIRED correspondences.

    Iteratively reweighted least squares: fit, drop residuals > z_thr*MAD, refit.
    """
    src = np.asarray(src, float); dst = np.asarray(dst, float)
    for _ in range(n_rounds):
        # linear least squares: [x, y, 1] @ [A; t] = q   (per output coordinate)
        n = len(src)
        M = np.column_stack([src, np.ones(n)])
        sol, *_ = np.linalg.lstsq(M, dst, rcond=None)
        A, t = sol[:2, :], sol[2, :]
        pred = src @ A.T + t
        resid = np.linalg.norm(pred - dst, axis=1)
        mad = np.median(np.abs(resid - np.median(resid)))
        keep = resid < (z_thr * 1.4826 * mad + 1e-6)
        if keep.sum() < 10 or keep.all():
            break
        src, dst = src[keep], dst[keep]
    M = np.column_stack([src, np.ones(len(src))])
    sol, *_ = np.linalg.lstsq(M, dst, rcond=None)
    return sol[:2, :], sol[2, :]


def apply_affine(A, t, pts):
    return pts @ A.T + np.asarray(t)


def b2_crosspanel(arr, n_null=3000):
    pairs_by_sample = {}
    for sample, panel, img1, img2, d in wt_pairs():
        pairs_by_sample.setdefault(sample, {})[panel] = (img1, img2, load_labels(d))

    rows = []
    for sample, panels in pairs_by_sample.items():
        anchor = next(iter(panels.values()))[0]
        aA = arr[anchor]
        aff, matched = {}, {}
        for panel, (img1, img2, labels) in panels.items():
            a2 = arr[img2]
            p2 = a2["pos"]
            ok = [(l1, l2) for l1, l2 in labels if l1 in aA["pos"] and l2 in p2]
            if len(ok) < 50:
                continue
            p1 = aA["cent"][[aA["pos"][l1] for l1, _ in ok]]
            p2c = a2["cent"][[p2[l2] for _, l2 in ok]]
            A, t = fit_affine(p1, p2c)
            invA = np.linalg.inv(A)
            aff[panel] = (invA, -invA @ t, a2)
            # map each anchor label -> panel label (for this panel)
            mp = {l1: l2 for l1, l2 in ok}
            matched[panel] = mp
        panel_names = list(aff.keys())
        if len(panel_names) < 2:
            continue
        # real: for each anchor fiber matched in >=2 panels, cross-panel spread
        real = []
        for l1 in aA["labels"]:
            i1 = aA["pos"][l1]
            pts = []
            for panel in panel_names:
                invA, ti, a2 = aff[panel]
                l2 = matched[panel].get(l1)
                if l2 is None:
                    continue
                c2 = a2["cent"][a2["pos"][l2]]
                pts.append(apply_affine(invA, ti, c2))
            if len(pts) >= 2:
                real.append(pdist(np.array(pts)).mean())
        # null: random panel counterparts
        null = []
        n_fib = len(aA["labels"])
        for _ in range(n_null):
            pts = []
            for panel in panel_names:
                invA, ti, a2 = aff[panel]
                c2 = a2["cent"][RNG.integers(0, len(a2["cent"]))]
                pts.append(apply_affine(invA, ti, c2))
            null.append(pdist(np.array(pts)).mean())
        real = np.array(real); null = np.array(null)
        rows.append({
            "sample": sample, "n_anchor_fibers": int(len(real)), "n_panels": len(panel_names),
            "real_spread_mean": float(real.mean()), "real_spread_med": float(np.median(real)),
            "null_spread_mean": float(null.mean()), "null_spread_med": float(np.median(null)),
            "fold_reduction": float(null.mean() / max(real.mean(), 1e-9)),
        })
    res = pd.DataFrame(rows)
    res.to_csv(EVAL / "crosspanel_consistency.csv", index=False)
    print("\n[B2b] cross-panel spatial consistency (star-topology cycle check)")
    print(res.round(2).to_string(index=False))
    pooled = {
        "real_spread_mean": float(res["real_spread_mean"].mean()),
        "null_spread_mean": float(res["null_spread_mean"].mean()),
        "fold_reduction": float((res["null_spread_mean"] / res["real_spread_mean"]).mean()),
    }
    print(f"    pooled: real spread {pooled['real_spread_mean']:.1f}px vs null "
          f"{pooled['null_spread_mean']:.1f}px  ({pooled['fold_reduction']:.1f}x tighter)")
    return res, pooled


# --------------------------------------------------------------------------- #
# B3 — baselines vs F2FMatcher
# --------------------------------------------------------------------------- #
def b3_baselines(arr, n_null=5000):
    """Score F2FMatcher and two geometric baselines on the SAME anchor fibers
    (those F2FMatcher matched), using the morphological-consistency metric."""
    rows = []
    for sample, panel, img1, img2, d in wt_pairs():
        a1, a2 = arr[img1], arr[img2]
        c1, c2 = a1["cent"], a2["cent"]
        p1, p2 = a1["pos"], a2["pos"]
        tree2 = cKDTree(c2)
        # F2FMatcher correspondences
        labels = load_labels(d)
        i1 = np.array([p1[l1] for l1, l2 in labels if l1 in p1 and l2 in p2])
        i2 = np.array([p2[l2] for l1, l2 in labels if l1 in p1 and l2 in p2])
        dm_f2f = _si_dist(a1["si"][i1], a2["si"][i2])
        # Baseline A: global affine (robust fit on correspondences) + kNN
        A, t = fit_affine(c1[i1], c2[i2])
        pred2 = apply_affine(A, t, c1[i1])
        nn_aff = tree2.query(pred2, k=1)[1]
        dm_aff = _si_dist(a1["si"][i1], a2["si"][nn_aff])
        # Baseline B: no-alignment kNN (raw coords)
        nn_raw = tree2.query(c1[i1], k=1)[1]
        dm_raw = _si_dist(a1["si"][i1], a2["si"][nn_raw])
        rows.append({
            "sample": sample, "panel": panel, "f2f_n": int(len(dm_f2f)),
            "f2f_mean": float(dm_f2f.mean()),
            "affine_mean": float(dm_aff.mean()), "raw_mean": float(dm_raw.mean()),
        })
    res = pd.DataFrame(rows)
    res.to_csv(EVAL / "baselines.csv", index=False)
    # shared random null
    nulls = []
    for sample, panel, img1, img2, d in wt_pairs():
        a1, a2 = arr[img1], arr[img2]
        r1 = RNG.integers(0, len(a1["labels"]), n_null)
        r2 = RNG.integers(0, len(a2["labels"]), n_null)
        nulls.append(_si_dist(a1["si"][r1], a2["si"][r2]).mean())
    null_mean = float(np.mean(nulls))
    cmp = pd.DataFrame({
        "method": ["F2FMatcher", "global-affine + kNN", "no-align kNN", "random null"],
        "mean_morph_dist": [res["f2f_mean"].mean(), res["affine_mean"].mean(),
                            res["raw_mean"].mean(), null_mean],
    })
    cmp["reduction_vs_random"] = (null_mean - cmp["mean_morph_dist"]) / null_mean
    cmp.to_csv(EVAL / "baseline_compare.csv", index=False)
    print("\n[B3] baselines vs F2FMatcher (morphological consistency, lower=better)")
    print(cmp.round(4).to_string(index=False))
    return res, cmp, null_mean


# --------------------------------------------------------------------------- #
def main():
    t0 = time.time()
    FIG.mkdir(parents=True, exist_ok=True)
    with open(MORPH, "rb") as f:
        morph = pickle.load(f)
    arr, mu, sd = build_arrays(morph)
    print(f"loaded morphology: {len(arr)} images, "
          f"{sum(len(a['labels']) for a in arr.values()):,} fibers")

    results = {}
    results["B1_coverage"] = b1_coverage()[2]
    results["B2a_morph_consistency"] = b2_morph_consistency(arr)[1]
    results["B2b_crosspanel"] = b2_crosspanel(arr)[1]
    _, cmp, _ = b3_baselines(arr)
    results["B3_baselines"] = {
        "F2FMatcher": float(cmp.iloc[0]["mean_morph_dist"]),
        "global_affine_kNN": float(cmp.iloc[1]["mean_morph_dist"]),
        "no_align_kNN": float(cmp.iloc[2]["mean_morph_dist"]),
        "random_null": float(cmp.iloc[3]["mean_morph_dist"]),
    }
    results["runtime_s"] = time.time() - t0
    with open(EVAL / "section_B_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {EVAL/'section_B_results.json'}  ({results['runtime_s']:.0f}s)")


if __name__ == "__main__":
    main()
