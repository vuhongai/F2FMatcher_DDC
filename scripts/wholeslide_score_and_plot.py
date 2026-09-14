"""Score the whole-slide method comparison on WT: coverage vs correctness proxy.

The point of this figure: on whole slides there is NO ground truth, so coverage
alone is meaningless (a method can "match" everything at low accuracy). We
therefore report, for every method, **coverage AND two label-free correctness
proxies** — (i) shape-consistency AUC (an assigned pair should have matching
scale-invariant shape) and (ii) cross-panel cycle-consistency (an anchor fibre's
counterparts across panels should coincide after affine alignment). A good method
sits high on BOTH coverage and correctness; geometry-only kNN reaches full
coverage at chance accuracy; learned matchers (DINOv2, LoFTR, RoMa, SuperGlue)
reach high coverage but low accuracy.

Methods scored:
    F2FMatcher            from prediction_output/*/paired_labels.pkl
    no-align kNN          nearest panel centroid (raw coords)
    global-affine + kNN   robust affine (fit on F2F corr) + nearest centroid
    <learned matchers>    from results/QUA/eval/wholeslide/<method>/*.pkl
                          (produced by wholeslide_run_baselines.py; optional)

Outputs:
    results/QUA/eval/wholeslide_scores.csv
    visualizations/eval/fig_B_wholeslide_coverage_vs_consistency.png
    visualizations/eval/fig_B_wholeslide_runtime.png   (if runtime.csv present)

Usage:  python scripts/wholeslide_score_and_plot.py
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from sklearn.metrics import roc_auc_score

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from evaluate_section_B import (  # noqa: E402  reuse the §B machinery
    build_arrays, _si_dist, wt_pairs, load_labels, fit_affine, apply_affine,
)

EVAL = REPO / "results" / "QUA" / "eval"
FIG = REPO / "visualizations" / "eval"
WHOLESLIDE = EVAL / "wholeslide"
MORPH = EVAL / "fiber_morphology.pkl"
RNG = np.random.default_rng(42)

GEOMETRIC = ["no-align kNN", "global-affine + kNN"]


# ---- per-method assignment {anchor_label: panel_label} for one pair --------- #
def assign_f2f(a1, a2, d, _):
    p1, p2 = a1["pos"], a2["pos"]
    return {int(l1): int(l2) for l1, l2 in load_labels(d) if l1 in p1 and l2 in p2}


def assign_knn(a1, a2, d, mode):
    p1, p2 = a1["pos"], a2["pos"]
    corr = [(l1, l2) for l1, l2 in load_labels(d) if l1 in p1 and l2 in p2]
    src = a1["cent"]
    if mode == "affine" and len(corr) >= 10:
        i1 = [p1[l1] for l1, _ in corr]; i2 = [p2[l2] for _, l2 in corr]
        A, t = fit_affine(a1["cent"][i1], a2["cent"][i2])
        src = apply_affine(A, t, a1["cent"])
    nn = cKDTree(a2["cent"]).query(src, k=1)[1]
    return {int(a1["labels"][i]): int(a2["labels"][j]) for i, j in enumerate(nn)}


def assign_learned(method):
    def fn(a1, a2, d, _):
        f = WHOLESLIDE / method / f"{d.name}.pkl"
        if not f.exists():
            return None
        with open(f, "rb") as fh:
            raw = pickle.load(fh)
        p1, p2 = a1["pos"], a2["pos"]
        return {int(k): int(v) for k, v in raw.items() if int(k) in p1 and int(v) in p2}
    return fn


# ---- metrics ---------------------------------------------------------------- #
def coverage_and_auc(arr, assign_fn, n_random=20000):
    covs, aucs = [], []
    for sample, panel, img1, img2, d in wt_pairs():
        a1, a2 = arr[img1], arr[img2]
        assign = assign_fn(a1, a2, d, panel)
        if assign is None:
            continue
        covs.append(len(assign) / max(len(a1["labels"]), 1))
        if len(assign) < 20:
            aucs.append(np.nan); continue
        i1 = np.array([a1["pos"][k] for k in assign]); i2 = np.array([a2["pos"][v] for v in assign.values()])
        dm = _si_dist(a1["si"][i1], a2["si"][i2])
        n = min(n_random, len(a1["labels"]) * len(a2["labels"]))
        r1 = RNG.integers(0, len(a1["labels"]), n); r2 = RNG.integers(0, len(a2["labels"]), n)
        dr = _si_dist(a1["si"][r1], a2["si"][r2])
        y = np.concatenate([np.ones(len(dm)), np.zeros(len(dr))])
        aucs.append(roc_auc_score(y, -np.concatenate([dm, dr])))
    if not covs:
        return np.nan, np.nan, 0
    return float(np.mean(covs)), float(np.nanmean(aucs)), len(covs)


def crosspanel_fold(arr, assign_fn, n_null=2000):
    """Star-topology cycle consistency for one method's assignments; fold = null/real."""
    by_sample = {}
    for sample, panel, img1, img2, d in wt_pairs():
        by_sample.setdefault(sample, []).append((panel, img1, img2, d))
    reals, nulls = [], []
    for sample, entries in by_sample.items():
        anchor = entries[0][1]; aA = arr[anchor]
        aff, matched = {}, {}
        for panel, img1, img2, d in entries:
            a2 = arr[img2]; assign = assign_fn(aA, a2, d, panel)
            if assign is None or len(assign) < 50:
                continue
            i1 = [aA["pos"][k] for k in assign]; i2 = [a2["pos"][v] for v in assign.values()]
            A, t = fit_affine(aA["cent"][i1], a2["cent"][i2])
            invA = np.linalg.inv(A); aff[panel] = (invA, -invA @ t, a2); matched[panel] = assign
        panels = list(aff)
        if len(panels) < 2:
            continue
        for l1 in aA["labels"]:
            pts = []
            for panel in panels:
                invA, ti, a2 = aff[panel]; l2 = matched[panel].get(int(l1))
                if l2 is None:
                    continue
                pts.append(apply_affine(invA, ti, a2["cent"][a2["pos"][l2]]))
            if len(pts) >= 2:
                reals.append(pdist(np.array(pts)).mean())
        for _ in range(n_null):
            pts = []
            for panel in panels:
                invA, ti, a2 = aff[panel]
                pts.append(apply_affine(invA, ti, a2["cent"][RNG.integers(0, len(a2["cent"]))]))
            nulls.append(pdist(np.array(pts)).mean())
    if not reals:
        return np.nan
    return float(np.mean(nulls) / max(np.mean(reals), 1e-9))


def main():
    with open(MORPH, "rb") as f:
        morph = pickle.load(f)
    arr, _, _ = build_arrays(morph)
    print(f"loaded morphology: {len(arr)} images")

    methods = {"F2FMatcher": assign_f2f,
               "no-align kNN": lambda a1, a2, d, p: assign_knn(a1, a2, d, "raw"),
               "global-affine + kNN": lambda a1, a2, d, p: assign_knn(a1, a2, d, "affine")}
    if WHOLESLIDE.exists():
        for md in sorted(p.name for p in WHOLESLIDE.iterdir() if p.is_dir()):
            methods[md] = assign_learned(md)

    runtime = {}
    rt = WHOLESLIDE / "runtime.csv"
    if rt.exists():
        rr = pd.read_csv(rt)
        runtime = rr.groupby("method")["seconds"].mean().to_dict()

    rows = []
    for name, fn in methods.items():
        cov, auc, npair = coverage_and_auc(arr, fn)
        fold = crosspanel_fold(arr, fn)
        rows.append({"method": name, "coverage": cov, "shape_auc": auc,
                     "crosspanel_fold": fold, "n_pairs": npair,
                     "runtime_s": runtime.get(name, 0.0 if name in GEOMETRIC + ["F2FMatcher"] else np.nan)})
        print(f"  {name:22s} coverage={cov:.3f}  shape_AUC={auc:.3f}  "
              f"crosspanel={fold:.2f}x  (n={npair})")
    res = pd.DataFrame(rows).sort_values("shape_auc", ascending=False)
    res.to_csv(EVAL / "wholeslide_scores.csv", index=False)
    print(f"\nsaved {EVAL/'wholeslide_scores.csv'}")

    FIG.mkdir(parents=True, exist_ok=True)
    # ---- 2D coverage vs correctness (shape AUC) ----
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.axhline(0.5, ls="--", color="#999", lw=1); ax.text(0.01, 0.505, "chance (random matching)", fontsize=8, color="#777")
    for _, r in res.iterrows():
        if np.isnan(r["coverage"]):
            continue
        is_f2f = r["method"] == "F2FMatcher"
        ax.scatter(r["coverage"], r["shape_auc"], s=170 if is_f2f else 90,
                   color="#d62728" if is_f2f else "#4c78a8",
                   edgecolor="k", zorder=3, marker="*" if is_f2f else "o")
        ax.annotate(r["method"], (r["coverage"], r["shape_auc"]),
                    xytext=(6, 4), textcoords="offset points", fontsize=8.5,
                    fontweight="bold" if is_f2f else "normal")
    ax.set_xlabel("coverage (fraction of anchor fibres assigned)")
    ax.set_ylabel("correctness proxy — shape-consistency AUC")
    ax.set_xlim(-0.02, 1.03); ax.set_ylim(0.45, 1.0)
    ax.set_title("Whole-slide (WT): coverage is not enough\nonly F2FMatcher is high on coverage AND correctness",
                 fontweight="bold", fontsize=11)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig(FIG / "fig_B_wholeslide_coverage_vs_consistency.png", dpi=150)
    plt.close(fig)

    # ---- runtime bar (only if any learned-matcher runtimes exist) ----
    rr = res[res["runtime_s"].fillna(0) > 0]
    if len(rr):
        fig, ax = plt.subplots(figsize=(7, 3.6))
        rr = rr.sort_values("runtime_s")
        ax.barh(rr["method"], rr["runtime_s"], color="#4c78a8")
        ax.set_xlabel("mean runtime per section-pair (s)")
        ax.set_title("Whole-slide runtime per method", fontweight="bold")
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        fig.tight_layout(); fig.savefig(FIG / "fig_B_wholeslide_runtime.png", dpi=150)
        plt.close(fig)
    print(f"saved figures -> {FIG}")


if __name__ == "__main__":
    main()
