"""Accumulate the FINAL complete §B comparison (no re-runs) from the two existing
result CSVs, and render one figure.

Sources:
  - results/benchmark/benchmark_scores.csv   : F2FMatcher + geometric baselines (edge-corrected)
  - results/benchmark/vismatch/vismatch_scores.csv : DINOv2 + VisMatch (loftr/roma/superpoint/superglue)

Recall/F1 use the matchable (non-edge) GT denominator (the fair set) for the
F2FMatcher/geometric side; the dense/sparse matchers are scored on the full GT
(they can map any pixel, so no edge discard applies to them).

Outputs:
  - results/benchmark/final_comparison.csv
  - visualizations/eval/fig_B4_final_comparison.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path("/DATA/F2FMatcher_DDC")
BENCH = REPO / "results" / "benchmark"
VISM = BENCH / "vismatch" / "vismatch_scores.csv"
GEOM = BENCH / "benchmark_scores.csv"
OUT_CSV = BENCH / "final_comparison.csv"
OUT_FIG = REPO / "visualizations" / "eval" / "fig_B4_final_comparison.png"
OUT_FIG.parent.mkdir(parents=True, exist_ok=True)


def prf(tp, np_, ng):
    p = tp / np_ if np_ else 0.0
    r = tp / ng if ng else 0.0
    f1 = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f1


def main():
    rows = []

    # --- F2FMatcher + geometric baselines (edge-corrected, matchable GT) ---
    g = pd.read_csv(GEOM)
    for m in ["F2FMatcher", "no-align kNN", "affine+kNN", "random"]:
        s = g[g.method == m]
        if m == "F2FMatcher":
            s = s[s.f2f_ran]
        tp, np_, ng = s.tp.sum(), s.n_pred.sum(), s.n_gt_matchable.sum()
        p, r, f1 = prf(tp, np_, ng)
        rows.append({"method": m, "family": "F2FMatcher" if m == "F2FMatcher" else "geometric",
                     "pairs": f"{len(s)}/38", "P": p, "R": r, "F1": f1,
                     "denom": "matchable GT"})

    # --- DINOv2 + VisMatch (full GT; dense/sparse matchers can map any pixel) ---
    v = pd.read_csv(VISM)
    name_map = {"dinov2": "DINOv2 (ViT-S/14)", "superpoint-lightglue": "SuperPoint-LightGlue",
                "superglue": "SuperGlue", "loftr": "LoFTR", "roma": "RoMa"}
    fam_map = {"dinov2": "appearance", "superpoint-lightglue": "sparse",
               "superglue": "sparse", "loftr": "semi-dense", "roma": "dense"}
    for m in ["dinov2", "superpoint-lightglue", "superglue", "loftr", "roma"]:
        s = v[v.method == m]
        tp, np_, ng = s.tp.sum(), s.n_pred.sum(), s.n_gt.sum()
        p, r, f1 = prf(tp, np_, ng)
        rows.append({"method": name_map[m], "family": fam_map[m], "pairs": "38/38",
                     "P": p, "R": r, "F1": f1, "denom": "full GT"})

    res = pd.DataFrame(rows)
    res = res.sort_values("F1", ascending=False).reset_index(drop=True)
    res.to_csv(OUT_CSV, index=False)
    print("\n=== FINAL §B comparison (sorted by F1) ===")
    print(res.to_string(index=False,
          formatters={"P": "{:.3f}".format, "R": "{:.3f}".format, "F1": "{:.3f}".format}))

    # --- figure: grouped bar chart of P / R / F1 ---
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(res))
    w = 0.26
    for i, (col, c) in enumerate([("P", "#4C72B0"), ("R", "#DD8452"), ("F1", "#55A868")]):
        ax.bar(x + (i - 1) * w, res[col], w, label=col, color=c)
    ax.set_xticks(x)
    ax.set_xticklabels(res["method"], rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("score")
    ax.set_ylim(0, 1.05)
    ax.set_title("F2FMatcher vs baselines — curated GT (TA crops)\n"
                 "P/R/F1; F2FMatcher+geometric on matchable GT, dense/sparse on full GT")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)
    print(f"\nsaved {OUT_CSV}\nsaved {OUT_FIG}")


if __name__ == "__main__":
    main()
