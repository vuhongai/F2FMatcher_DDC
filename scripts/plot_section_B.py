"""Plot the Section B evaluation figures (WT cohort).

Reads the CSV/JSON written by scripts/evaluate_section_B.py and writes
publication-style PNGs to visualizations/eval/:
    fig_B1_coverage.png          per-panel coverage
    fig_B2_morph_consistency.png matched-vs-random shape consistency + AUC/panel
    fig_B2b_crosspanel.png       cross-panel spatial consistency (real vs null)
    fig_B3_baselines.png         F2FMatcher vs geometric baselines

Usage:  python scripts/plot_section_B.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path("/DATA/F2FMatcher_DDC")
EVAL = REPO / "results" / "QUA" / "eval"
FIG = REPO / "visualizations" / "eval"
FIG.mkdir(parents=True, exist_ok=True)

PANEL_ORDER = ["10X_NADH", "10X_HE", "10X_DAPI_LAM_IgG_CD11B", "10X_COX",
               "10X_WGA_Myh7_Myh2_Myh4", "10X_DAPI_LAMP2_LGALS3_SQSTM1"]
PANEL_SHORT = {"10X_NADH": "NADH", "10X_HE": "HE",
               "10X_DAPI_LAM_IgG_CD11B": "IgG/CD11b", "10X_COX": "COX",
               "10X_WGA_Myh7_Myh2_Myh4": "WGA/Myh",
               "10X_DAPI_LAMP2_LGALS3_SQSTM1": "LAMP2/LGALS3/SQSTM1"}
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 200, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})


def _order(df, col="panel"):
    df = df.copy()
    df["_o"] = df[col].map({p: i for i, p in enumerate(PANEL_ORDER)})
    return df.sort_values("_o")


def fig_b1():
    bp = _order(pd.read_csv(EVAL / "coverage_WT_by_panel.csv"))
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    x = np.arange(len(bp))
    ax.bar(x, bp["cov_mean"] * 100, color="#4C72B0", width=0.62)
    for xi, (c, npair) in enumerate(zip(bp["cov_mean"], bp["mean_pairs"])):
        ax.text(xi, c * 100 + 1.5, f"{c*100:.0f}%", ha="center", fontsize=8)
        ax.text(xi, 3, f"n={int(npair)}", ha="center", color="white", fontsize=7)
    ax.set_xticks(x, [PANEL_SHORT[p] for p in bp["panel"]], rotation=18, ha="right")
    ax.set_ylabel("mean coverage (both sides), %")
    ax.set_ylim(0, 100)
    ax.set_title("B1 — F2FMatcher mapping coverage per panel (WT, n=5 samples/panel)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_B1_coverage.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_B1_coverage.png")


def fig_b2():
    mc = _order(pd.read_csv(EVAL / "morph_consistency.csv"))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.6, 3.6))
    # left: per-pair mean matched vs random shape distance (dot plot)
    y = np.arange(len(mc))[::-1]
    a1.scatter(mc["matched_mean"], y, s=14, color="#C44E52", label="matched (F2FMatcher)")
    a1.scatter(mc["random_mean"], y, s=14, color="#888888", label="random null")
    for yi, m, r, auc in zip(y, mc["matched_mean"], mc["random_mean"], mc["auc"]):
        a1.plot([m, r], [yi, yi], color="#cccccc", lw=0.8, zorder=0)
    a1.set_yticks(y, [PANEL_SHORT[p] for p in mc["panel"]])
    a1.set_xlabel("mean shape distance (scale-invariant, std. units)")
    a1.set_title("B2a — shape consistency per panel")
    a1.legend(fontsize=8, loc="lower right")
    # right: AUC per panel
    x = np.arange(len(mc))
    a2.bar(x, mc["auc"], color="#4C72B0", width=0.62)
    a2.axhline(0.5, color="k", ls=":", lw=1)
    a2.text(len(mc) - 0.4, 0.51, "chance", fontsize=7, ha="right")
    for xi, v in zip(x, mc["auc"]):
        a2.text(xi, v + 0.005, f"{v:.2f}", ha="center", fontsize=7)
    a2.set_xticks(x, [PANEL_SHORT[p] for p in mc["panel"]], rotation=18, ha="right")
    a2.set_ylabel("AUC (matched vs random)")
    a2.set_ylim(0, 1.0)
    a2.set_title("B2a — AUC of shape-consistency separation")
    fig.suptitle("Morphological consistency: F2FMatcher matches vs random (WT)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_B2_morph_consistency.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_B2_morph_consistency.png")


def fig_b2b():
    cp = pd.read_csv(EVAL / "crosspanel_consistency.csv")
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    x = np.arange(len(cp))
    w = 0.38
    ax.bar(x - w / 2, cp["real_spread_mean"], w, color="#4C72B0", label="F2FMatcher matches")
    ax.bar(x + w / 2, cp["null_spread_mean"], w, color="#bbbbbb", label="random null")
    for xi, r, n, f in zip(x, cp["real_spread_mean"], cp["null_spread_mean"],
                           cp["fold_reduction"]):
        ax.text(xi + w / 2, n + 60, f"{f:.1f}×", ha="center", fontsize=8, color="#333")
    ax.set_xticks(x, cp["sample"])
    ax.set_ylabel("cross-panel spread (px)")
    ax.set_title("B2b — cross-panel spatial consistency (star-topology cycle check)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_B2b_crosspanel.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_B2b_crosspanel.png")


def fig_b3():
    cmp = pd.read_csv(EVAL / "baseline_compare.csv")
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    colors = ["#C44E52", "#4C72B0", "#8172B2", "#bbbbbb"]
    x = np.arange(len(cmp))
    ax.bar(x, cmp["mean_morph_dist"], color=colors, width=0.6)
    for xi, v, r in zip(x, cmp["mean_morph_dist"], cmp["reduction_vs_random"]):
        ax.text(xi, v + 0.05, f"{v:.2f}\n({r*100:+.0f}% vs random)",
                ha="center", fontsize=8)
    ax.axhline(cmp.iloc[3]["mean_morph_dist"], color="k", ls=":", lw=1)
    ax.set_xticks(x, ["F2FMatcher", "global-affine\n+ kNN", "no-align\nkNN", "random\nnull"])
    ax.set_ylabel("mean shape distance (lower = more consistent)")
    ax.set_ylim(0, cmp["mean_morph_dist"].max() * 1.25)
    ax.set_title("B3 — F2FMatcher vs geometric baselines (shape consistency, WT)")
    fig.tight_layout()
    fig.savefig(FIG / "fig_B3_baselines.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_B3_baselines.png")


if __name__ == "__main__":
    fig_b1()
    fig_b2()
    fig_b2b()
    fig_b3()
    print("done")
