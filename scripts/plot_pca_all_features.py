import os, sys, pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from sklearn.decomposition import PCA

sys.path.append("/DATA/F2FMatcher_DDC")
from config.ddc_config import *

MUSCLE = "QUA"
DIR_VISUALIZATIONS = Path("/DATA/F2FMatcher_DDC/visualizations")
DIR_COMBINED = Path(f"/DATA/F2FMatcher_DDC/results/{MUSCLE}/features_combined")

COMPARTMENT_ORDER = ["whole", "mem", "cyto1", "cyto2"]
GROUP_COLORS = {"WT": "tab:blue", "mdx": "tab:red", "AAV9": "tab:green", "LICA1": "tab:orange"}
MIN_SLIDES = 4  # min slides that fibers detected


def feature_columns():
    """Column names of the combined vector: 15 mask + 36 per channel, in stored order."""
    cols = list(MASK_FEATURES)
    for slide in SLIDES:
        for ch, stain in SLIDES[slide]["stainings"].items():
            for comp in COMPARTMENT_ORDER:
                for stat in FEATURE_STATISTICS:
                    cols.append(f"S{slide}_{stain}_{comp}_{stat}")
    return cols


def slide_channel_blocks():
    """slide -> (first channel index, n channels) in the 22-channel combined vector."""
    blocks, offset = {}, 0
    for slide in SLIDES:
        n = len(SLIDES[slide]["stainings"])
        blocks[slide] = (offset, n)
        offset += n
    return blocks


def load_all():
    frames = []
    for sample in QUA_SAMPLES:
        with open(DIR_COMBINED / f"{sample}.pkl", "rb") as f:
            combined = pickle.load(f)
        for label_id, vec in combined.items():
            frames.append(list(vec) + [label_id, sample, get_group(sample)])
    cols = feature_columns() + ["label_id_x", "sample", "group"]
    return pd.DataFrame(frames, columns=cols)


def n_slides_detected(X, blocks):
    """Per fiber: number of slides with at least one observed (non-NaN) value."""
    n = np.zeros(X.shape[0], dtype=int)
    for slide in SLIDES:
        st, nc = blocks[slide]
        block = X[:, st * 36:(st + nc) * 36]
        n += (~np.isnan(block)).any(axis=1)
    return n


def group_median_impute(X, groups):
    """Impute NaNs with the median of the same group (fallback: overall median)."""
    X = X.copy()
    overall = np.nanmedian(X, axis=0)
    overall[np.isnan(overall)] = 0.0
    n_imp = 0
    for g in np.unique(groups):
        med = np.nanmedian(X[groups == g], axis=0)
        bad = np.isnan(med)
        med[bad] = overall[bad]
        mask = np.isnan(X[groups == g])
        n_imp += mask.sum()
        X[groups == g] = np.where(mask, med, X[groups == g])
    return X, n_imp


def run_pca(X, tag):
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    pca = PCA(n_components=2, random_state=ANALYSIS_CONFIG["random_state"])
    scores = pca.fit_transform(Xs)
    print(f"[{tag}] explained variance: PC1 {pca.explained_variance_ratio_[0]:.1%}, "
          f"PC2 {pca.explained_variance_ratio_[1]:.1%}")
    return scores


def plot_combined(scores, groups, tag):
    group_order = [g for g in ["WT", "mdx", "AAV9", "LICA1"] if g in set(groups)]
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_axes([0.08, 0.12, 0.78, 0.78])
    for g in group_order:
        idx = np.where(groups == g)[0]
        ax.scatter(scores[idx, 0], scores[idx, 1], s=8, alpha=0.3,
                   color=GROUP_COLORS.get(g, "gray"), label=f"{g} (n={len(idx)})", edgecolors="none")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(title="Group")
    fig.suptitle(f"{MUSCLE} — PCA all features [{tag}]")

    ax_top = fig.add_axes([0.08, 0.92, 0.78, 0.06])
    ax_top.hist(scores[:, 0], bins=60, density=True, histtype="step", color="black")
    ax_top.set_xticks([])
    ax_top.set_title("PC1", fontsize=9)
    ax_right = fig.add_axes([0.88, 0.12, 0.06, 0.78])
    ax_right.hist(scores[:, 1], bins=60, density=True, histtype="step", color="black")
    ax_right.set_yticks([])
    ax_right.set_title("PC2", fontsize=9)

    out = DIR_VISUALIZATIONS / f"{MUSCLE}_pca_all_{tag}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out.name}")


def plot_per_group(scores, groups, tag):
    group_order = [g for g in ["WT", "mdx", "AAV9", "LICA1"] if g in set(groups)]
    fig2, axes = plt.subplots(2, 2, figsize=(12, 11), sharex=True, sharey=True)
    for axi, g in zip(axes.flat, group_order):
        idx = np.where(groups == g)[0]
        axi.scatter(scores[idx, 0], scores[idx, 1], s=8, alpha=0.3,
                    color=GROUP_COLORS.get(g, "gray"), edgecolors="none")
        axi.set_title(f"{g} (n={len(idx)})", fontsize=12)
        axi.set_xlabel("PC1")
        axi.set_ylabel("PC2")
    fig2.suptitle(f"{MUSCLE} — PCA all features [{tag}], per group")
    fig2.tight_layout(rect=[0, 0, 1, 0.96])
    out2 = DIR_VISUALIZATIONS / f"{MUSCLE}_pca_all_{tag}_groups.png"
    fig2.savefig(out2, dpi=200)
    plt.close(fig2)
    print(f"  saved {out2.name}")


def group_stats(scores, groups, tag):
    group_order = [g for g in ["WT", "mdx", "AAV9", "LICA1"] if g in set(groups)]
    print(f"\n[{tag}] per-group PC statistics:")
    for pc in (0, 1):
        groups_vals = [scores[groups == g, pc] for g in group_order]
        h, p = kruskal(*groups_vals)
        print(f"PC{pc+1}: Kruskal-Wallis H={h:.1f}, p={p:.2e}")
        for g in group_order:
            v = scores[groups == g, pc]
            print(f"  {g:6s} mean {v.mean():8.3f}  std {v.std():7.3f}  n {len(v)}")


def main():
    DIR_VISUALIZATIONS.mkdir(parents=True, exist_ok=True)
    df = load_all()
    feat_cols = feature_columns()
    blocks = slide_channel_blocks()
    X = df[feat_cols].to_numpy(dtype=np.float64)
    groups = df["group"].to_numpy()
    print(f"Total fibers: {len(df)} ({len(feat_cols)} features)")

    # filter: fibers detected in at least MIN_SLIDES of the 7 slides
    n_slides = n_slides_detected(X, blocks)
    keep = n_slides >= MIN_SLIDES
    print(f"Slide-detection distribution: {pd.Series(n_slides).value_counts().sort_index().to_dict()}")
    print(f"Keeping {keep.sum()} fibers detected in >= {MIN_SLIDES}/7 slides "
          f"({keep.mean():.1%})")
    X, groups, samples = X[keep], groups[keep], df["sample"].to_numpy()[keep]

    # imputation: group-conditional medians (KNN is O(n^2 d), infeasible at this scale)
    X, n_imp = group_median_impute(X, groups)
    print(f"Group-conditional median imputation: {n_imp} values ({n_imp / X.size:.2%})")

    scores = run_pca(X, "grp_median_4of7")
    plot_combined(scores, groups, "grp_median_4of7")
    plot_per_group(scores, groups, "grp_median_4of7")
    group_stats(scores, groups, "grp_median_4of7")

    out_csv = DIR_VISUALIZATIONS / f"{MUSCLE}_pca_all_grp_median_4of7_scores.csv"
    pd.DataFrame({
        "sample": samples, "group": groups,
        "PC1": scores[:, 0], "PC2": scores[:, 1],
    }).to_csv(out_csv, index=False)
    print(f"\nSaved scores to {out_csv}")


if __name__ == "__main__":
    main()
