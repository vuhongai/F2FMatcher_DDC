import sys, argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from sklearn.manifold import TSNE

sys.path.insert(0, "/DATA/F2FMatcher_DDC")
sys.path.insert(0, "/DATA/F2FMatcher_DDC/scripts")
from config.ddc_config import *
from cluster_fibers import (
    load_cache, load_all, n_included_slides_detected, group_median_impute,
    MIN_SLIDES, MUSCLE, DIR_VISUALIZATIONS, GROUP_COLORS, KEY_FEATURES, RANDOM_STATE,
)

try:
    from matplotlib import colormaps
    def get_cmap(name):
        return colormaps[name]
except ImportError:
    def get_cmap(name):
        return plt.cm.get_cmap(name)

CLUSTER_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                  "#a65628", "#f781bf", "#66c2a5"]
GROUPS = ["WT", "mdx", "AAV9", "LICA1"]
DOT_S = 14
DOT_ALPHA = 0.35


def scatter_cluster(ax, xy, cluster_ids, code, order, n_per_cluster, show_legend=False):
    for t, c in enumerate(order):
        m = code == t
        ax.scatter(xy[m, 0], xy[m, 1], s=DOT_S, alpha=DOT_ALPHA,
                   color=CLUSTER_COLORS[t % len(CLUSTER_COLORS)],
                   label=f"{c} (n={n_per_cluster[c]:,})", edgecolors="none")
    if show_legend:
        ax.legend(title="cluster", fontsize=9, title_fontsize=10,
                  framealpha=0.9, loc="upper right")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.tick_params(labelsize=9)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--norm", default="zscore")
    parser.add_argument("--n-pca", type=int, default=30)
    parser.add_argument("--mindist", type=float, default=0.0)
    parser.add_argument("--k", type=int, default=6)
    parser.add_argument("--perplexity", type=int, default=30)
    args = parser.parse_args()

    DIR_VISUALIZATIONS.mkdir(parents=True, exist_ok=True)
    data = load_cache(args.norm, args.n_pca, args.mindist)
    X_pca, X_umap = data["X_pca"], data["X_umap"]
    groups, samples, feat_cols = data["groups"], data["samples"], data["feat_cols"]
    n = len(X_pca)
    print(f"{n:,} fibers, PCA{args.n_pca}, UMAP min_dist={args.mindist}")

    print(f"t-SNE (perplexity={args.perplexity}) ...")
    tsne = TSNE(n_components=2, perplexity=args.perplexity, init="pca",
                learning_rate="auto", n_iter=1000, random_state=RANDOM_STATE, n_jobs=8)
    X_tsne = tsne.fit_transform(X_pca)

    print(f"GMM k={args.k} ...")
    gm = GaussianMixture(n_components=args.k, covariance_type="full",
                         n_init=1, max_iter=200, random_state=RANDOM_STATE)
    gm.fit(X_pca)
    cluster_ids = gm.predict(X_pca)

    sizes = pd.Series(cluster_ids).value_counts()
    order = list(sizes.index)
    code = np.array([order.index(c) for c in cluster_ids])
    n_per_cluster = {c: int(sizes[c]) for c in order}

    # raw (imputed) features for the feature-coloring plot
    df = load_all()
    Xraw = df[feat_cols].to_numpy(dtype=np.float64)
    keep = n_included_slides_detected(Xraw) >= MIN_SLIDES
    Xraw = Xraw[keep]
    Xraw, _ = group_median_impute(Xraw, groups)

    tag = f"{args.norm}_k{args.k}"

    # ---- Figure 1: UMAP | tSNE side by side, colored by cluster ----
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    scatter_cluster(axes[0], X_umap, cluster_ids, code, order, n_per_cluster, show_legend=True)
    axes[0].set_title(f"UMAP (min_dist={args.mindist})", fontsize=13)
    scatter_cluster(axes[1], X_tsne, cluster_ids, code, order, n_per_cluster)
    axes[1].set_title(f"t-SNE (perplexity={args.perplexity})", fontsize=13)
    fig.suptitle(f"{MUSCLE} — GMM k={args.k} clusters on UMAP vs t-SNE "
                 f"(norm={args.norm}, PCA{args.n_pca}, {n:,} fibers)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out1 = DIR_VISUALIZATIONS / f"{MUSCLE}_umap_vs_tsne_{tag}.png"
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"Saved {out1.name}")

    # ---- Figure 2: 4 groups x 2 embeddings, colored by cluster ----
    fig2, axes2 = plt.subplots(4, 2, figsize=(14, 22), sharex="col", sharey="col")
    for r, g in enumerate(GROUPS):
        m_g = groups == g
        scatter_cluster(axes2[r, 0], X_umap[m_g], cluster_ids, code[m_g], order, n_per_cluster)
        axes2[r, 0].set_title(f"{g} — UMAP (n={int(m_g.sum()):,})",
                              fontsize=11, color=GROUP_COLORS[g])
        scatter_cluster(axes2[r, 1], X_tsne[m_g], cluster_ids, code[m_g], order, n_per_cluster)
        axes2[r, 1].set_title(f"{g} — t-SNE (n={int(m_g.sum()):,})",
                              fontsize=11, color=GROUP_COLORS[g])
    axes2[0, 0].legend(title="cluster", fontsize=8, title_fontsize=9,
                       framealpha=0.9, loc="upper right")
    for ax in axes2[-1]:
        ax.set_xlabel("x")
    for ax in axes2[:, 0]:
        ax.set_ylabel("UMAP y")
    for ax in axes2[:, 1]:
        ax.set_ylabel("t-SNE y")
    fig2.suptitle(f"{MUSCLE} — per-group cluster composition (UMAP | t-SNE, k={args.k})",
                  fontsize=14)
    fig2.tight_layout(rect=[0, 0, 1, 0.98])
    out2 = DIR_VISUALIZATIONS / f"{MUSCLE}_groups_umap_tsne_{tag}.png"
    fig2.savefig(out2, dpi=150)
    plt.close(fig2)
    print(f"Saved {out2.name}")

    # ---- Figure 3: raw feature values on UMAP (3x3) ----
    fig3, axes3 = plt.subplots(3, 3, figsize=(18, 15))
    for ax, (name, col) in zip(axes3.flat, KEY_FEATURES.items()):
        vals = np.log2(Xraw[:, feat_cols.index(col)] + 1)
        sc = ax.scatter(X_umap[:, 0], X_umap[:, 1], c=vals, cmap="viridis",
                        s=DOT_S, alpha=DOT_ALPHA, edgecolors="none")
        cb = plt.colorbar(sc, ax=ax, fraction=0.046)
        cb.set_label(name, fontsize=10)
        ax.set_title(name, fontsize=12)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.tick_params(labelsize=8)
    fig3.suptitle(f"{MUSCLE} — raw feature values on UMAP (k={args.k} clusters for reference)",
                  fontsize=14)
    fig3.tight_layout(rect=[0, 0, 1, 0.97])
    out3 = DIR_VISUALIZATIONS / f"{MUSCLE}_features_umap_{tag}.png"
    fig3.savefig(out3, dpi=150)
    plt.close(fig3)
    print(f"Saved {out3.name}")

    print("\nCluster sizes:")
    for c in order:
        m = cluster_ids == c
        fr = {g: f"{(groups[m] == g).mean():.0%}" for g in GROUPS}
        print(f"  cluster {c}: n={n_per_cluster[c]:,}  " + "  ".join(f"{g} {v}" for g, v in fr.items()))
    print("\nDone.")


if __name__ == "__main__":
    main()
