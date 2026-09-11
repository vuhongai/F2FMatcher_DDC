"""Regenerate UMAP plots from the precomputed cache:
1) colored by group, 2) colored by GMM clusters (Cls1..Clsk),
3) 2x2 per-group subplots colored by cluster,
4) optional 1x4 per-group subplots colored by --plot-feature.
No PCA/UMAP re-run — reads the cache (built by cluster_fibers.py).
"""
import sys, argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.mixture import GaussianMixture

sys.path.insert(0, "/DATA/F2FMatcher_DDC")
sys.path.insert(0, "/DATA/F2FMatcher_DDC/scripts")
from config.ddc_config import *
from plot_pca_all_features import load_all, feature_columns, group_median_impute
from cluster_fibers import (
    load_cache, set_preset, resolve_feature,
    n_included_slides_detected_full, min_slides_effective,
    MUSCLE, DIR_VISUALIZATIONS, RANDOM_STATE,
)

GROUPS = ["WT", "mdx", "AAV9", "LICA1"]
GROUP_COLORS = {
    "WT": "#0000C0",
    "mdx": "#FF6000",
    "AAV9": "#C0C000",
    "LICA1": "#008000",
}
CLUSTER_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                  "#a65628", "#f781bf", "#66c2a5"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--norm", default="zscore")
    parser.add_argument("--preset", default="ch3")
    parser.add_argument("--n-pca", type=int, default=30)
    parser.add_argument("--mindist", type=float, default=0.0)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--k", type=int, default=6)
    parser.add_argument("--plot-feature", default=None,
                        help="extra 1x4 per-group plot colored by this feature "
                             "(raw column or marker name, e.g. area, S1_Dystrophin_mem_mean, COX)")
    parser.add_argument("--vmin", type=float, default=None,
                        help="cap color scale low value (default: data min)")
    parser.add_argument("--vmax", type=float, default=None,
                        help="cap color scale high value (default: data max)")
    parser.add_argument("--no-morph", action="store_true")
    args = parser.parse_args()

    feat_tag = set_preset(args.preset) + ("_nomorph" if args.no_morph else "")
    DIR_VISUALIZATIONS.mkdir(parents=True, exist_ok=True)
    data = load_cache(args.norm, args.n_pca, args.mindist, feat_tag, args.n_neighbors)
    X_pca, X_umap = data["X_pca"], data["X_umap"]
    groups = data["groups"]
    n = len(X_umap)
    print(f"{n:,} fibers, PCA{args.n_pca}, UMAP min_dist={args.mindist}, nn={args.n_neighbors}")

    fig, ax = plt.subplots(figsize=(6, 5))
    handles = []
    for g in GROUPS:
        m = groups == g
        ax.scatter(X_umap[m, 0], X_umap[m, 1], s=6, alpha=0.25,
                   color=GROUP_COLORS[g], edgecolors="none")
        handles.append(Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=GROUP_COLORS[g], markersize=8,
                              alpha=1.0, label=f"{g}"))
    ax.legend(handles=handles, fontsize=10, loc="best")
    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(f"UMAP (Morphology+HPS)", fontsize=14, fontweight="bold")
    ax.tick_params(labelsize=10)
    fig.tight_layout()
    out = DIR_VISUALIZATIONS / f"{MUSCLE}_umap_{args.norm}_{feat_tag}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Saved {out.name}")
    for g in GROUPS:
        print(f"  {g:6s} {int((groups == g).sum()):>7,}")

    # ---- Plot 2: UMAP colored by GMM clusters (Cls1..Clsk) ----
    print(f"GMM k={args.k} (labeling only) ...")
    gm = GaussianMixture(n_components=args.k, covariance_type="full",
                         n_init=1, max_iter=200, random_state=RANDOM_STATE)
    gm.fit(X_pca)
    cluster_ids = gm.predict(X_pca)

    fig2, ax2 = plt.subplots(figsize=(6, 5))
    handles2 = []
    for c in range(args.k):
        m = cluster_ids == c
        col = CLUSTER_COLORS[c % len(CLUSTER_COLORS)]
        ax2.scatter(X_umap[m, 0], X_umap[m, 1], s=6, alpha=0.25,
                    color=col, edgecolors="none")
        handles2.append(Line2D([0], [0], marker="o", color="w",
                               markerfacecolor=col, markersize=8,
                               alpha=1.0, label=f"Cls{c + 1}"))
    ax2.legend(handles=handles2, fontsize=10, loc="best")
    ax2.set_xlabel("UMAP 1", fontsize=12)
    ax2.set_ylabel("UMAP 2", fontsize=12)
    ax2.set_title(f"UMAP (Morphology+HPS)", fontsize=14,
                  fontweight="bold")
    ax2.tick_params(labelsize=10)
    fig2.tight_layout()
    out2 = DIR_VISUALIZATIONS / f"{MUSCLE}_umap_clusters_{args.norm}_{feat_tag}_k{args.k}.png"
    fig2.savefig(out2, dpi=300)
    plt.close(fig2)
    print(f"Saved {out2.name}")
    for c in range(args.k):
        m = cluster_ids == c
        fr = "  ".join(f"{g} {(groups[m] == g).mean():.0%}" for g in GROUPS)
        print(f"  Cls{c + 1}: n={int(m.sum()):>7,}  {fr}")

    # ---- Plot 3: 2x2 per-group subplots, colored by cluster, minimal ----
    fig3, axes3 = plt.subplots(2, 2, figsize=(7, 5), sharex=True, sharey=True)
    for ax, g in zip(axes3.ravel(), GROUPS):
        m = groups == g
        for c in range(args.k):
            mc = m & (cluster_ids == c)
            ax.scatter(X_umap[mc, 0], X_umap[mc, 1], s=6, alpha=0.25,
                       color=CLUSTER_COLORS[c % len(CLUSTER_COLORS)], edgecolors="none")
        ax.text(0.05, 0.95, g, transform=ax.transAxes, fontsize=12, fontweight="bold",
                ha="left", va="top", color="black")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig3.tight_layout()
    out3 = DIR_VISUALIZATIONS / f"{MUSCLE}_umap_groups_clusters_{args.norm}_{feat_tag}_k{args.k}.png"
    fig3.savefig(out3, dpi=300)
    plt.close(fig3)
    print(f"Saved {out3.name}")

    # ---- Plot 4 (optional): 1x4 per-group subplots colored by feature ----
    if args.plot_feature:
        df = load_all()
        full_cols = feature_columns()
        Xfull = df[full_cols].to_numpy(dtype=np.float64)
        keep = n_included_slides_detected_full(Xfull) >= min_slides_effective()
        Xfull = Xfull[keep]
        Xfull, _ = group_median_impute(Xfull, groups)
        vals = resolve_feature(args.plot_feature, Xfull, full_cols)
        vmin = args.vmin if args.vmin is not None else float(vals.min())
        vmax = args.vmax if args.vmax is not None else float(vals.max())

        fig4, axes4 = plt.subplots(1, 4, figsize=(16, 4), sharex=True, sharey=True)
        for ax, g in zip(axes4, GROUPS):
            m = groups == g
            ax.scatter(X_umap[m, 0], X_umap[m, 1], c=vals[m], cmap="viridis",
                       vmin=vmin, vmax=vmax, s=6, alpha=0.25, edgecolors="none")
            ax.text(0.05, 0.95, g, transform=ax.transAxes, fontsize=12, fontweight="bold",
                    ha="left", va="top", color="black")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
        # colorbar outside, right of the row, no title
        cb = fig4.colorbar(plt.cm.ScalarMappable(cmap="viridis",
                                                 norm=plt.Normalize(vmin, vmax)),
                           ax=axes4[-1], fraction=0.03, pad=0.06)
        cb.ax.tick_params(labelsize=8)
        fig4.tight_layout()
        cap = ""
        if args.vmin is not None or args.vmax is not None:
            lo = f"{args.vmin:g}" if args.vmin is not None else "min"
            hi = f"{args.vmax:g}" if args.vmax is not None else "max"
            cap = f"_cap{lo}_{hi}"
        out4 = DIR_VISUALIZATIONS / f"{MUSCLE}_umap_feature_{args.plot_feature}{cap}_{args.norm}_{feat_tag}.png"
        fig4.savefig(out4, dpi=300)
        plt.close(fig4)
        print(f"Saved {out4.name}")


if __name__ == "__main__":
    main()
