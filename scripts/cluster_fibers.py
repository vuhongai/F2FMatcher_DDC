import sys, argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import igraph as ig
import leidenalg
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
from umap import UMAP

sys.path.insert(0, "/DATA/F2FMatcher_DDC")
sys.path.insert(0, "/DATA/F2FMatcher_DDC/scripts")
from config.ddc_config import *
from plot_pca_all_features import (
    load_all, feature_columns, slide_channel_blocks,
    n_slides_detected, group_median_impute,
    COMPARTMENT_ORDER, FEATURE_STATISTICS,
    MUSCLE, DIR_VISUALIZATIONS, GROUP_COLORS,
)

try:
    from matplotlib import colormaps
    def get_cmap(name):
        return colormaps[name]
except ImportError:
    def get_cmap(name):
        return plt.cm.get_cmap(name)

N_NEIGHBORS = 15
MIN_SLIDES = 4
RANDOM_STATE = ANALYSIS_CONFIG["random_state"]
CACHE_DIR = Path(f"/DATA/F2FMatcher_DDC/results/{MUSCLE}/clustering_cache")

# Channels used to build UMAP/tSNE (slide -> stainings, None = all).
INCLUDED_SLIDES = {
    1: ["DAPI"],
    2: ["IgG", "CD11b"],
    4: ["Myh7", "Myh2", "Myh4"],
    6: None,
}


def included_channels():
    chs = []
    for slide in SLIDES:
        if slide not in INCLUDED_SLIDES:
            continue
        sel = INCLUDED_SLIDES[slide]
        for stain in SLIDES[slide]["stainings"].values():
            if sel is None or stain in sel:
                chs.append((slide, stain))
    return chs


def included_feature_columns():
    cols = list(MASK_FEATURES)
    for slide, stain in included_channels():
        for comp in COMPARTMENT_ORDER:
            for stat in FEATURE_STATISTICS:
                cols.append(f"S{slide}_{stain}_{comp}_{stat}")
    return cols


def included_slide_blocks():
    """slide -> (col_start, col_end) in the 15-mask + 36-per-channel subset layout."""
    blocks = {}
    for i, (slide, stain) in enumerate(included_channels()):
        lo, hi = 15 + i * 36, 15 + (i + 1) * 36
        if slide in blocks:
            blocks[slide] = (blocks[slide][0], hi)
        else:
            blocks[slide] = (lo, hi)
    return blocks


def n_included_slides_detected(X):
    """Per fiber: number of included slides with at least one observed value."""
    n = np.zeros(X.shape[0], dtype=int)
    for lo, hi in included_slide_blocks().values():
        n += (~np.isnan(X[:, lo:hi])).any(axis=1)
    return n


FEAT_TAG = f"ch{len(included_channels())}"

KEY_FEATURES = {
    "area": "area",
    "DAPI": "S1_DAPI_cyto2_mean",
    "Laminin": "S1_Laminin_mem_mean",
    "Dystrophin": "S1_Dystrophin_mem_mean",
    "IgG": "S2_IgG_whole_mean",
    "CD11b": "S2_CD11b_mem_mean",
    "Myh7": "S4_Myh7_whole_mean",
    "Myh2": "S4_Myh2_whole_mean",
    "Myh4": "S4_Myh4_whole_mean",
}


def normalize(X, mode, feat_cols, groups=None):
    """mode: none | zscore | groupzscore (per-feature, per-group)"""
    if mode == "none":
        return X
    X = X.copy()
    if mode == "groupzscore":
        for g in np.unique(groups):
            mu = np.nanmean(X[groups == g], axis=0)
            sd = np.nanstd(X[groups == g], axis=0)
            sd[sd == 0] = 1.0
            X[groups == g] = (X[groups == g] - mu) / sd
        return X
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def build_cache(norm, n_pca, mindist):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / f"{FEAT_TAG}_{norm}_pca{n_pca}_mindist{mindist}.npz"
    df = load_all()
    feat_cols = included_feature_columns()
    X = df[feat_cols].to_numpy(dtype=np.float64)
    groups = df["group"].to_numpy()
    samples = df["sample"].to_numpy()
    print(f"Total fibers: {len(df)} ({len(feat_cols)} features, {FEAT_TAG})")

    n_incl = len(included_slide_blocks())
    keep = n_included_slides_detected(X) >= MIN_SLIDES
    print(f"Keeping {keep.sum()} fibers detected in >= {MIN_SLIDES}/{n_incl} included slides")
    X, groups, samples = X[keep], groups[keep], samples[keep]

    X, n_imp = group_median_impute(X, groups)
    print(f"Group-conditional median imputation: {n_imp} values ({n_imp / X.size:.2%})")
    X = normalize(X, norm, feat_cols, groups)
    print(f"Normalization: {norm}")

    print(f"\nPCA ({n_pca} components) ...")
    pca = PCA(n_components=n_pca, svd_solver="randomized", random_state=RANDOM_STATE)
    X_pca = pca.fit_transform(X)
    print("  explained variance (top 5):", np.round(pca.explained_variance_ratio_[:5], 3))
    for pc in (0, 1, 2):
        load = np.abs(pca.components_[pc])
        top = np.argsort(load)[::-1][:5]
        print(f"  PC{pc+1} top loadings: " +
              ", ".join(f"{feat_cols[j]}({load[j]:.2f})" for j in top))

    print(f"\nUMAP (n_neighbors={N_NEIGHBORS}, min_dist={mindist}) ...")
    reducer = UMAP(n_neighbors=N_NEIGHBORS, min_dist=mindist, metric="euclidean",
                   random_state=RANDOM_STATE)
    X_umap = reducer.fit_transform(X_pca)

    np.savez_compressed(cache_path, X_pca=X_pca, X_umap=X_umap,
                        samples=samples, groups=groups, feat_cols=np.array(feat_cols))
    print(f"Saved cache to {cache_path}")
    return cache_path


def load_cache(norm, n_pca, mindist):
    cache_path = CACHE_DIR / f"{FEAT_TAG}_{norm}_pca{n_pca}_mindist{mindist}.npz"
    if not cache_path.exists():
        build_cache(norm, n_pca, mindist)
    print(f"Loading cache {cache_path}")
    z = np.load(cache_path, allow_pickle=True)
    return {
        "X_pca": z["X_pca"], "X_umap": z["X_umap"],
        "samples": z["samples"], "groups": z["groups"],
        "feat_cols": [str(c) for c in z["feat_cols"]],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--norm", default="zscore", choices=["none", "zscore", "groupzscore", "user"])
    parser.add_argument("--method", default="gmm", choices=["gmm", "leiden"])
    parser.add_argument("--kmin", type=int, default=6)
    parser.add_argument("--kmax", type=int, default=8)
    parser.add_argument("--resolutions", default="0.05,0.1,0.15,0.2,0.3,0.5,1.0")
    parser.add_argument("--n-pca", type=int, default=30)
    parser.add_argument("--mindist", type=float, default=0.0)
    parser.add_argument("--rebuild-cache", action="store_true")
    args = parser.parse_args()

    DIR_VISUALIZATIONS.mkdir(parents=True, exist_ok=True)
    if args.rebuild_cache:
        build_cache(args.norm, args.n_pca, args.mindist)
        data = load_cache(args.norm, args.n_pca, args.mindist)
    else:
        data = load_cache(args.norm, args.n_pca, args.mindist)
    X_pca, X_umap = data["X_pca"], data["X_umap"]
    groups, samples, feat_cols = data["groups"], data["samples"], data["feat_cols"]

    # raw (imputed) features for annotation
    df = load_all()
    Xraw = df[feat_cols].to_numpy(dtype=np.float64)
    keep = n_included_slides_detected(Xraw) >= MIN_SLIDES
    Xraw = Xraw[keep]
    Xraw, _ = group_median_impute(Xraw, groups)

    def save_outputs(cluster_ids, tag):
        n_cl = len(np.unique(cluster_ids))
        out = pd.DataFrame({
            "sample": samples, "group": groups,
            "cluster": cluster_ids,
            "UMAP_1": X_umap[:, 0], "UMAP_2": X_umap[:, 1],
        })

        ann_rows = []
        for c in np.unique(cluster_ids):
            m = cluster_ids == c
            n_c = int(m.sum())
            row = {"cluster": c, "n": n_c, "frac_of_total": n_c / len(cluster_ids)}
            for g in ["WT", "mdx", "AAV9", "LICA1"]:
                row[f"frac_{g}"] = (groups[m] == g).mean()
            top_samples = pd.Series(samples[m]).value_counts()
            row["top_samples"] = "; ".join(f"{s}({k})" for s, k in top_samples.head(3).items())
            for name, col in KEY_FEATURES.items():
                row[f"mean_{name}"] = Xraw[m, feat_cols.index(col)].mean()
            ann_rows.append(row)
        ann = pd.DataFrame(ann_rows).sort_values("n", ascending=False)
        out_ann = DIR_VISUALIZATIONS / f"{MUSCLE}_cluster_annotation_{tag}.csv"
        ann.to_csv(out_ann, index=False)

        def scatter_umap(ax, color_arr, title="", categorical=None, cmap_name=None):
            if categorical is not None:
                for val in categorical:
                    m = color_arr == val
                    ax.scatter(X_umap[m, 0], X_umap[m, 1], s=3, alpha=0.25,
                               color=GROUP_COLORS.get(val, "gray"), label=f"{val} (n={int(m.sum())})",
                               edgecolors="none")
                ax.legend(title=None, fontsize=8, loc="best")
            else:
                sc = ax.scatter(X_umap[:, 0], X_umap[:, 1], c=color_arr, cmap=get_cmap(cmap_name),
                                s=1, alpha=0.3, edgecolors="none")
                plt.colorbar(sc, ax=ax, fraction=0.046)
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            ax.set_title(title, fontsize=10)
            ax.tick_params(labelsize=8)

        order = list(ann["cluster"])
        code = {c: i for i, c in enumerate(order)}
        disp_plot = np.array([code[c] for c in cluster_ids])
        cmap = get_cmap("tab20")

        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        scatter_umap(axes[0, 0], disp_plot, cmap_name="tab20",
                     title=f"clusters (n={n_cl}, {args.method}, norm={args.norm})")
        for t, c in enumerate(order):
            axes[0, 0].text(0.02, 0.98 - 0.045 * t, f"{c} (n={int(ann.iloc[t]['n'])})",
                            transform=axes[0, 0].transAxes, fontsize=7,
                            color=cmap(t / max(n_cl - 1, 1)), va="top")
        scatter_umap(axes[0, 1], groups, categorical=["WT", "mdx", "AAV9", "LICA1"], title="Group")
        sample_codes = {s: i for i, s in enumerate(QUA_SAMPLES)}
        scatter_umap(axes[1, 0], np.array([sample_codes[s] for s in samples]),
                     cmap_name="tab20", title="Sample (batch check)")
        scatter_umap(axes[1, 1], disp_plot, cmap_name="tab20", title="clusters (reference)")
        fig.suptitle(f"{MUSCLE} — UMAP of fiber features (norm={args.norm}, {args.method}; "
                      f"{len(X_pca):,} fibers, >= {MIN_SLIDES}/{len(included_slide_blocks())} included slides)")
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out_png = DIR_VISUALIZATIONS / f"{MUSCLE}_umap_clusters_{tag}.png"
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        print(f"Saved {out_png.name} and {out_ann.name}")

        out_csv = DIR_VISUALIZATIONS / f"{MUSCLE}_umap_clusters_{tag}.csv"
        out.to_csv(out_csv, index=False)

        print(f"\n[{tag}] {n_cl} clusters:")
        print(ann[["cluster", "n", "frac_of_total", "frac_WT", "frac_mdx", "frac_AAV9",
                   "frac_LICA1", "mean_area", "mean_Laminin", "mean_Myh7",
                   "mean_Myh2", "mean_Myh4", "mean_IgG"]].to_string(index=False,
                                                                     float_format=lambda v: f"{v:8.3f}"))

    if args.method == "gmm":
        print(f"\nGMM clustering (k={args.kmin}..{args.kmax}, one figure per k) ...")
        for k in range(args.kmin, args.kmax + 1):
            gm = GaussianMixture(n_components=k, covariance_type="full",
                                 n_init=1, max_iter=200, random_state=RANDOM_STATE)
            gm.fit(X_pca)
            print(f"  GMM k={k}: BIC {gm.bic(X_pca):.0f}")
            save_outputs(gm.predict(X_pca), f"{args.norm}_gmm_k{k}")
    else:
        resolutions = [float(r) for r in args.resolutions.split(",")]
        nn = NearestNeighbors(n_neighbors=N_NEIGHBORS, metric="euclidean", n_jobs=8)
        nn.fit(X_pca)
        _, idx = nn.kneighbors(X_pca)
        src = [(i, int(j)) for i in range(len(X_pca)) for j in idx[i] if j > i]
        graph = ig.Graph(n=len(X_pca))
        graph.add_edges(src)
        print(f"  graph: {graph.vcount()} nodes, {graph.ecount()} edges")
        results = {}
        for res in resolutions:
            part = leidenalg.find_partition(graph, leidenalg.RBConfigurationVertexPartition,
                                            weights=None, resolution_parameter=res,
                                            seed=RANDOM_STATE)
            n = len(np.unique(part.membership))
            results[res] = (part, n)
            print(f"  resolution {res:5.2f} -> {n} clusters")
        in_range = {r: n for r, (_, n) in results.items() if args.kmin <= n <= args.kmax}
        if in_range:
            best_res = min(in_range, key=lambda r: abs(in_range[r] - (args.kmin + args.kmax) / 2))
        else:
            best_res = min(results, key=lambda r: abs(results[r][1] - (args.kmin + args.kmax) / 2))
        part, n_cl = results[best_res]
        save_outputs(np.array(part.membership), f"{args.norm}_leiden_res{best_res}")

    print("\nDone.")


if __name__ == "__main__":
    main()
