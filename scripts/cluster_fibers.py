"""Build/load the UMAP cache: features -> filter -> impute -> normalize -> PCA -> UMAP.
No plotting here — see plot_umap.py (reads the precomputed cache).
"""
import sys, argparse
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from umap import UMAP

sys.path.insert(0, "/DATA/F2FMatcher_DDC")
sys.path.insert(0, "/DATA/F2FMatcher_DDC/scripts")
from config.ddc_config import *
from plot_pca_all_features import (
    load_all, feature_columns, slide_channel_blocks,
    group_median_impute,
    COMPARTMENT_ORDER, FEATURE_STATISTICS,
    MUSCLE, DIR_VISUALIZATIONS,
)

N_NEIGHBORS = 15
MIN_SLIDES = 4
RANDOM_STATE = ANALYSIS_CONFIG["random_state"]
CACHE_DIR = Path(f"/DATA/F2FMatcher_DDC/results/{MUSCLE}/clustering_cache")

# Channel presets for UMAP construction (slide -> stainings, None = all).
CHANNEL_PRESETS = {
    # full spec: everything except S1 Dystrophin and slide 8
    "ch18": {
        1: ["DAPI", "Laminin", "Collagen4"],
        2: ["IgG", "CD11b"],
        3: None,
        4: None,
        6: None,
        7: None,
    },
    # curated: nuclei, immune, fiber type, HE brightfield
    "ch9": {
        1: ["DAPI"],
        2: ["IgG", "CD11b"],
        4: ["Myh7", "Myh2", "Myh4"],
        6: None,
    },
    # unbiased: morphology (MASK_FEATURES) + HE_10x brightfield only
    "ch3": {
        6: None,
    },
}
INCLUDED_SLIDES = dict(CHANNEL_PRESETS["ch18"])


def set_preset(name):
    """Switch the channel preset; returns the new feat tag (callers must use the
    return value — a module-level `from ... import FEAT_TAG` would go stale)."""
    global INCLUDED_SLIDES, FEAT_TAG
    INCLUDED_SLIDES = dict(CHANNEL_PRESETS[name])
    FEAT_TAG = f"ch{len(included_channels())}"
    return FEAT_TAG


def min_slides_effective():
    """Detection filter threshold: MIN_SLIDES, capped at the number of included slides."""
    return min(MIN_SLIDES, len(included_slide_blocks()))


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
    """Per fiber: number of included slides with at least one observed value (subset layout)."""
    n = np.zeros(X.shape[0], dtype=int)
    for lo, hi in included_slide_blocks().values():
        n += (~np.isnan(X[:, lo:hi])).any(axis=1)
    return n


def n_included_slides_detected_full(X):
    """Same filter, but X in the full 807-column layout (15 mask + 22 channels x 36)."""
    blocks = slide_channel_blocks()
    n = np.zeros(X.shape[0], dtype=int)
    for slide in INCLUDED_SLIDES:
        st, nc = blocks[slide]
        block = X[:, 15 + st * 36:15 + (st + nc) * 36]
        n += (~np.isnan(block)).any(axis=1)
    return n


FEAT_TAG = f"ch{len(included_channels())}"

# Markers for feature-coloring plots. str = single channel value;
# list = brightfield RGB, value = mean of the 3 channel whole-ROI means.
MARKERS = {
    "area": "area",
    "DAPI": "S1_DAPI_cyto2_mean",
    "Laminin": "S1_Laminin_mem_mean",
    "Dystrophin": "S1_Dystrophin_mem_mean",
    "Collagen4": "S1_Collagen4_mem_mean",
    "IgG": "S2_IgG_whole_mean",
    "CD11b": "S2_CD11b_mem_mean",
    "NADH": ["S3_R_whole_mean", "S3_G_whole_mean", "S3_B_whole_mean"],
    "WGA": "S4_WGA_mem_mean",
    "Myh7": "S4_Myh7_whole_mean",
    "Myh2": "S4_Myh2_whole_mean",
    "Myh4": "S4_Myh4_whole_mean",
    "HE_10x": ["S6_R_whole_mean", "S6_G_whole_mean", "S6_B_whole_mean"],
    "COX": ["S7_R_whole_mean", "S7_G_whole_mean", "S7_B_whole_mean"],
    "LAMP2": "S8_LAMP2_whole_mean",
    "LGALS3": "S8_LGALS3_whole_mean",
    "SQSTM1": "S8_SQSTM1_whole_mean",
}


def marker_values(Xfull, cols):
    """Per-marker value for each fiber (RGB markers = mean of R/G/B whole means)."""
    idx = {c: i for i, c in enumerate(cols)}
    out = {}
    for name, spec in MARKERS.items():
        if isinstance(spec, str):
            out[name] = Xfull[:, idx[spec]]
        else:
            out[name] = np.mean([Xfull[:, idx[c]] for c in spec], axis=0)
    return out


def resolve_feature(name, Xfull, cols):
    """Value array for a raw column name or a MARKERS name."""
    idx = {c: i for i, c in enumerate(cols)}
    if name in idx:
        return Xfull[:, idx[name]]
    if name in MARKERS:
        spec = MARKERS[name]
        if isinstance(spec, str):
            return Xfull[:, idx[spec]]
        return np.mean([Xfull[:, idx[c]] for c in spec], axis=0)
    raise SystemExit(f"unknown feature/marker: {name}")


def feature_observed_mask(name, Xfull, cols):
    """Per-fiber: True if the feature's channel(s) have any observed (non-NaN) value.
    Xfull must be the PRE-imputation matrix in the full 807-column layout."""
    blocks = slide_channel_blocks()

    def chan_block(col):
        slide, stain = col.split("_")[0][1:], col.split("_")[1]  # strip "S" prefix
        st0 = blocks[int(slide)][0]
        ci = list(SLIDES[int(slide)]["stainings"].values()).index(stain)
        base = 15 + (st0 + ci) * 36
        return base, base + 36

    if name in MASK_FEATURES:
        return np.ones(len(Xfull), dtype=bool)
    if name in cols:
        cols_list = [name]
    elif name in MARKERS:
        spec = MARKERS[name]
        cols_list = spec if isinstance(spec, list) else [spec]
    else:
        raise SystemExit(f"unknown feature/marker: {name}")
    obs = np.zeros(len(Xfull), dtype=bool)
    for c in cols_list:
        lo, hi = chan_block(c)
        obs |= (~np.isnan(Xfull[:, lo:hi])).any(axis=1)
    return obs


def normalize(X, mode, feat_cols, groups=None):
    """mode: none | zscore | groupzscore (per-feature, per-group) | pooled
    (center per group, whiten with the pooled covariance) | wtshift
    (shift WT onto the treated-group mean, keep treated mean structure)"""
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
    if mode == "pooled":
        Xc = X.copy()
        for g in np.unique(groups):
            Xc[groups == g] -= np.nanmean(X[groups == g], axis=0)
        cov = np.cov(Xc, rowvar=False)
        w, V = np.linalg.eigh(cov)
        w = np.clip(w, 1e-8, None)
        return Xc @ (V / np.sqrt(w))
    if mode == "wtshift":
        # shift WT onto the mean of the 3 treated groups (removes the WT-vs-rest
        # axis) while KEEPING the mdx/AAV9/LICA1 mean structure (treatment effect)
        Xc = X.copy()
        rest = groups != "WT"
        shift = np.nanmean(X[groups == "WT"], axis=0) - np.nanmean(X[rest], axis=0)
        Xc[groups == "WT"] -= shift
        mu = np.nanmean(Xc, axis=0)
        sd = np.nanstd(Xc, axis=0)
        sd[sd == 0] = 1.0
        return (Xc - mu) / sd
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def build_cache(norm, n_pca, mindist, feat_tag=FEAT_TAG, n_neighbors=N_NEIGHBORS):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / f"{feat_tag}_{norm}_pca{n_pca}_mindist{mindist}_nn{n_neighbors}.npz"
    df = load_all()
    feat_cols = included_feature_columns()
    if feat_tag.endswith("_nomorph"):
        feat_cols = [c for c in feat_cols if c not in MASK_FEATURES]
    X = df[feat_cols].to_numpy(dtype=np.float64)
    groups = df["group"].to_numpy()
    samples = df["sample"].to_numpy()
    print(f"Total fibers: {len(df)} ({len(feat_cols)} features, {feat_tag})")

    n_incl = len(included_slide_blocks())
    keep = n_included_slides_detected(X) >= min_slides_effective()
    print(f"Keeping {keep.sum()} fibers detected in >= {min_slides_effective()}/{n_incl} included slides")
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

    print(f"\nUMAP (n_neighbors={n_neighbors}, min_dist={mindist}) ...")
    reducer = UMAP(n_neighbors=n_neighbors, min_dist=mindist, metric="euclidean",
                   random_state=RANDOM_STATE)
    X_umap = reducer.fit_transform(X_pca)

    np.savez_compressed(cache_path, X_pca=X_pca, X_umap=X_umap,
                        samples=samples, groups=groups, feat_cols=np.array(feat_cols))
    print(f"Saved cache to {cache_path}")
    return cache_path


def load_cache(norm, n_pca, mindist, feat_tag=FEAT_TAG, n_neighbors=N_NEIGHBORS):
    cache_path = CACHE_DIR / f"{feat_tag}_{norm}_pca{n_pca}_mindist{mindist}_nn{n_neighbors}.npz"
    if not cache_path.exists():
        build_cache(norm, n_pca, mindist, feat_tag, n_neighbors)
    print(f"Loading cache {cache_path}")
    z = np.load(cache_path, allow_pickle=True)
    return {
        "X_pca": z["X_pca"], "X_umap": z["X_umap"],
        "samples": z["samples"], "groups": z["groups"],
        "feat_cols": [str(c) for c in z["feat_cols"]],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--norm", default="zscore",
                        choices=["none", "zscore", "groupzscore", "pooled", "wtshift"])
    parser.add_argument("--preset", default="ch18", choices=list(CHANNEL_PRESETS))
    parser.add_argument("--n-pca", type=int, default=30)
    parser.add_argument("--mindist", type=float, default=0.0)
    parser.add_argument("--n-neighbors", type=int, default=N_NEIGHBORS)
    parser.add_argument("--no-morph", action="store_true",
                        help="exclude the 15 morphology features from UMAP construction")
    parser.add_argument("--rebuild-cache", action="store_true")
    args = parser.parse_args()

    set_preset(args.preset)
    feat_tag = FEAT_TAG + ("_nomorph" if args.no_morph else "")
    if args.rebuild_cache:
        build_cache(args.norm, args.n_pca, args.mindist, feat_tag, args.n_neighbors)
    else:
        load_cache(args.norm, args.n_pca, args.mindist, feat_tag, args.n_neighbors)
    print("Done.")


if __name__ == "__main__":
    main()
