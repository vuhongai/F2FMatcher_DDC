"""Build the analysis bundle used by the Q1/Q2 notebooks.

Everything is derived from PRECOMPUTED artifacts (no re-clustering, no re-UMAP):
  1. the clustering cache  results/QUA/clustering_cache/ch3_zscore_pca30_mindist0.0_nn15.npz
     (unbiased setup: 15 morphology + HE_10x RGB = 123 features, zscore, PCA30, UMAP nn15/mindist0)
  2. the combined per-fiber features  results/QUA/features_combined/*.pkl (807 cols)

The GMM (k=6) is re-fit on the cached PCA scores with the exact parameters used by
plot_umap.py (full covariance, n_init=1, max_iter=200, random_state=42) -> deterministic.

Output: results/QUA/analysis_bundle.npz
  Xpre       float32  (n, 807)  raw features, pre-imputation (NaN = fiber not mapped on that channel)
  X_umap     float32  (n, 2)
  groups     (n,)     WT / mdx / AAV9 / LICA1
  samples    (n,)     QUAGxx
  cluster_ids (n,)    GMM k=6 labels (0..5)
  cols       (807,)   feature column names

Usage:  python scripts/build_analysis_bundle.py
"""
import sys, time
from pathlib import Path

import numpy as np
from sklearn.mixture import GaussianMixture

REPO = Path("/DATA/F2FMatcher_DDC")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from plot_pca_all_features import load_all, feature_columns  # noqa: E402
from cluster_fibers import (  # noqa: E402
    CACHE_DIR, set_preset, n_included_slides_detected_full, min_slides_effective, RANDOM_STATE,
)

MUSCLE = "QUA"
K_GMM = 6


def main():
    t0 = time.time()
    set_preset("ch3")
    cache_path = CACHE_DIR / "ch3_zscore_pca30_mindist0.0_nn15.npz"
    assert cache_path.exists(), f"missing precomputed cache: {cache_path}"
    z = np.load(cache_path, allow_pickle=True)
    X_pca, X_umap = z["X_pca"], z["X_umap"]
    groups, samples = z["groups"], z["samples"]
    print(f"[1/4] cache loaded: {len(X_umap):,} fibers  ({time.time()-t0:.0f}s)")

    t1 = time.time()
    gm = GaussianMixture(n_components=K_GMM, covariance_type="full",
                         n_init=1, max_iter=200, random_state=RANDOM_STATE)
    gm.fit(X_pca)
    cluster_ids = gm.predict(X_pca)
    sizes = np.bincount(cluster_ids)
    print(f"[2/4] GMM k={K_GMM} on cached PCA30  ({time.time()-t1:.0f}s)")
    print("      cluster sizes:", dict(enumerate(sizes)))

    t2 = time.time()
    df = load_all()
    cols = feature_columns()
    X = df[cols].to_numpy(dtype=np.float64)
    g_all = df["group"].to_numpy()
    keep = n_included_slides_detected_full(X) >= min_slides_effective()
    X, g_all, s_all = X[keep], g_all[keep], df["sample"].to_numpy()[keep]
    assert len(X) == len(X_umap), f"row mismatch: {len(X)} vs {len(X_umap)}"
    assert np.array_equal(g_all, groups) and np.array_equal(s_all, samples), \
        "row order mismatch with cache"
    print(f"[3/4] combined features aligned to cache: {len(X):,} fibers  ({time.time()-t2:.0f}s)")

    out = REPO / "results" / MUSCLE / "analysis_bundle.npz"
    np.savez_compressed(out,
                        Xpre=X.astype(np.float32),
                        X_umap=X_umap.astype(np.float32),
                        groups=groups, samples=samples,
                        cluster_ids=cluster_ids.astype(np.int32),
                        cols=np.array(cols))
    print(f"[4/4] saved {out}  ({time.time()-t0:.0f}s total)")
    for c in range(K_GMM):
        m = cluster_ids == c
        fr = "  ".join(f"{g} {(groups[m] == g).mean():5.1%}" for g in ["WT", "mdx", "AAV9", "LICA1"])
        print(f"      C{c+1}: n={int(m.sum()):>7,}  {fr}")


if __name__ == "__main__":
    main()
