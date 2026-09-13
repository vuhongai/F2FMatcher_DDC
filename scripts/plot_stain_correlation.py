"""Density scatter of two stainings vs each other, one panel per group.

Usage:
    python scripts/plot_stain_correlation.py --stain1 COX --stain2 NADH

Stain names accept a MARKERS name (COX, NADH, Myh7, Dystrophin, Laminin, ...) or a
raw feature column (e.g. S7_R_whole_mean). RGB markers (COX/NADH/HE_10x) use the mean
of the R/G/B whole-ROI means. Only fibers where BOTH stainings are observed (mapped on
both slides) are plotted. Self-contained (density scatter inlined, no libLICA import).
"""
import sys, argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interpn
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, "/DATA/F2FMatcher_DDC")
sys.path.insert(0, "/DATA/F2FMatcher_DDC/scripts")
from config.ddc_config import *
from plot_pca_all_features import load_all, feature_columns
from cluster_fibers import (
    set_preset, resolve_feature, feature_observed_mask,
    n_included_slides_detected_full, min_slides_effective, DIR_VISUALIZATIONS, MUSCLE,
)

GROUPS = ["WT", "mdx", "AAV9", "LICA1"]
GROUP_COLORS = {"WT": "#0000C0", "mdx": "#FF6000", "AAV9": "#C0C000", "LICA1": "#008000"}
DIR_VISUALIZATIONS = DIR_VISUALIZATIONS / "QUA_corr"

def density_scatter(ax, x, y, bins=100, subsample=200000, s=1.0, **kwargs):
    """Scatter colored by local 2D-histogram density (spline-interpolated to points).
    `s` = marker size (points^2); raise it to make points bigger."""
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 10:
        return
    xf, yf = x[m], y[m]
    if n > subsample:
        rng = np.random.default_rng(42)
        idx = rng.choice(n, subsample, replace=False)
        xf, yf = xf[idx], yf[idx]
    data, xe, ye = np.histogram2d(xf, yf, bins=bins, density=True)
    z = interpn(
        (0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])),
        data, np.vstack([xf, yf]).T,
        method="splinef2d", bounds_error=False,
    )
    z[np.isnan(z)] = 0.0
    order = z.argsort()
    ax.scatter(xf[order], yf[order], s=s, c=z[order], **kwargs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stain1", required=True, help="x-axis staining (MARKERS name or raw column)")
    ap.add_argument("--stain2", required=True, help="y-axis staining (MARKERS name or raw column)")
    ap.add_argument("--preset", default="ch3")
    ap.add_argument("--bins", type=int, default=50)
    ap.add_argument("--subsample", type=int, default=200000)
    ap.add_argument("--size", type=float, default=10.0, help="marker size (points^2)")
    args = ap.parse_args()

    set_preset(args.preset)
    df = load_all()
    cols = feature_columns()
    Xfull = df[cols].to_numpy(dtype=np.float64)
    groups = df["group"].to_numpy()

    # same fiber set as the clustering cache (ch3 = detected on the included slides)
    keep = n_included_slides_detected_full(Xfull) >= min_slides_effective()
    Xfull, groups = Xfull[keep], groups[keep]
    print(f"fibers: {len(Xfull):,}")

    v1 = resolve_feature(args.stain1, Xfull, cols)
    v2 = resolve_feature(args.stain2, Xfull, cols)
    o1 = feature_observed_mask(args.stain1, Xfull, cols)
    o2 = feature_observed_mask(args.stain2, Xfull, cols)
    both = o1 & o2
    print(f"both stainings observed: {int(both.sum()):,} "
          f"({args.stain1} {int(o1.sum()):,} / {args.stain2} {int(o2.sum()):,})")

    # shared axis limits from the union of observed data (0.5-99.5 pct)
    lx = np.nanpercentile(v1[both], [0.5, 99.5])
    ly = np.nanpercentile(v2[both], [0.5, 99.5])

    DIR_VISUALIZATIONS.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(14,4), sharex=True, sharey=True)
    sm = None
    print(f"\n{'group':6s} {'n':>7s} {'pearson r':>10s} {'spearman rho':>13s}")
    for ax, g in zip(axes, GROUPS):
        m = (groups == g) & both
        x, y = v1[m], v2[m]
        density_scatter(ax, x, y, s=args.size, bins=args.bins, subsample=args.subsample,
                        cmap="viridis", edgecolors="none")
        if len(x) >= 3:
            r, _ = pearsonr(x, y)
            rho, _ = spearmanr(x, y)
        else:
            r = rho = np.nan
        print(f"{g:6s} {len(x):>7,} {r:>10.3f} {rho:>13.3f}")
        ax.set_xlim(lx); ax.set_ylim(ly)
        ax.set_title(g, fontsize=13, fontweight="bold", color="black")
        ax.text(0.04, 0.96, f"n={len(x):,}\nr={r:.2f}\nrho={rho:.2f}",
                transform=ax.transAxes, fontsize=12, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.4))
        ax.set_xlabel(args.stain1, fontsize=12)
        ax.tick_params(labelsize=11)
    axes[0].set_ylabel(args.stain2, fontsize=12)

    fig.suptitle(f"{args.stain1} vs {args.stain2}", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = DIR_VISUALIZATIONS / f"{MUSCLE}_corr_{args.stain1}_vs_{args.stain2}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
