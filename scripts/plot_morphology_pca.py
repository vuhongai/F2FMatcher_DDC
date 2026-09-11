import os, sys, pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.append("/DATA/F2FMatcher_DDC")
from config.ddc_config import *

MUSCLE = "QUA"
SLIDE_REF = 1
SLIDE_HE = 6

DIR_VISUALIZATIONS = Path("/DATA/F2FMatcher_DDC/visualizations")
DIR_COMBINED = Path(f"/DATA/F2FMatcher_DDC/results/{MUSCLE}/features_combined")
DIR_PAIR_OUTPUT = Path(PAIR_DIRS_BASE_QUA)

COMPARTMENT_ORDER = ["whole", "mem", "cyto1", "cyto2"]
GROUP_COLORS = {"WT": "tab:blue", "mdx": "tab:red", "AAV9": "tab:green", "LICA1": "tab:orange"}


def get_filename(name_sample, list_names):
    names = [n for n in list_names if name_sample.upper() in n.upper()]
    if len(names) == 0:
        return None
    return names[0].split(".")[0]


def he_channel_offset():
    """Channel index of the first HE_10x channel in the combined vector."""
    offset = 0
    for slide in SLIDES:
        if slide == SLIDE_HE:
            return offset
        offset += len(SLIDES[slide]["stainings"])


def load_sample(sample):
    """Fibers detected in both slide 1 and slide 6: mask features + HE features."""
    dir_czi_source = CZI_BASE_DIR_QUA
    img1 = get_filename(sample, [f.split(".czi")[0] for f in os.listdir(dir_czi_source / SLIDES[SLIDE_REF]["czi_dir"]) if f.endswith(".czi")])
    img6 = get_filename(sample, [f.split(".czi")[0] for f in os.listdir(dir_czi_source / SLIDES[SLIDE_HE]["czi_dir"]) if f.endswith(".czi")])
    if img1 is None or img6 is None:
        print(f"[skip] {sample}: image not found")
        return None

    pair_dir = DIR_PAIR_OUTPUT / f"{img1}___vs___{img6}"
    path_pairs = pair_dir / "paired_labels.pkl"
    if not path_pairs.exists():
        print(f"[skip] {sample}: no slide 1->6 mapping ({pair_dir.name})")
        return None
    with open(path_pairs, "rb") as f:
        matched_ref_labels = {p[0] for p in pickle.load(f)}

    with open(DIR_COMBINED / f"{sample}.pkl", "rb") as f:
        combined = pickle.load(f)

    n_ch_he = len(SLIDES[SLIDE_HE]["stainings"])
    he_start = 15 + he_channel_offset() * 36
    he_end = he_start + n_ch_he * 36

    rows = []
    for label_id, vec in combined.items():
        if label_id not in matched_ref_labels:
            continue
        rows.append([label_id] + list(vec[:15]) + list(vec[he_start:he_end]))

    df = pd.DataFrame(rows)
    he_cols = []
    for ch, stain in SLIDES[SLIDE_HE]["stainings"].items():
        for comp in COMPARTMENT_ORDER:
            for stat in FEATURE_STATISTICS:
                he_cols.append(f"HE_{stain}_{comp}_{stat}")
    df.columns = ["label_id"] + MASK_FEATURES + he_cols
    df["sample"] = sample
    df["group"] = get_group(sample)
    return df


def main():
    DIR_VISUALIZATIONS.mkdir(parents=True, exist_ok=True)

    frames = []
    for sample in QUA_SAMPLES:
        df = load_sample(sample)
        if df is not None:
            frames.append(df)
            print(f"{sample}: {len(df)} fibers (detected in slide 1+6)")
    data = pd.concat(frames, ignore_index=True)
    print(f"\nTotal: {len(data)} fibers from {data['sample'].nunique()} samples")

    feat_cols = MASK_FEATURES + [c for c in data.columns if c.startswith("HE_")]
    X = data[feat_cols].to_numpy(dtype=np.float64)
    n_nan = np.isnan(X).sum()
    if n_nan:
        print(f"Imputing {n_nan} NaN value(s) ({n_nan / X.size:.2%}) with column medians")
        col_med = np.nanmedian(X, axis=0)
        col_med[np.isnan(col_med)] = 0.0
        for j in range(X.shape[1]):
            nan_rows = np.isnan(X[:, j])
            X[nan_rows, j] = col_med[j]

    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=ANALYSIS_CONFIG["random_state"])
    scores = pca.fit_transform(Xs)
    print(f"Explained variance: PC1 {pca.explained_variance_ratio_[0]:.1%}, PC2 {pca.explained_variance_ratio_[1]:.1%}")

    out_csv = DIR_VISUALIZATIONS / f"{MUSCLE}_morph_HE_pca_scores.csv"
    pd.DataFrame({
        "sample": data["sample"],
        "group": data["group"],
        "PC1": scores[:, 0],
        "PC2": scores[:, 1],
    }).to_csv(out_csv, index=False)
    print(f"Saved scores to {out_csv}")

    group_order = [g for g in ["WT", "mdx", "AAV9", "LICA1"] if g in set(data["group"])]
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_axes([0.08, 0.12, 0.78, 0.78])
    for g in group_order:
        dg = data[data["group"] == g]
        ax.scatter(scores[dg.index, 0], scores[dg.index, 1], s=8, alpha=0.3,
                   color=GROUP_COLORS.get(g, "gray"), label=f"{g} (n={len(dg)})", edgecolors="none")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(title="Group")
    fig.suptitle(f"{MUSCLE} — PCA on morphology + HE_10x (fibers detected in slide 1+6)")

    ax_top = fig.add_axes([0.08, 0.92, 0.78, 0.06])
    ax_top.hist(scores[:, 0], bins=60, density=True, histtype="step", color="black")
    ax_top.set_xticks([])
    ax_top.set_title("PC1", fontsize=9)
    ax_right = fig.add_axes([0.88, 0.12, 0.06, 0.78])
    ax_right.hist(scores[:, 1], bins=60, density=True, histtype="step", color="black")
    ax_right.set_yticks([])
    ax_right.set_title("PC2", fontsize=9)

    out_png = DIR_VISUALIZATIONS / f"{MUSCLE}_morph_HE_pca.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot to {out_png}")

    # 4 subplots: one per group, shared axes
    fig2, axes = plt.subplots(2, 2, figsize=(12, 11), sharex=True, sharey=True)
    for axi, g in zip(axes.flat, group_order):
        dg = data[data["group"] == g]
        axi.scatter(scores[dg.index, 0], scores[dg.index, 1], s=8, alpha=0.3,
                    color=GROUP_COLORS.get(g, "gray"), edgecolors="none")
        axi.set_title(f"{g} (n={len(dg)})", fontsize=12)
        axi.set_xlabel("PC1")
        axi.set_ylabel("PC2")
    fig2.suptitle(f"{MUSCLE} — PCA on morphology + HE_10x, per group (fibers detected in slide 1+6)")
    fig2.tight_layout(rect=[0, 0, 1, 0.96])
    out_png2 = DIR_VISUALIZATIONS / f"{MUSCLE}_morph_HE_pca_groups.png"
    fig2.savefig(out_png2, dpi=200)
    plt.close(fig2)
    print(f"Saved per-group plot to {out_png2}")

    print("\nPer-group PC statistics:")
    for pc in (0, 1):
        groups_vals = [scores[data["group"] == g, pc] for g in group_order]
        h, p = kruskal(*groups_vals)
        print(f"PC{pc+1}: Kruskal-Wallis H={h:.1f}, p={p:.2e}")
        for g in group_order:
            v = scores[data["group"] == g, pc]
            print(f"  {g:6s} mean {v.mean():8.3f}  std {v.std():7.3f}  n {len(v)}")


if __name__ == "__main__":
    main()
