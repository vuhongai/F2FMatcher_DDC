import os, sys, pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.append("/DATA/F2FMatcher_DDC")
from config.ddc_config import *

MUSCLE = "QUA"
DIR_FEATURES = Path(f"/DATA/F2FMatcher_DDC/results/{MUSCLE}/features")
DIR_VISUALIZATIONS = Path("/DATA/F2FMatcher_DDC/visualizations")
GROUP_COLORS = {"WT": "tab:blue", "mdx": "tab:red", "AAV9": "tab:green", "LICA1": "tab:orange"}


def get_filename(name_sample, list_names):
    names = [n for n in list_names if name_sample.upper() in n.upper()]
    if len(names) == 0:
        return None
    return names[0].split(".")[0]


def main():
    DIR_VISUALIZATIONS.mkdir(parents=True, exist_ok=True)
    dir_czi_source = CZI_BASE_DIR_QUA

    rows = []
    for slide in SLIDES:
        list_images = [f.split(".czi")[0] for f in os.listdir(dir_czi_source / SLIDES[slide]["czi_dir"]) if f.endswith(".czi")]
        for sample in QUA_SAMPLES:
            img = get_filename(sample, list_images)
            if img is None:
                continue
            for ch, stain in SLIDES[slide]["stainings"].items():
                path = DIR_FEATURES / f"slide_{slide}" / f"{img}_c{ch}.pkl"
                if not path.exists():
                    continue
                with open(path, "rb") as f:
                    feats = pickle.load(f)
                vals = np.array([ft["intensity"]["whole"]["mean"] for ft in feats], dtype=np.float64)
                vals = vals[~np.isnan(vals)]
                rows.append({
                    "channel": f"S{slide}_{stain}",
                    "sample": sample,
                    "group": get_group(sample),
                    "n_fibers": len(vals),
                    "median": np.median(vals),
                    "p10": np.percentile(vals, 10),
                    "p90": np.percentile(vals, 90),
                })
    df = pd.DataFrame(rows)

    # per channel: how much do sample medians differ?
    summary = df.groupby("channel").agg(
        med_min=("median", "min"),
        med_max=("median", "max"),
        med_mean=("median", "mean"),
        med_std=("median", "std"),
    )
    summary["max_min_ratio"] = summary["med_max"] / summary["med_min"]
    summary["cv"] = summary["med_std"] / summary["med_mean"]
    summary = summary.sort_values("max_min_ratio", ascending=False)
    print("Per-channel spread of per-sample medians (whole-fiber mean intensity):")
    print(summary.to_string(float_format=lambda v: f"{v:10.3f}"))

    # heatmap: rows = channels, cols = samples, value = log10(median)
    samples = list(QUA_SAMPLES)
    channels = list(SLIDES.keys()) and [f"S{s}_{st}" for s in SLIDES for st in SLIDES[s]["stainings"].values()]
    piv = df.pivot(index="channel", columns="sample", values="median").reindex(index=channels, columns=samples)
    logpiv = np.log10(piv.replace(0, np.nan))

    fig, ax = plt.subplots(figsize=(13, 7))
    im = ax.imshow(logpiv.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=8)
    ax.set_yticks(range(len(channels)))
    ax.set_yticklabels(channels, fontsize=8)
    for i, s in enumerate(samples):
        g = get_group(s)
        ax.plot([-0.5, len(channels) - 0.5 + 0.5], [i + 0.35, i + 0.35], color="none")
    # group color bar on top
    for i, s in enumerate(samples):
        ax.add_patch(plt.Rectangle((i - 0.5, len(channels) - 0.32), 1, 0.25,
                                   color=GROUP_COLORS[get_group(s)], clip_on=False))
    ax.set_ylim(len(channels) - 0.5, -0.5)
    ax.set_title(f"{MUSCLE} — per-sample median whole-fiber mean intensity per channel (log10)\n"
                 f"rows = channels, cols = samples (top bar = group)")
    fig.colorbar(im, ax=ax, label="log10 median intensity")
    fig.tight_layout()
    out = DIR_VISUALIZATIONS / f"{MUSCLE}_slide_variability_heatmap.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"\nSaved heatmap to {out}")

    # per-slide (across channels of the slide) consistency
    print("\nWorst channels (largest between-sample median spread):")
    for ch, r in summary.head(8).iterrows():
        print(f"  {ch:35s} max/min {r['max_min_ratio']:6.2f}  cv {r['cv']:5.2f}")


if __name__ == "__main__":
    main()
