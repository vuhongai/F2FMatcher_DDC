"""Training-illustration figures for §A.3 (VAE) and §A.4 (classifier).

Sources (old F2FMatcher project on CIFS + local production checkpoints):
  - VAE reconstructions: training_joinedVAE2/logs/plot_model_epoch*_{original,augmented}.png
  - VAE loss log:        training_joinedVAE2/logs/training_log.txt
  - classifier log:      training_fibermatcher_cls_2/3_train_fibermatcher.ipynb outputs
  - test pairs:          training_fibermatcher_cls_2/dataset_pair_test (pickled)
  - embeddings:          training_fibermatcher_cls_2/embed/*.npy (256-d)
  - classifier ckpt:     /DATA/F2FMatcher/models/fibermatcher_cls_2.pth

Outputs (manuscript/figures/):
  figA3a_vae_reconstruction.{png,pdf}       reconstruction at epochs 1/10/30, original + augmented
  figA3b_vae_training_curves.{png,pdf}      train/val loss: total + recon/KL/latent
  figA3c_vae_latent_umap.{png,pdf}          UMAP of 54 held-out ROIs coloured by fibre identity,
                                            latent-distance evidence, example same-fibre crops
  figA4a_classifier_training_curves.{png,pdf}  train/val loss + val F1 (best epoch marked)
  figA4b_classifier_score_distribution.{png,pdf} score density + ROC on held-out test pairs
"""
import argparse
import json
import os
import pickle
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "manuscript" / "figures"
OLD = Path("/media/DATABRUT/DB_DDC/serverGPU/Cache_GPU_Ai/fiber_matcher/F2FMatcher")
VAE_LOGS = OLD / "training_joinedVAE2" / "logs"
CLS_DIR = OLD / "training_fibermatcher_cls_2"
CKPT = Path("/DATA/F2FMatcher/models/fibermatcher_cls_2.pth")

EPOCHS_SHOWN = [1, 10, 30]
COLTITLES = ["flow_x", "flow_y", "mask", "flow_x (rec)", "flow_y (rec)", "mask (rec)"]


def save_fig(fig, name):
    fig.savefig(OUT_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"saved {name}.png/.pdf", flush=True)


# ---------------------------------------------------------------- VAE reconstruction
def fig_vae_reconstruction():
    rows = []
    for ep in EPOCHS_SHOWN:
        for kind in ["original", "augmented"]:
            p = VAE_LOGS / f"plot_model_epoch{ep}_{kind}.png"
            rows.append((f"epoch {ep} — {kind}", np.array(Image.open(p))))
    n_r, n_c = len(rows), 6
    fig, axes = plt.subplots(n_r, n_c, figsize=(15, 2.6 * n_r), dpi=150)
    fig.subplots_adjust(hspace=0.55, wspace=0.12, left=0.115, right=0.99, top=0.93, bottom=0.02)
    for i, (label, arr) in enumerate(rows):
        for j in range(n_c):
            ax = axes[i, j]
            ax.imshow(arr)
            ax.set_xticks([])
            ax.set_yticks([])
        axes[i, 0].set_ylabel(label, fontsize=11, rotation=0, ha="right", va="center", labelpad=90)
    for j, t in enumerate(COLTITLES):
        axes[0, j].set_title(t, fontsize=11)
    fig.suptitle("VAE reconstruction during training (held-out crop; original and augmented views)",
                 fontsize=13, y=0.985)
    save_fig(fig, "figA3a_vae_reconstruction")


# ---------------------------------------------------------------- VAE curves
def parse_vae_log():
    pat = re.compile(
        r"Train Loss: ([\d.]+) \(Recon: ([\d.]+), KL: ([\d.]+), Latent: ([\d.]+)\) \| "
        r"Val Loss: ([\d.]+) \(Recon: ([\d.]+), KL: ([\d.]+), Latent: ([\d.]+)\)")
    data = []
    with open(VAE_LOGS / "training_log.txt") as f:
        for line in f:
            m = pat.search(line)
            if m:
                data.append([float(x) for x in m.groups()])
    return np.array(data)


def fig_vae_curves():
    d = parse_vae_log()
    epochs = np.arange(1, len(d) + 1)
    cols = ["total", "recon", "KL", "latent consistency"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.4), dpi=150)
    fig.subplots_adjust(wspace=0.28, left=0.06, right=0.985, top=0.88, bottom=0.16)
    for j, ax in enumerate(axes):
        tr, va = d[:, j], d[:, j + 4]
        ax.plot(epochs, tr, "tab:blue", lw=2, label="train")
        ax.plot(epochs, va, "tab:orange", lw=2, label="val")
        ax.set_title(cols[j], fontsize=12)
        ax.set_xlabel("epoch")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("loss")
    axes[0].legend(fontsize=10)
    fig.suptitle("VAE training (144,911 train / 7,945 val augmented views, 30 epochs)", fontsize=13)
    save_fig(fig, "figA3b_vae_training_curves")


# ---------------------------------------------------------------- classifier curves
def parse_cls_log():
    nb = json.load(open(CLS_DIR / "3_train_fibermatcher.ipynb"))
    pat = re.compile(
        r"Epoch (\d+), loss_train=([\d.]+), loss_val=([\d.]+), acc_val=([\d.]+), "
        r"precision_val=([\d.]+), recall_val=([\d.]+), F1_val=([\d.]+)")
    rows = []
    for c in nb["cells"]:
        if c["cell_type"] != "code":
            continue
        for o in c.get("outputs", []):
            for line in "".join(o.get("text", [])).splitlines():
                m = pat.search(line)
                if m:
                    rows.append([int(m.group(1))] + [float(x) for x in m.groups()[1:]])
    return np.array(rows)


def fig_cls_curves():
    d = parse_cls_log()
    ep = d[:, 0]
    best_i = int(np.argmax(d[:, 6]))
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), dpi=150)
    fig.subplots_adjust(wspace=0.25, left=0.07, right=0.985, top=0.86, bottom=0.16)
    ax = axes[0]
    ax.plot(ep, d[:, 1], "tab:blue", lw=2, label="train")
    ax.plot(ep, d[:, 2], "tab:orange", lw=2, label="val")
    ax.set_title("BCE loss", fontsize=12)
    ax.set_xlabel("epoch")
    ax.set_ylabel("loss")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=10)
    ax = axes[1]
    ax.plot(ep, d[:, 3], "tab:gray", lw=1.5, ls="--", label="acc")
    ax.plot(ep, d[:, 4], "tab:green", lw=1.5, ls="--", label="precision")
    ax.plot(ep, d[:, 5], "tab:purple", lw=1.5, ls="--", label="recall")
    ax.plot(ep, d[:, 6], "tab:blue", lw=2.5, label="F1")
    ax.plot(ep[best_i], d[best_i, 6], "o", color="tab:red", ms=9, zorder=5)
    ax.annotate(f"best: epoch {int(ep[best_i])}, F1 = {d[best_i, 6]:.4f}",
                xy=(ep[best_i], d[best_i, 6]), xytext=(ep[best_i] + 2.5, d[best_i, 6] - 0.02),
                arrowprops=dict(arrowstyle="->", color="tab:red"), color="tab:red", fontsize=10)
    ax.set_title("validation metrics", fontsize=12)
    ax.set_xlabel("epoch")
    ax.set_ylabel("score")
    ax.set_ylim(0.85, 1.0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, ncol=2)
    fig.suptitle("Same-fibre classifier training (early stopping, patience 10)", fontsize=13)
    save_fig(fig, "figA4a_classifier_training_curves")


# ---------------------------------------------------------------- shared: pairs + embeddings
def load_pairs(name):
    with open(CLS_DIR / name, "rb") as f:
        return pickle.load(f)


def emb_path(pair_path):
    return str(CLS_DIR / "embed" / (os.path.basename(pair_path).replace(".npz", ".npy")))


def load_pair_embeddings(pairs):
    e1 = np.stack([np.load(emb_path(p[0])) for p in pairs])
    e2 = np.stack([np.load(emb_path(p[1])) for p in pairs])
    y = np.array([p[2] for p in pairs])
    return e1, e2, y


def run_classifier(e1, e2):
    import torch
    sys.path.insert(0, "/DATA/F2FMatcher/src")
    from f2fmatcher.classifier.model import PairClassifier
    model = PairClassifier(embedding_dim=256)
    model.load_state_dict(torch.load(CKPT, map_location="cpu"))
    model.eval()
    with torch.no_grad():
        s1 = torch.from_numpy(e1).float()
        s2 = torch.from_numpy(e2).float()
        scores = np.concatenate([
            model(s1[i:i + 4096], s2[i:i + 4096]).squeeze().numpy()
            for i in range(0, len(e1), 4096)
        ])
    return scores


# ---------------------------------------------------------------- classifier scores
def fig_cls_scores():
    pairs = load_pairs("dataset_pair_test")
    e1, e2, y = load_pair_embeddings(pairs)
    scores = run_classifier(e1, e2)
    pos, neg = scores[y == 1], scores[y == 0]
    from sklearn.metrics import roc_auc_score, roc_curve
    auc = roc_auc_score(y, scores)
    fpr, tpr, _ = roc_curve(y, scores)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    fig.subplots_adjust(wspace=0.28, left=0.07, right=0.985, top=0.86, bottom=0.17)
    ax = axes[0]
    lo, hi = 0.0, 1.0
    for s, c, lab in [(neg, "tab:orange", f"different fibre (n={len(neg):,})"),
                      (pos, "tab:blue", f"same fibre (n={len(pos):,})")]:
        ax.hist(s, bins=60, range=(lo, hi), alpha=0.55, color=c, label=lab, density=True)
    ax.axvline(0.5, color="k", ls=":", lw=1.2)
    ax.text(0.51, ax.get_ylim()[1] * 0.9, "threshold 0.5", fontsize=9)
    ax.set_title(f"same-fibre score, held-out test pairs (AUC = {auc:.3f})", fontsize=12)
    ax.set_xlabel("classifier score")
    ax.set_ylabel("density")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3, axis="y")
    ax = axes[1]
    ax.plot(fpr, tpr, "tab:blue", lw=2.2)
    ax.plot([0, 1], [0, 1], "k:", lw=1)
    ax.fill_between(fpr, tpr, alpha=0.08, color="tab:blue")
    ax.set_title(f"ROC (AUC = {auc:.3f})", fontsize=12)
    ax.set_xlabel("1 − specificity")
    ax.set_ylabel("sensitivity")
    ax.grid(alpha=0.3)
    fig.suptitle("Classifier decision on held-out serial-section pairs "
                 f"({len(pairs):,} pairs, 4 sections, 2 samples)", fontsize=13)
    save_fig(fig, "figA4b_classifier_score_distribution")
    print(f"  AUC={auc:.4f}, n_pos={len(pos)}, n_neg={len(neg)}")


# ---------------------------------------------------------------- latent UMAP
def build_fibre_classes(pairs, roi_names):
    """Union-find over label-1 pairs -> fibre identity per ROI."""
    parent = {r: r for r in roi_names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for p in pairs:
        if p[2] != 1:
            continue
        r1 = os.path.basename(p[0]).replace(".npz", "").rsplit("_aug", 1)[0]
        r2 = os.path.basename(p[1]).replace(".npz", "").rsplit("_aug", 1)[0]
        if r1 in parent and r2 in parent:
            union(r1, r2)
    cls = {}
    for r in roi_names:
        cls[r] = find(r)
    return cls


def fig_latent_umap():
    import umap
    pairs = load_pairs("dataset_pair_test")
    rois = set()
    for p in pairs:
        for path in (p[0], p[1]):
            rois.add(os.path.basename(path).replace(".npz", "").rsplit("_aug", 1)[0])
    rois = sorted(rois)
    X = np.stack([np.load(CLS_DIR / "embed" / (r + ".npy")) for r in rois])
    cls = build_fibre_classes(pairs, rois)
    classes = sorted(set(cls.values()))
    cmap = {c: i for i, c in enumerate(classes)}
    colors = np.array([cmap[cls[r]] for r in rois])
    n_cls = len(classes)

    # latent distances: same-fibre vs random pairs
    idx = {r: i for i, r in enumerate(rois)}
    same_d, rand_d = [], []
    pos_pairs = [(os.path.basename(p[0]).replace(".npz", "").rsplit("_aug", 1)[0],
                  os.path.basename(p[1]).replace(".npz", "").rsplit("_aug", 1)[0])
                 for p in pairs if p[2] == 1]
    for a, b in pos_pairs:
        if a in idx and b in idx and a != b:
            same_d.append(np.linalg.norm(X[idx[a]] - X[idx[b]]))
    rng = np.random.default_rng(42)
    n = len(same_d)
    for _ in range(n):
        i, j = rng.integers(0, len(rois), 2)
        if i != j:
            rand_d.append(np.linalg.norm(X[i] - X[j]))
    same_d, rand_d = np.array(same_d), np.array(rand_d)

    emb = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42).fit_transform(X)

    fig = plt.figure(figsize=(15, 4.8), dpi=150)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.85, 1.0], wspace=0.32)
    ax = fig.add_subplot(gs[0])
    n_unmatched = sum(1 for r in rois if sum(1 for p in pos_pairs if r in p) == 0)
    cc = plt.cm.tab20(colors % 20) if n_cls <= 20 else plt.cm.hsv(colors / max(n_cls, 1))
    ax.scatter(emb[:, 0], emb[:, 1], c=cc, s=90, edgecolors="k", linewidths=0.5, alpha=0.9)
    ax.set_title(f"UMAP of 256-d latents — {len(rois)} ROIs, 4 held-out sections\n"
                 f"colour = fibre identity ({n_cls} fibres, {n_unmatched} unmatched)", fontsize=11)
    ax.set_xlabel("umap-1")
    ax.set_ylabel("umap-2")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1])
    # 1-NN class accuracy: how often is the nearest latent the annotated same fibre?
    from scipy.spatial.distance import cdist
    Dnn = cdist(X, X)
    np.fill_diagonal(Dnn, np.inf)
    nn_same, chance = 0, 0.0
    for i, r in enumerate(rois):
        j = int(np.argmin(Dnn[i]))
        if cls[r] == cls[rois[j]]:
            nn_same += 1
        s = sum(1 for r2 in rois if cls[r2] == cls[r])
        chance += (s - 1) / (len(rois) - 1)
    chance /= len(rois)
    obs = nn_same / len(rois)
    bars = ax.bar(["observed", "chance"], [100 * obs, 100 * chance],
                  color=["tab:blue", "tab:orange"], alpha=0.7, width=0.5)
    for b, v in zip(bars, [100 * obs, 100 * chance]):
        ax.text(b.get_x() + b.get_width() / 2, v + 1, f"{v:.1f}%",
                ha="center", fontsize=11, fontweight="bold")
    ax.set_ylim(0, 60)
    ax.set_title(f"nearest latent is the annotated same fibre\n"
                 f"({nn_same}/{len(rois)} ROIs; median ‖Δμ‖ same {np.median(same_d):.2f} "
                 f"vs random {np.median(rand_d):.2f})", fontsize=11)
    ax.set_ylabel("% of ROIs")
    ax.grid(alpha=0.3, axis="y")

    # example: one fibre seen in >= 2 sections
    ax = fig.add_subplot(gs[2])
    ax.axis("off")
    best = None
    for c in classes:
        members = [r for r in rois if cls[r] == c]
        imgs = {r.split("_roi_")[0] for r in members}
        if len(imgs) >= 2 and (best is None or len(members) > len(best[1])):
            best = (c, members, imgs)
    if best is not None:
        c, members, imgs = best
        ex = sorted(members)[:2]
        ims = []
        for r in ex:
            z = np.load(CLS_DIR / "npz_256" / (r + ".npz"))
            fx = z["flow_x"]
            lo, hi = np.percentile(fx, [1, 99])
            ims.append((r, np.clip((fx - lo) / (hi - lo), 0, 1)))
        w, h = ims[0][1].shape
        canvas = np.zeros((h, 2 * w + 8, 3))
        for k, (_, arr) in enumerate(ims):
            canvas[:, k * (w + 8):k * (w + 8) + w] = arr[..., None] * np.array([0.83, 0.66, 0.27])
        ax.imshow(canvas)
        for k, (r, _) in enumerate(ims):
            img_name = r.split("_roi_")[0]
            ax.text(k * (w + 8) + w // 2, h + 10, img_name, ha="center", fontsize=9)
        ax.set_title(f"same fibre (flow_x input) in two sections\n"
                     f"({len(members)} ROIs of this fibre across {len(imgs)} sections)", fontsize=11)
    fig.suptitle("Stain-invariant latent space: the same fibre lands together regardless of section",
                 fontsize=13, y=1.02)
    save_fig(fig, "figA3c_vae_latent_umap")
    print(f"  n_rois={len(rois)}, n_classes={n_cls}, n_unmatched={n_unmatched}, "
          f"median dist same={np.median(same_d):.3f} random={np.median(rand_d):.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="comma list: recon,curves,clscurves,scores,umap")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    which = set(args.only.split(",")) if args.only else {"recon", "curves", "clscurves", "scores", "umap"}
    if "recon" in which:
        fig_vae_reconstruction()
    if "curves" in which:
        fig_vae_curves()
    if "clscurves" in which:
        fig_cls_curves()
    if "scores" in which:
        fig_cls_scores()
    if "umap" in which:
        fig_latent_umap()


if __name__ == "__main__":
    main()
