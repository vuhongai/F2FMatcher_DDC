# QUA deep-dive notebooks

Single-fiber analysis of the Quadriceps (QUA) samples of study 22-082, answering the
biological questions on top of the **precomputed, unbiased fiber clustering**
(morphology + HE, `ch3` preset).

| Notebook | Question |
|---|---|
| `Q1_mdx_vs_WT.ipynb` | Disease signature: **mdx (DMD model) vs WT** — globally, per morphological cluster, and in marker–marker correlation structure. |
| `Q2_AAV9_vs_LICA1.ipynb` | Treatment effect: **AAV9 vs LICA1** micro-dystrophin gene therapy (same dose) — cluster-composition shift, sarcolemmal dystrophin restoration, inflammation/necrosis, lysosomal burden, and restored correlations. |
| `Q3_restoration_mechanisms.ipynb` | **Mechanisms**: (1) what underlies restoration in responsive clusters — basal-lamina/laminin scaffold (dominant), oxidative identity, myotropic capsid; (2) what defines the resistant Cls5 — fast-glycolytic IIb remodeled architecture, dissociating molecular (dystrophin restored) from morphological rescue. Plus a **literature synthesis** (confirmed vs newly discovered) with verified PubMed citations. |
| `Q4_centronuclei_laminin_review.ipynb` | **Regeneration & scaffold** (English): (A) cross-review of the independent Qwen analysis on `main`; (B) centronucleated (regenerated) fibres — prevalence by group/cluster, and the finding that regenerated fibres carry **less transgene dystrophin** (a two-factor model: laminin scaffold × regeneration/vector-retention); (C) laminin is **dystrophin-independent** (permissive scaffold, not restored by dystrophin). Verified PubMed citations. |

## Reproducibility design
- **No re-clustering.** Both notebooks load the committed UMAP cache
  `results/QUA/clustering_cache/ch3_zscore_pca30_mindist0.0_nn15.npz` (PCA+UMAP built by
  `scripts/cluster_fibers.py`) and only reproduce the **deterministic GMM k=6 labelling**
  (`random_state=42`), exactly as `scripts/plot_umap.py` does. Clusters are byte-stable
  across runs.
- **6 clusters (Cls1–Cls6)** are named biologically in Q1 §3. Framework used in Q2:
  *Cls2/3/4/6 = treatment-responsive, Cls5 = treatment-resistant.*
- **Statistics at the animal level** (per-animal medians → Mann–Whitney U + Cohen's d;
  n = 5 WT / 5 mdx / 4 AAV9 / 5 LICA1). Only **observed** (matched) marker values are
  used; NaNs are never imputed for statistics.
- **Data-quality caveats enforced:** slide-8 (LAMP2/LGALS3/SQSTM1) **excludes AAV9**
  (mapping ≈16%); slide-7 (COX) for LICA1 (~51%) is flagged. Slide-1 (Dystrophin, anchor)
  and slide-6 (HE) are 100% covered.

## How to run
```bash
pip install numpy pandas scikit-learn scipy matplotlib seaborn umap-learn statsmodels \
            jupyter nbconvert ipykernel
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=900 notebooks/Q1_mdx_vs_WT.ipynb
```
The first run reads the combined-feature pickles (`results/QUA/features_combined/`, ~1 GB)
and writes a convenience frame `results/QUA/clustering_cache/analysis_frame_ch3_k6.npz`
(~315 MB, git-ignored); later runs reload it in seconds. Figures are embedded in the
executed `.ipynb`; loose PNGs written during a run are git-ignored.

The setup cell auto-detects the repo root, so the notebooks run from `notebooks/`, the
repo root, or anywhere the repo is checked out (no `/DATA` symlink needed).
