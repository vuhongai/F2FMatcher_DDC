# Manuscript — completion roadmap

*Snapshot of what exists and what remains to turn the drafted sections into a submittable manuscript.
Updated after §B (evaluation) landed on `main`.*

## Where we are (drafted ✓)

| Piece | File | Figure(s) | Status |
|---|---|---|---|
| §A Method | `Results_A_method.md` | Fig 1 (a pipeline, b worked example, c feature table) | ✓ draft + Methods |
| §B Evaluation | `section_B_evaluation.md` | Fig 2 (B1 coverage, B2 shape-AUC, B2b cross-panel, B3 baselines, B4 final P/R/F1) | ✓ draft; B1–B4 done, B5 pending |
| §C Dystrophic phenotypes | `Results_C_dystrophic_phenotypes.md` | Fig 3–5 | ✓ draft + legends |
| §D Gene-therapy evaluation | `Results_D_gene_therapy.md` | Fig 6–8 | ✓ draft + legends |

**Headline results already in hand:** F2FMatcher F1 = 0.913 vs best baseline 0.381 (§B); six reproducible
fibre states + fibre-intrinsic mdx signature (§C); LICA1 ≫ AAV9 (3%/28%/65% reach WT), laminin scaffold +
regeneration two-factor model, Cls5 architectural resistance, residual galectin-3 (§D).

---

## What remains — prioritised

### Tier 1 — Writing I can do now (no new compute; biggest gaps)
1. **Title + Abstract** — one-paragraph abstract covering method → evaluation → DMD phenotypes → therapy
   ranking. *(→ `MANUSCRIPT.md`)*
2. **Introduction** — the multiplex-serial-section problem; why instance correspondence (not registration/
   keypoint matching); DMD + micro-dystrophin gene therapy context; what F2FMatcher enables; contributions.
3. **Discussion** — interpretation + novelty (two-factor restoration model; molecular≠morphological rescue;
   AAV9-worse-necrosis; residual lysosomal damage reproducing Jaber 2025), limitations, translational
   outlook (laminin/ECM support, oxidative conversion, early treatment, lysosome-directed combos), method
   generalisation beyond muscle.
4. **Consolidated Methods** — merge the per-section Methods (mapping+features from §A; evaluation from §B;
   clustering+stats from §C/§D) into one Methods section with a reproducibility statement.
5. **References / bibliography** — assemble one list: the 8 verified biology refs + §B tool citations
   (Cellpose, DINOv2, LoFTR, RoMa, SuperGlue, LightGlue/SuperPoint, image-matching-models/VisMatch,
   scikit-image, UMAP, scikit-learn). Verify the tool citations.
6. **Assemble `MANUSCRIPT.md`** — single master document stitching Abstract → Intro → Results §A–§D →
   Discussion → Methods → References, with consistent figure numbering and cross-references.
7. **Turn `section_B_evaluation.md` into clean §B Results prose** in the A/C/D manuscript voice (it is
   currently a working draft with tables) and **compose Figure 2** as one multi-panel figure
   (coverage · shape-consistency AUC · cross-panel cycle · final P/R/F1 bar · ablation).

### Tier 2 — Experiments to finish §B (GPU / code; the PI is running some)
8. **B5 robustness** — accuracy vs inter-section distortion, fibre density, segmentation error (synthetic
   perturbation of one side). *Deliverable: Fig 2 robustness panel.*
9. **10×-HE ↔ 5×-HE near-GT** — same stain, two magnifications: a label-free accuracy anchor on full WT
   sections (pair dirs exist, never matched).
10. **Elastic-registration baseline** (ANTs/bUnwarpJ/SIFT-flow + NN) — a stronger geometric baseline than the
    global affine, to quantify the local-distortion argument.
11. **Fix the 15 small-crop crashes** so all 38 GT pairs score (currently 23/38; NaN-safe guard applied,
    propagation edge cases remain).

### Tier 3 — Analyses that strengthen the biology (optional, high value)
12. **TA cross-muscle replication** — run Q1–Q7 on the TA cohort (data present, 19 samples, 663 features);
    show the six states, LICA1≫AAV9, laminin scaffold and regeneration penalty replicate. *Strong for review.*
13. **Clustering robustness** — BIC/silhouette vs k; bootstrap/seed stability; ch3 vs ch9 vs ch18 presets.
14. **Mixed-effects statistics** — fibre nested in animal, for the headline claims (rebuttal-proofing).
15. **Spatial / niche analysis** — do necrotic/regenerating/low-transgene fibres form focal lesions;
    macrophage adjacency. (Centroids available in the mapping output.)

### Tier 4 — Figure & repo housekeeping
16. **Final figure assembly** — composite multi-panel main figures with panel letters, unified style/DPI;
    consolidate §B figures (`visualizations/eval/`) with §A/§C/§D figures (`manuscript/figures/`).
17. **Resolve duplicate notebooks** — `notebooks/Q2_AAV9_vs_LICA1.ipynb` (this work) vs
    `notebooks/Q2_treatment_effect_AAV9_vs_LICA1.ipynb` (independent Qwen analysis): keep one, or keep both
    clearly labelled; likewise decide on `visualizations/Q1`,`Q2` (Qwen) vs `manuscript/figures`.
18. **Author list, affiliations, significance statement, cover letter; Supplementary plan** (supp figures/
    tables: architecture & training curves, per-pair evaluation tables, coverage full-cohort, robustness).

---

## Suggested next step
Tier 1 (items 1–7) is pure writing on results already in hand and produces a **complete first-draft
manuscript**; Tier 2 finishes §B (some GPU-gated, in progress by the PI); Tiers 3–4 strengthen and polish.
Recommended order: **1–7 now → integrate Tier 2 results as they land → 12 (TA) if time → 16–18 for submission.**
