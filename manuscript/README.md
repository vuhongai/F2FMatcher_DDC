# manuscript/

Planning materials for the F2FMatcher publication(s).

- **`OUTLINE.md`** — full argument structure (4 themes → sections §A–§D), figure list per section,
  a prioritised list of suggested additional analyses (including the method benchmark vs DINOv2 /
  `image-matching-models`, and TA cross-muscle replication), a novelty summary, and verified references.
- **`section_B_evaluation.md`** — **§B working draft (WT-only)**: B1 coverage + B2 self-supervised
  accuracy (shape consistency, cross-panel cycle-consistency) + B3 geometric baselines, all with
  computed numbers, figures, methods, and a prioritised remaining-work list (expert GT, DINOv2 /
  keypoint matchers, elastic registration, ablations, robustness). Reproducible via
  `scripts/precompute_fiber_morphology.py` → `scripts/evaluate_section_B.py` → `scripts/plot_section_B.py`.

## Provenance of the analysis figures (already generated & embedded in the notebooks)

| Section | Source | Key figures |
|---|---|---|
| §B Evaluation (WT) | `scripts/evaluate_section_B.py` → `visualizations/eval/` | Fig 2a coverage; Fig 2 shape-consistency AUC; Fig 2d cross-panel cycle-consistency; Fig 2b baselines |
| §C Dystrophic phenotypes | `Q1_mdx_vs_WT`, `Q5_Q7_fibre_states` | UMAP/cluster identity, mdx–WT effects, within-cluster, lost correlations; regeneration niche, necrosis, lysosomal |
| §D Gene-therapy evaluation | `Q2_AAV9_vs_LICA1`, `Q3_restoration_mechanisms`, `Q4_centronuclei_laminin_review` | composition shift, dystrophin restoration, laminin scaffold, regeneration penalty, Cls5 paradox |

§A (method schematic) is **to be produced**. §B (evaluation/benchmark) has a **WT-only working draft**
(`section_B_evaluation.md`): B1 coverage, B2 self-supervised accuracy, and B3 geometric baselines are
**done**; B3 learned baselines (DINOv2 / keypoint matchers), ablations, and B4 robustness remain (see
*Remaining work* in `section_B_evaluation.md`).

All analysis notebooks are fully reproducible from the committed precomputed cluster cache (no re-clustering;
deterministic GMM seed=42) — see `notebooks/README.md`.
