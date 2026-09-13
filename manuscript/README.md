# manuscript/

Planning materials for the F2FMatcher publication(s).

- **`OUTLINE.md`** — full argument structure (4 themes → sections §A–§D), figure list per section,
  a prioritised list of suggested additional analyses (including the method benchmark vs DINOv2 /
  `image-matching-models`, and TA cross-muscle replication), a novelty summary, and verified references.

## Provenance of the analysis figures (already generated & embedded in the notebooks)

| Section | Notebook(s) | Key figures (embedded in the .ipynb) |
|---|---|---|
| §C Dystrophic phenotypes | `Q1_mdx_vs_WT`, `Q5_Q7_fibre_states` | UMAP/cluster identity, mdx–WT effects, within-cluster, lost correlations; regeneration niche, necrosis, lysosomal |
| §D Gene-therapy evaluation | `Q2_AAV9_vs_LICA1`, `Q3_restoration_mechanisms`, `Q4_centronuclei_laminin_review` | composition shift, dystrophin restoration, laminin scaffold, regeneration penalty, Cls5 paradox |

§A (method schematic) and §B (evaluation/benchmark) are **to be produced** — §B1 coverage numbers are
already computed from `results/QUA/mapping_summary.csv` (see OUTLINE §B); §B2–B4 (ground truth + baselines +
ablations) are the main remaining methods work.

All analysis notebooks are fully reproducible from the committed precomputed cluster cache (no re-clustering;
deterministic GMM seed=42) — see `notebooks/README.md`.
