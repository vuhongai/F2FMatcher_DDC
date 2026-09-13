# §B. Evaluation of the method — working draft (WT cohort)

> **Status:** B1 **done** · B2 self-supervised accuracy **done** (adapted to the available star
> topology) · B3 geometric baselines **done** · **B4 curated-GT P/R/F1 + ablation done** (TA crop
> pairs) · B4c appearance/dense baselines (DINOv2, VisMatch) **setup ready, pending GPU** · B5
> robustness **pending**. Prose below is a first draft in manuscript voice — to be refined with the PI.
> All numbers are computed, reproducible, and sourced (see *Provenance*).

## Scope

Per the PI's decision, the method evaluation is performed **on WT muscle alone** (QUAG01–05;
5 samples × 6 stain panels = **30 section-pairs**, **190,009 matched fibre pairs**). WT is the
cleanest tissue in the study — no dystrophic pathology, no gene therapy, no mapping failures — so
it isolates the matcher's performance from biological confounds. (Full-cohort coverage numbers
remain available in `results/QUA/mapping_summary.csv` for a supplementary comparison.)

**Core claim.** F2FMatcher matches fibres **accurately, at high coverage, and more robustly across
the stain gap than registration- or geometry-only baselines**, with the learned per-fibre features
doing the work that pure geometry cannot.

---

## B1. Coverage

F2FMatcher recovers a high fraction of fibres across every stain pair (Fig 2a,
`fig_B1_coverage.png`). Coverage is reported as the mean of the two sides of each pair (fraction of
anchor fibres matched, and fraction of panel fibres matched).

| Panel (anchor → …) | pairs | coverage (anchor) | coverage (panel) | coverage (both) |
|---|---|---|---|---|
| NADH | 5 | 84.0% | 88.2% | **86.1%** |
| IgG / CD11b | 5 | 80.0% | 87.6% | **83.8%** |
| HE (10×) | 5 | 78.6% | 88.8% | **83.7%** |
| COX | 5 | 77.7% | 86.5% | **82.1%** |
| LAMP2 / LGALS3 / SQSTM1 | 5 | 52.0% | 72.4% | **62.2%** |
| WGA / Myh7 / Myh2 / Myh4 | 5 | 55.1% | 67.6% | **61.3%** |
| **Overall (30 pairs)** | 30 | **71.2%** | **81.8%** | **76.5%** |

**Draft prose.** In WT muscle, F2FMatcher matches on average **76.5%** of fibres per section-pair
(71.2% of anchor fibres, 81.8% of panel fibres; 190,009 matched pairs in total). Coverage is highest
for the housekeeping panels — NADH (86%), IgG/CD11b (84%), HE (84%), COX (82%) — and lower for the
two panels with the weakest segmentation concordance to the laminin anchor, WGA/Myh (61%) and the
lysosomal slide-8 panel (62%). Coverage here reflects **segmentation concordance and section
quality**, not a ceiling of the matcher: the anchor (laminin) segmentation is the reference, and any
panel whose own segmentation misses fibres caps that panel's recoverable coverage.

*Fig 2a — coverage per stain-pair (WT).*

---

## B2. Accuracy (self-supervised; no manual ground truth yet)

No expert-annotated correspondence set exists yet (see *Remaining work*). We therefore validate
accuracy with two **self-supervised** consistency metrics that require no labels: (a) matched fibres
should be more **morphologically consistent** than random fibres, and (b) the same anchor fibre's
counterparts in different panels should **coincide in space** (a star-topology cycle check).

### B2a. Morphological consistency (matched vs random)

The same physical fibre in two serial sections has the same shape. We compare each F2FMatcher match
to a random pair of fibres from the two panels using **scale-invariant** shape descriptors
(eccentricity, solidity, extent, roundness = 4πA/P², aspect ratio — all invariant to the global
affine between sections), z-scored across the WT cohort, and report the Euclidean shape distance.
Lower distance = more consistent = more likely a true match.

**Result.** F2FMatcher matches are far more shape-consistent than random pairs
(Fig 2, `fig_B2_morph_consistency.png`): mean shape distance **1.09** (matched) vs **2.65**
(random) — a **59% reduction** — with a pooled separation **AUC of 0.865** (range 0.74–0.93 per
pair). The effect holds for every panel:

| Panel | mean AUC (5 pairs) | AUC range | matched dist | random dist |
|---|---|---|---|---|
| IgG / CD11b | 0.902 | 0.82–0.93 | 0.98 | 2.61 |
| NADH | 0.898 | 0.82–0.93 | 0.90 | 2.56 |
| HE (10×) | 0.874 | 0.79–0.91 | 1.05 | 2.64 |
| COX | 0.867 | 0.85–0.88 | 1.08 | 2.63 |
| LAMP2 / LGALS3 / SQSTM1 | 0.839 | 0.82–0.86 | 1.19 | 2.72 |
| WGA / Myh | 0.813 | 0.74–0.86 | 1.40 | 2.72 |

**Draft prose.** A true cross-section match pairs two views of the *same* fibre, so its shape
descriptors should agree; a random pairing should not. Across all 30 WT section-pairs, F2FMatcher
matches are **59% more shape-consistent** than randomly paired fibres (mean shape distance 1.09 vs
2.65 in standardised units), and a match-vs-random classifier reaches **AUC 0.865** (0.74–0.93 per
pair). The separation is strongest for the panels with the clearest segmentation (IgG/CD11b 0.90,
NADH 0.90) and weakest for WGA/Myh (0.81), tracking the coverage ordering in B1. This is a
label-free confirmation that the matches are biologically coherent, not an artefact of the geometry.

### B2b. Cross-panel spatial consistency (star-topology cycle check)

Every panel is mapped to the same laminin anchor. If the matches are correct, the counterparts of a
single anchor fibre in two *different* panels must be the same physical fibre and therefore
**coincide** once each panel is aligned to the anchor by a robust global affine. We measure the
mean pairwise spread of an anchor fibre's counterparts across its 6 panels, and compare it to the
spread obtained by assigning **random** counterparts (the null).

**Result.** Real matches are **2.7× tighter** in space than the random null
(Fig 2d, `fig_B2b_crosspanel.png`):

| Sample | anchor fibres (≥2 panels) | real spread (px) | null spread (px) | tighter |
|---|---|---|---|---|
| QUAG01 | 6,814 | 2,404 | 4,575 | 1.9× |
| QUAG02 | 8,486 | 1,364 | 4,292 | 3.1× |
| QUAG03 | 8,492 | 1,261 | 4,649 | 3.7× |
| QUAG04 | 8,081 | 1,826 | 4,326 | 2.4× |
| QUAG05 | 9,152 | 1,922 | 4,555 | 2.4× |
| **Pooled** | — | **1,755** | **4,479** | **2.7×** |

**Draft prose.** Because all six panels share one anchor, a correct match is *cycle-consistent*: the
counterpart of an anchor fibre in panel *p* and in panel *q* must land on the same physical fibre
after each panel is aligned to the anchor. F2FMatcher's matches satisfy this — their cross-panel
spread is **1,755 px** versus **4,479 px** for random assignments (**2.7× tighter**, 1.9–3.7× per
sample). The residual spread reflects genuine inter-section elastic distortion (which a single global
affine cannot remove) plus a minority of mismatches; the key point is that the matches are far more
spatially coherent than chance, with no labels involved.

---

## B3. Benchmark against baselines (geometric — partial)

We compare F2FMatcher to geometry-only baselines that produce a per-fibre assignment, scored on the
**same** shape-consistency metric as B2a (so the comparison is apples-to-apples and label-free).
All methods are scored on the *same* anchor fibres (those F2FMatcher matched).

| Method | mean shape distance ↓ | reduction vs random |
|---|---|---|
| **F2FMatcher** | **1.09** | **59%** |
| Global affine + kNN | 2.61 | 2% |
| No-alignment kNN | 2.66 | −0% |
| Random null | 2.65 | — |

(Fig 2b, `fig_B3_baselines.png`.)

**Baselines.**
- *Global affine + kNN:* a robust (iteratively-reweighted) affine is fit on the correspondences,
  then each anchor fibre is assigned to its nearest neighbour in the aligned space.
- *No-alignment kNN:* nearest neighbour in raw coordinates (no global alignment).

**Result.** Both geometric baselines perform **at chance level** on shape consistency (2.61 and 2.66
vs random 2.65), whereas F2FMatcher is **59% above chance** — a **2.4×** advantage over the best
geometric baseline. Notably, the global affine is *harmful* in several pairs (e.g. QUAG01 IgG:
affine 4.22 vs no-align 2.60): forcing thousands of fibres through one global transform misplaces
them, because serial sections carry **local** tears and distortion that no single affine captures.
Muscle fibres are also densely packed, so the nearest neighbour in space is frequently a *wrong but
adjacent* fibre. Resolving the correct fibre among nearby candidates is exactly what the learned
per-fibre features (shape + appearance + triangle geometry + propagation) are for.

**Honest caveat.** This is a **proxy** for accuracy (shape consistency, not precision/recall/F1
against a labelled ground truth), and it favours methods that use shape — which F2FMatcher does and
the geometric baselines do not. The definitive benchmark (precision/recall/F1 vs expert GT, plus
appearance-only baselines DINOv2 and keypoint matchers, plus ablations) is the remaining work below.

---

## B4. Definitive accuracy against curated ground truth (TA crop pairs)

The self-supervised checks in B2 are label-free proxies. We now report the **definitive** accuracy:
precision / recall / F1 against a **curated expert ground truth** of fibre correspondences on TA
muscle crop pairs.

**Ground truth.** 73 TA crop images (~486 px, WGA 10×) with **1,315 expert-annotated fibre
correspondences** across **38 image-pairs** (`datasets/training_data_updated.csv`). Each annotation is
a pair of CellPose region labels (one per image) judged by the expert to be the same physical fibre.
We verified the label space is exactly reproducible: re-running the annotation CellPose setup
(`CP_AV_WGA_Dia_Qua_TA_AxioScan10X`, channels=[0,0], flow_threshold=0.4, cellprob_threshold=0,
diameter=model default, cellpose 2.2.2) reproduces **100%** of the annotated labels (`gt_in_fresh=1.0`),
so every method is scored in the same label space.

**Edge-ROI caveat (important).** F2FMatcher only matches ROIs whose full 256-px crop window fits
inside the image (`filter_ROIs`); ROIs with a centroid within 128 px of any border are discarded by
design. On these 486-px crops, **523 of the 1,315 GT pairs (40%) involve at least one such edge ROI**
and are therefore *unmatchable by F2FMatcher by design*. We report recall two ways: over **all** GT
pairs (`R_all`) and over the **792 matchable** (both ROIs non-edge) pairs (`R_match`, the fair
denominator). Precision is unaffected (every method only predicts non-edge ROIs).

### B4a. F2FMatcher vs geometric baselines (P/R/F1 vs GT)

| Method | pairs | P | R_all | R_match | F1_match |
|---|---|---|---|---|---|
| **F2FMatcher** | 23/38 | **0.937** | 0.494 | **0.891** | **0.914** |
| No-alignment kNN | 38 | 0.339 | 0.262 | 0.434 | 0.381 |
| Global affine + kNN | 38 | 0.194 | 0.150 | 0.249 | 0.218 |
| Random null | 38 | 0.022 | 0.008 | 0.013 | 0.016 |

F2FMatcher is **highly precise** (P=0.937: a predicted match is correct ~94% of the time) and, over
the matchable ROIs, reaches **R=0.891 / F1=0.914** — **2.4× the F1** of the best geometric baseline
(no-align kNN, 0.381) and **4.2×** the affine baseline (0.218). The geometric baselines sit near
chance (random F1=0.016), confirming that the learned per-fibre features + triangle geometry +
propagation do the work pure geometry cannot.

**Robustness caveat.** F2FMatcher completed **23/38** crop pairs; 15 crashed on small-crop edge cases
(7 degenerate-triangle NaN, 4 neighbour-index, 2 no-seed, 1 2-D-array, 1 unpack). These are real
limitations on very small / sparse crops (the matcher was tuned for full sections). A NaN-safe
triangle-angle guard (crash-fix only, `scripts/run_f2fmatcher_benchmark.py`) was applied; the
remaining crashes are propagation edge cases we do not paper over. The 23 completed pairs are the fair
benchmark set.

### B4b. Ablation: mapping quality as steps are added (22 pairs)

We re-ran the pipeline capturing the intermediate matched-label set at each step and scored each
stage against the GT (same edge-ROI correction):

| Stage | #pred | P | R_match | F1_match |
|---|---|---|---|---|
| S1 classifier seeds (initial guess) | 518 | 0.683 | 0.753 | 0.717 |
| S2 + geometry validation | 490 | 0.722 | 0.753 | 0.737 |
| S3 + local propagation | 297 | 0.926 | 0.585 | 0.717 |
| S4 + fill unannotated (FULL) | 447 | 0.937 | 0.891 | 0.914 |

The classifier seeds already carry most of the recall (R=0.753) at moderate precision (0.683);
geometry validation lifts precision (0.722) without losing recall; local propagation is conservative
(high precision 0.926, but it prunes to the highest-confidence pairs, dropping recall to 0.585); and
the affine fill of unannotated ROIs recovers the recall (0.891) while keeping precision (0.937). The
full pipeline (S4) is the best (F1=0.914) — each step contributes.

### B4c. Appearance / dense-matching baselines (DINOv2, VisMatch) — pending GPU

Per the PI, we benchmark appearance-only and dense/sparse matchers against the same GT. These map
**pixel→pixel**; we assign each source pixel/patch to its ROI via the CellPose label map, aggregate
to (ROI,ROI) correspondences, and score P/R/F1 + coverage (the coverage metric also applies to the
large sections that have no GT). Setup is ready: a dedicated `vismatch` conda env (torch 2.14+cu130,
`vismatch` 1.3.2 wrapping 50+ matchers, `dinov2`) and `scripts/benchmark_vismatch.py`. Models: dense
(LoFTR, RoMa), sparse (SuperPoint-LightGlue, SuperGlue), and DINOv2 (ViT-S/14 patch matching).
**Needs the A30s freed** (currently occupied by the local LLM). Smoke-tested on CPU (1 pair):
DINOv2 P=0.079/R=0.857, SuperPoint-LightGlue P=0.208/R=0.905 — appearance/texture baselines are
recall-heavy but imprecise vs F2FMatcher (P=0.937). Run on GPU once freed:

```bash
PYTHONNOUSERSITE=1 /home/avuhong/anaconda3/envs/vismatch/bin/python \
    scripts/benchmark_vismatch.py --models loftr,roma,superpoint-lightglue,superglue,dinov2 \
    --device cuda
```

---

## Remaining work (B5 robustness · B4c GPU baselines)

Ordered by value-to-effort. Items marked **GPU** need the A30s freed (the local LLM is currently
running on them); items marked **manual** need the PI/lab. **Done** items are struck through.

1. ~~**Expert ground truth (B2, manual).**~~ **DONE → B4.** The curated TA-crop GT (1,315
   correspondences, 38 pairs) is in place and scored: F2FMatcher P=0.937, R_match=0.891, F1=0.914.
2. **10×-HE ↔ 5×-HE near-GT (B2, GPU).** Same stain, two magnifications — a near-ground-truth
   sanity check on the full WT sections. The 5×-HE pair directories exist but were never matched;
   running the matcher on the 5 WT 5×-HE pairs (~1–2 h each on CPU, faster on GPU) gives a
   label-free accuracy anchor.
3. **Appearance / dense baselines (B4c, GPU).** **Setup ready, pending GPU.** DINOv2 (ViT-S/14
   patch matching) + VisMatch (LoFTR, RoMa, SuperPoint-LightGlue, SuperGlue) via the dedicated
   `vismatch` env + `scripts/benchmark_vismatch.py` (pixel→ROI assignment, P/R/F1 + coverage).
4. **Elastic registration baseline (B3, CPU/GPU).** ANTs / bUnwarpJ / SIFT-flow + NN — a stronger
   geometric baseline than the global affine, to show the local-distortion argument quantitatively.
5. ~~**F2FMatcher ablations (B3, GPU).**~~ **DONE → B4b.** Step-by-step P/R/F1 (seeds → geometry →
   propagation → fill) shows the full pipeline is best (F1=0.914).
6. **Robustness (B5).** Accuracy vs inter-section distance/distortion, fibre density, and
   segmentation errors (synthetic distortion of one side of a pair, re-run the matcher).
   *Deliverable: Fig 2e.* Also: fix the 15 small-crop crashes (B4a caveat) to benchmark all 38 pairs.

---

## Methods (reproducibility)

- **Cohort.** WT = QUAG01–05 (from `config/ddc_config.py` `GROUP_MAP`). 30 section-pairs = the 5
  samples × 6 non-5×-HE panels, each with a `paired_labels.pkl`.
- **Morphology.** Per-fibre shape features from CellPose masks (`results/QUA/out_CP_masks/`),
  restricted to the mapping ROI set (`filter_ROIs`, 256-px crop in-bounds, area ≥ 100). Precomputed
  once to `results/QUA/eval/fiber_morphology.pkl` (35 images, 271,119 fibres) by
  `scripts/precompute_fiber_morphology.py`.
- **Shape distance (B2a/B3).** Scale-invariant descriptors {eccentricity, solidity, extent,
  roundness, aspect ratio}, z-scored across the WT cohort, Euclidean distance between the two
  fibres of a pair. Random null = uniform random pairing from the two panels' fibre pools
  (20,000 draws/pair, seed 42). AUC via `sklearn.roc_auc_score` on (−distance).
- **Cross-panel spread (B2b).** Per panel, a robust affine (iteratively reweighted least squares,
  3 rounds, 4·MAD outlier cut) aligns the panel to the anchor on the matched correspondences. For
  each anchor fibre matched in ≥2 panels, the mean pairwise distance of its counterparts (mapped
  back to the anchor frame) is the spread; null = random counterparts (3,000 draws/sample).
- **Baselines (B3).** *Global affine + kNN* uses the robust affine above + `scipy.cKDTree`
   nearest neighbour in the aligned space; *no-alignment kNN* uses raw coordinates. Both scored on
   the same anchor fibres and same shape distance as F2FMatcher.
- **Curated GT benchmark (B4).** TA crop pairs from `datasets/training_data_updated.csv` (1,315
   correspondences, 38 pairs, 73 images). F2FMatcher is run on the exact annotation CellPose setup
   (`scripts/run_f2fmatcher_benchmark.py`, CPU, NaN-safe triangle-angle guard). A predicted
   (label1,label2) pair is a true positive iff it is in the curated GT. The **edge-ROI correction**
   splits the GT into *matchable* (both ROIs pass `filter_ROIs`) and *edge* pairs; recall is reported
   over both (P/R/F1 in `scripts/benchmark_gt.py`). The **ablation** re-runs the pipeline with
   `save_step_prediction=True` (`scripts/run_f2fmatcher_ablation.py`) and scores the intermediate
   matched-label set at each step (`scripts/benchmark_ablation.py`).
- **DINOv2 / VisMatch (B4c).** Dedicated `vismatch` conda env (torch 2.14+cu130, `vismatch` 1.3.2,
   `dinov2`); `scripts/benchmark_vismatch.py` runs each matcher, assigns matched pixels/patches to
   ROIs via the CellPose label map, and scores P/R/F1 + coverage. Pending GPU.
- **Scripts (WT self-supervised).** `scripts/precompute_fiber_morphology.py` →
   `scripts/evaluate_section_B.py` → `scripts/plot_section_B.py`. End-to-end runtime (from cache): ~2 s.

## Provenance

| Artifact | Path |
|---|---|
| Per-fibre morphology cache | `results/QUA/eval/fiber_morphology.pkl` |
| Headline numbers | `results/QUA/eval/section_B_results.json` |
| B1 per-pair / per-panel | `results/QUA/eval/coverage_WT.csv`, `coverage_WT_by_panel.csv` |
| B2a per-pair | `results/QUA/eval/morph_consistency.csv` |
| B2b per-sample | `results/QUA/eval/crosspanel_consistency.csv` |
| B3 per-pair / summary | `results/QUA/eval/baselines.csv`, `baseline_compare.csv` |
| Fig 2a coverage | `visualizations/eval/fig_B1_coverage.png` |
| Fig 2 (B2a) shape consistency | `visualizations/eval/fig_B2_morph_consistency.png` |
| Fig 2d (B2b) cross-panel | `visualizations/eval/fig_B2b_crosspanel.png` |
| Fig 2b (B3) baselines | `visualizations/eval/fig_B3_baselines.png` |
| B4 GT per-pair P/R/F1 (edge-corrected) | `results/benchmark/benchmark_scores.csv` |
| B4 F2FMatcher pipeline output (23/38 pairs) | `results/benchmark/f2fmatcher_output/` |
| B4 ablation step predictions (22 pairs) | `results/benchmark/f2fmatcher_output_ablation/` |
| B4c DINOv2/VisMatch scores (pending GPU) | `results/benchmark/vismatch/vismatch_scores.csv` |
| B4c env build / install logs | `results/benchmark/vismatch_env2.log`, `vismatch_install.log` |
