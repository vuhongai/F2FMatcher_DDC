# Results §A — F2FMatcher: fibre-to-fibre matching across histological stains

*Draft Results/Methods text for theme 1 ("Development of fibre-to-fibre mapping"). This is the opening
section of the paper; it introduces the method and produces the per-fibre multiplex table on which §B
(evaluation) and §C–§D (biology) depend. Figure 1 is the method overview. This version details how each
component is **trained**, the exact **matching steps**, and the **rationale** behind every design and
parameter choice.*

## A.1 The problem: single-fibre correspondence across stains

Skeletal muscle phenotyping is multi-parametric — dystrophin, membrane integrity, immune infiltration,
metabolism, fibre type, autophagy — but no single section can be stained for all of these markers at once.
Markers are instead distributed across **serial sections**, each carrying a different immunofluorescence or
histochemical panel. Relating them requires knowing which fibre in one section corresponds to which fibre in
the next. This is hard for three reasons: (i) the **appearance gap** between modalities (a laminin
immunofluorescence border versus an NADH brightfield field) means the *same* fibre looks completely different
across stains; (ii) serial sections undergo **local tears, folds and non-rigid distortion**, so a single
global transform does not align them; and (iii) a whole muscle section contains **thousands of near-identical
fibres**, so appearance alone is ambiguous. Classical image registration returns a warp field, not fibre
identities; general keypoint matchers and dense visual features are trained for same-modality scenes and do
not assign individual cells across a stain gap. F2FMatcher solves the *instance-correspondence* problem
directly.

**F2FMatcher** links corresponding fibres across pairs of differently-stained serial sections by combining a
**shape-based deep embedding** that is invariant to stain, a **pairwise classifier** that scores whether two
fibres are the same, and a **geometry-aware matching algorithm** that enforces spatial consistency (Fig 1a).
Applied to a seven-panel serial stack of murine quadriceps anchored on one reference section, it produced a
single **per-fibre multiplex table of 807 features for 111,403 fibres** across 19 animals — the substrate for
every downstream analysis in this paper.

![Figure 1a](figures/fig1a_pipeline.png)

## A.2 A stain-invariant input: the Cellpose flow field, not the stain

The design principle that makes cross-stain matching tractable is to **describe each fibre by its shape, not
by its staining intensity**. Every section is segmented with a fine-tuned Cellpose 2 model, which returns for
each fibre both a binary mask and the underlying **flow field** — the 2-D vector field `(flow_x, flow_y)`
that Cellpose integrates to define cell boundaries, pointing from every pixel toward its cell centre. Around
each fibre centroid we cut a fixed **256×256 window** (fibres whose window would fall outside the slide, or
whose area is < 100 px², are excluded) and keep three channels: the two flow components and the binary ROI
mask.

Crucially, the 256×256 window is deliberately larger than a single fibre: the two flow channels are cut from
the **whole-section** flow field, so they retain the flow of the **neighbouring fibres** that fall inside the
window, while the binary mask marks only the central target. The embedding therefore describes each fibre
**together with its immediate surroundings**, not in isolation — a property we exploit below (§A.3).

The flow field is re-encoded before it enters the network. Rather than feed the signed Cartesian components
`(flow_x, flow_y)`, we convert them to **polar form** — magnitude `mag = clip(‖flow‖/10, 0, 1)` and direction
`angle = (atan2(flow_y, flow_x)+π)/2π` — stack them with the mask, and resize the 256×256 window to
**3×128×128**.

> *Rationale (input).* (1) *Stain-invariance:* the flow field is derived from the segmentation geometry, so a
> given fibre yields essentially the same input whether it was stained for laminin or for NADH — the property
> that lets fibres be compared across the modality gap. (2) *Why polar, why bounded:* the training/augmentation
> pipeline uses PIL/torchvision transforms, which require non-negative, bounded channels; the signed,
> unbounded `(flow_x, flow_y)` cannot pass through them, whereas `mag, angle ∈ [0,1]` can. Polar form also
> **disentangles** the descriptor into a size/shape term (magnitude, largest at the centre) and a contour
> term (direction), which are cleaner, common-scaled features for a CNN than two entangled signed components.
> (3) *Why 128:* the flow field is smooth, so down-sampling 256→128 preserves fibre geometry while quartering
> compute.

## A.3 Training the embedding (VAE)

The 3×128×128 input is compressed by a variational auto-encoder (**`SharedMultiHeadVAE`**) into a
**256-dimensional latent vector**. A six-convolution encoder (channels 64-64-128-128-256-256 with two
max-pools and an adaptive 4×4 pool) maps the input to two 256-d heads, the latent mean **μ** and
log-variance; the reparameterization trick draws `z = μ + exp(½·logvar)·ε`. A single shared transposed-
convolution decoder expands `z` back to 128×128 and three linear heads reconstruct `flow_x`, `flow_y` and the
mask (≈15 M parameters total).

**Training.** For every fibre, **50 augmented views** are generated (rotation ≤ 90°, shear ≤ 5°, isotropic
scale ± 10%). The network is trained to reconstruct the flow field and mask under a variational prior:

```
L = Σ MSE(recon_{fx,fy,mask}, target)  +  β_KL · KL(q(z|x) ‖ N(0,I))          (β_KL = 0.001)
      + β_consistency · ‖ μ(view_a) − μ(view_b) ‖   for paired augmentations   (β_consistency = 1.0)
```

optimised with Adam (lr 1e-3), batch 32, up to 20 epochs with early stopping (patience 10). The embedding
used everywhere downstream is the **deterministic mean μ** (the stochastic `z` and the `logvar` head serve
only the KL term during training).

> *Rationale (embedding training).* The reconstruction objective forces the 256-d bottleneck to retain the
> full geometry of the fibre (it must regenerate the entire flow field and mask from `z`), so μ is a faithful
> shape descriptor. The small KL weight (β_KL = 0.001) keeps the latent space smooth and regularised without
> letting the prior wash out shape detail — a reconstruction-dominant VAE. **Invariance is engineered, not
> assumed:** the 50 geometric augmentations per fibre, together with the latent-consistency term that pulls
> augmented views of the *same* fibre to the same μ, make the embedding robust to the rotation, shear and
> scale differences that separate serial sections, so the *same* fibre lands at nearly the same point in
> latent space regardless of section or stain. Reconstructing the **Cartesian** `(flow_x, flow_y)` from a
> **polar** input adds a mild extra constraint (the decoder must learn the polar→Cartesian map) and avoids the
> angle wrap-around discontinuity in the loss. Finally, because the reconstruction targets are the **full**
> flow crop — which contains the neighbouring fibres inside the 256-px window — μ is *forced* to encode the
> target fibre **together with the morphology and arrangement of its neighbours**. The appearance term is thus
> already **context-aware**: it compares not two isolated silhouettes but two fibres each embedded in their
> local surroundings, a far more discriminative signal in a field of near-identical fibres.

## A.4 Training the same-fibre classifier

Shape similarity is turned into a decision by a **pairwise classifier** (`PairClassifier`). It concatenates
two 256-d embeddings (→ 512) and passes them through an MLP `512 → 128 → 64 → 1` (ReLU activations, **sigmoid**
output), returning the probability that the two fibres are the same physical cell. Scoring every cross-section
pair gives an **N₁×N₂ same-fibre score matrix S**.

**Training.** Positive pairs are curated same-fibre correspondences; negatives are sampled at **4× the number
of positives** (hard-negative mining, `negative_fold = 4`). Training uses binary cross-entropy, Adam
(lr 1e-4), batch 256, up to 50 epochs with early stopping on validation F1 (patience 10). The selected
checkpoint reached **validation F1 ≈ 0.944**.

> *Rationale (classifier).* Operating on the *frozen* 256-d embeddings (not on pixels) makes the classifier
> tiny, fast enough to score the full N₁×N₂ grid of a whole section, and dependent only on shape. The **4×
> negative oversampling** reflects the real matching regime — for any fibre there is at most one true partner
> and thousands of impostors — so the classifier is trained where it must operate: rejecting look-alikes. A
> low learning rate (1e-4) and F1-based early stopping tune the operating point for a heavily imbalanced,
> precision-critical task. The score is used as a soft term (and later thresholded), never as the sole
> criterion — geometry (§A.5–A.6) resolves the residual ambiguity between similar fibres that appearance
> cannot.

## A.5 The matching cost: appearance × geometry

Both ingredients of the cost are **local descriptors**: they characterise a fibre only by itself and its
immediate neighbourhood, with no knowledge of the section-wide layout. Appearance similarity alone cannot
separate thousands of near-identical fibres, so the matching cost pairs it with a purely geometric term. For
each cross-section pair `(i, j)` the two ingredients are:

- **Appearance — classifier score `S[i,j]` (§A.4):** does fibre *i* look like the same cell as fibre *j*?
- **Geometry — spatial-signature similarity `G[i,j]`:** does fibre *i* sit in the same local neighbourhood as
  fibre *j*? Each fibre's "geometric niche" is the sorted vector of distances to its **k nearest neighbours**;
  two fibres are compared by the **Wasserstein distance** between these distance distributions,
  `G_k = 1/(1+W)`. This is done at **three scales (k = 3, 5, 7)** and the three similarity matrices are
  combined by geometric mean (`G = (∏_k G_k)^{1/3}`).

> *Rationale (cost).* The two terms are **orthogonal and complementary**: the classifier is invariant to stain
> but confused by look-alike fibres; the spatial signature ignores appearance entirely but is invariant to
> global translation, rotation and scale (it uses only *relative* neighbour distances). Requiring a match to
> satisfy **both** breaks the ambiguity that defeats either alone. The Wasserstein distance compares the
> *distribution* of neighbour distances rather than element-wise values, so it tolerates a neighbour gained or
> lost between sections (segmentation noise); the multiscale `k = 3, 5, 7` captures both the tightest local
> arrangement and a wider context, and the geometric mean makes a pair pay a penalty on **every** scale to
> score well.

## A.6 The matching algorithm: seeds → propagation → affine fill

All three stages accept a candidate correspondence only if it passes the same **triangle-congruence test**:
for two matched fibres plus a nearby matched neighbour, the triangle formed by their centroids in section 1
must be congruent to the corresponding triangle in section 2 — i.e. the mean absolute difference of the three
**side lengths** is below `max_cost_geo_neighbors_sides` (30 px) *and* of the three **angles** below
`max_cost_geo_neighbors_angles` (0.15, ≈ 27°). Side lengths and angles are invariant to translation and
rotation, so a congruent triangle is strong evidence that the correspondences are jointly correct; a single
wrong pair distorts at least one triangle and is rejected.

**Step 1 — triangle-consistent seeds.** The combined cost `(S + G)`, gated by `S > min_cls_logit_init` (0.75),
is scanned greedily: the highest-cost pair is taken, its row and column are removed (each fibre used once),
and this repeats to collect the top `n_initial_guess` (80) mutually-exclusive candidate anchors. Among these,
every combination of `n_pair_selected` (4) candidates is tested and kept only if **all** of its centroid
triangles are congruent across the two sections; the union of surviving candidates forms the seed set. If
fewer than three seeds survive, `n_initial_guess` is increased and the search retried (up to 3 times).

> *Rationale (seeds).* Seeds must be almost certainly correct because everything propagates from them, so the
> gate is deliberately **strict** (appearance threshold 0.75, higher than the 0.5 used later) and correctness
> is verified by *joint* geometry, not just high individual score — testing combinations of four and requiring
> *all* their triangles to agree removes a high-scoring pair that happens to sit in the wrong place. A handful
> of such globally-consistent anchors is enough to bootstrap the rest.

**Step 2 — iterative local propagation.** From each matched pair, its neighbours within a 200-px radius
(`distance_neighbors_ref`) in each section are considered; the best-scoring local candidate is accepted only
if the triangles it forms with its **three nearest already-matched neighbours** are congruent. Accepted
matches become new anchors, and the wave-front grows. After each pass, duplicates are removed and every match
is re-validated against its three nearest matched neighbours (`n_neighbors_validation = 3`); candidates that
repeatedly fail are abandoned after `patience_label` (5) attempts. Iteration stops when a full pass adds fewer
than **0.25 % of fibres** as new matches.

> *Rationale (propagation).* The transform between serial sections is **not globally affine** — tissue
> stretches and tears unevenly across a slide — so each decision is validated only against the *nearest*
> confirmed matches, assuming rigidity only **locally**, which holds even where the global map is non-linear.
> Growing outward from trusted seeds turns matching into constraint propagation on a geometric graph: each new
> match both extends the field and tightens the context for its neighbours. The re-validation pass and the
> patience limit prevent a single early error from contaminating a whole region.

**Step 3 — affine fill.** Fibres still unmatched after propagation (isolated, sparse regions) are recovered
with a coarse global prior: a single affine transform is fitted to all matched centroids, each unmatched
section-1 centroid is projected into section 2, and candidates within `max_distance_affine` (150 px) of the
predicted location with classifier score > `min_cls_logit` (0.5) are accepted in score order **if** they pass
the same triangle-congruence check. A final duplicate-removal and re-validation pass closes the procedure.

> *Rationale (fill).* Propagation cannot reach fibres with no matched neighbours nearby; the global affine is
> only a **search prior** to localise the handful of candidates, never the decision itself — appearance and
> local-triangle geometry still gate every acceptance, so the coarse global fit does not introduce
> geometrically inconsistent matches.

Together the three stages implement a **coarse-to-fine confidence** strategy — a few globally-verified seeds,
a locally-verified propagation wave-front, then an affine-primed clean-up — in which geometry imposes the
global spatial consistency that appearance cannot, while local validation absorbs the non-rigid distortion
between serial sections.

**The two-scale logic of the method.** This is the conceptual division of labour between the cost (§A.5) and
the matcher (§A.6). Everything in the cost is **local** — the appearance embedding sees a fibre and its
immediate surroundings within one 256-px patch, and the spatial signature sees a fibre and its k nearest
neighbours — so on its own the cost can say *"these two fibres look and sit alike locally"* but cannot resolve
which of many locally-similar candidates is the true partner. The matcher supplies the missing **wider scale**:
triangle-verified seeds are consistent across the whole section, propagation stitches the local descriptors
into one coherent field, and the affine fill reasons over a section-wide transform. Local evidence proposes;
global geometry disposes.

## A.7 From a stained stack to a per-fibre multiplex table

Matching every panel to a common anchor section links each anchor fibre to its counterparts throughout the
stack. For every matched fibre, staining intensity is quantified in **four sub-cellular compartments** derived
from the mask by morphological operations — **whole** (full mask), **mem** (a membrane ring), **cyto1**
(peripheral cytoplasm) and **cyto2** (deep core) — capturing, for example, sarcolemmal versus cytoplasmic
versus central-nuclear signal (Fig 1c). Concatenating 15 morphological descriptors with nine intensity
statistics (mean, s.d., p10–p90, skewness, kurtosis) per channel per compartment across the 22 stained
channels yields **807 features per fibre** (15 + 36 × 22). Fibres not matched on a given panel receive NaN
for those channels and are never imputed for statistics. The result is a single table in which every fibre
carries its full multi-marker, sub-compartment phenotype — enabling the marker–marker and state analyses of
§C–§D that are impossible from any one stain.

![Figure 1c](figures/fig1c_features.png)

---

### Figure legend

**Figure 1. F2FMatcher matches muscle fibres across differently-stained serial sections into a per-fibre
multiplex table.**
(a) Pipeline: serial sections are segmented with fine-tuned Cellpose; 256×256 flow-field crops around each
fibre are converted to polar (magnitude, angle, mask), resized to 128×128 and encoded by a VAE into stain-
invariant 256-d embeddings; a pairwise classifier scores same-fibre probability; a geometry-aware algorithm
(cost = classifier × spatial signature; triangle-consistent seeds → iterative propagation → affine fill →
validation) assigns fibre correspondences; per-fibre staining is quantified in four compartments, giving an
807-feature vector per fibre. (b) Worked example: two serial sections stained for laminin (immunofluorescence,
left) and NADH (brightfield, right); Cellpose segmentation (coloured overlays) and F2FMatcher correspondence
link fibres across the modality gap. (c) Each fibre is partitioned into whole/mem/cyto1/cyto2 compartments
(left); the per-fibre feature vector concatenates 15 morphological features with 9 statistics × 4 compartments
for each of 22 channels = 807 features (right).

*Optional additional panels (to add): (i) a UMAP of the VAE latent space coloured by fibre, showing that the
embedding groups the same fibre across stains; (ii) a schematic of Step 3 (triangle-congruence check → seeds →
propagation wave-front → affine fill).*

---

## Methods (§A — training, matching and feature extraction)

**Imaging and segmentation.** Whole-slide CZI images were exported to PNG and resized to a common pixel
resolution (`f2fmatcher.io.czi_reader`). Fibres were segmented with fine-tuned Cellpose 2 models
(`CellPose2_finetuned`; cellprob threshold 0, flow threshold 0.4), producing per-fibre masks and flow fields;
objects below 100 px² and fibres whose 256-px crop window fell outside the slide were discarded.

**VAE embedding and training.** For each fibre a 256×256 crop of (flow_x, flow_y, mask) was converted to
(magnitude = clip(‖flow‖/10, 0, 1), angle = (atan2+π)/2π, mask), resized to 128×128, and encoded by
`SharedMultiHeadVAE` (six-conv encoder → 256-d latent mean μ and log-variance via the reparameterization
trick; a shared transposed-conv decoder with three heads reconstructing flow_x, flow_y and mask; ≈15 M
parameters). Both the input and the reconstruction targets were resized 256→128 by the same bilinear
transform, so training operated entirely at 128×128. The objective was reconstruction MSE + β_KL·KL
(β_KL = 0.001) with a latent-consistency term between paired augmented views (β_consistency = 1.0); 50
augmentations per fibre (rotation ≤ 90°, shear ≤ 5°, scale ± 10%); Adam (lr 1e-3), batch 32, up to 20 epochs,
early stopping (patience 10). The 256-d latent **mean μ** was used as the fibre embedding at inference.

**Pairwise classifier and training.** Two 256-d embeddings were concatenated (→ 512) and passed through an MLP
512 → 128 → 64 → 1 (ReLU, sigmoid; no dropout). Training used BCE loss, Adam (lr 1e-4), batch 256, negatives
sampled at 4× positives (`negative_fold = 4`), up to 50 epochs with early stopping on validation F1
(patience 10); the selected checkpoint reached validation F1 ≈ 0.944.

**Matching.** For each image pair the cost combined the classifier score `S` with a spatial-signature
similarity `G`, the latter the multiscale (k = 3, 5, 7) kNN-distance-vector similarity (per-scale
`1/(1+Wasserstein)`, geometrically averaged). **Seeds:** the 80 top-scoring pairs (`S + G` gated by
`S > 0.75`) were filtered by triangle geometry over combinations of 4 (side cost < 30 px, angle cost < 0.15),
retried up to 3× if fewer than 3 survived. **Propagation:** a 200-px neighbourhood in each section, each
candidate validated by triangle geometry against its 3 nearest matched neighbours, with per-pair patience 5,
to convergence (< 0.25 % new pairs/step). **Fill:** unmatched fibres were localised by an affine transform of
matched centroids (max distance 150 px, min classifier score 0.5) and accepted under the same geometry
validation; duplicates were removed and all matches re-validated. Full defaults are in `configs/default.yaml`.

**Feature extraction.** Each matched fibre's mask was eroded/dilated to define whole, mem (dilate 4 → erode 4),
cyto1 (erode 4 → erode 12) and cyto2 (erode 12) compartments. Per channel per compartment, nine intensity
statistics (mean, s.d., p10, p25, p50, p75, p90, skewness, kurtosis) were computed; with 15 morphological
descriptors this gave 15 + 36 × 22 = 807 features per QUA fibre (663 for TA, which lacks the WGA/Myh panel).
Unmatched channels were left as NaN and excluded from statistics without imputation.
