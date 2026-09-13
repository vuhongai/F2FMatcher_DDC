# Results §A — F2FMatcher: fibre-to-fibre matching across histological stains

*Draft Results/Methods text for theme 1 ("Development of fibre-to-fibre mapping"). This is the opening
section of the paper; it introduces the method and produces the per-fibre multiplex table on which §B
(evaluation) and §C–§D (biology) depend. Figure 1 is the method overview.*

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

## A.2 A stain-invariant shape embedding and a same-fibre classifier

The design principle that makes cross-stain matching tractable is to **encode fibre shape rather than stain
intensity**. Each section is segmented with a fine-tuned Cellpose model, which yields both per-fibre masks
and **flow fields** (the vector field Cellpose uses to define cell boundaries). A 256×256 crop around every
fibre centroid — its two flow components and mask, converted to magnitude/angle and resized to 128×128 — is
encoded by a variational auto-encoder (**VAE**, `SharedMultiHeadVAE`) into a **256-dimensional latent
vector**. Because the input is the flow-field *geometry* and not the staining channel, the embedding of a
given fibre is essentially the same whether that fibre was stained for laminin or for NADH — the property
that lets fibres be compared across modalities.

A **pairwise classifier** then takes two 256-d embeddings and outputs the probability that they are the same
physical fibre (concatenation → MLP → sigmoid), trained with hard negative sampling (4× negatives) and
reaching a validation F1 of ≈ 0.94. Scoring every cross-section pair gives an N₁×N₂ same-fibre score matrix.

Figure 1b shows a worked example on two whole serial sections stained for laminin (immunofluorescence) and
NADH (brightfield): Cellpose segments the thousands of fibres in each, and F2FMatcher assigns a shared
identity to corresponding fibres across the modality gap and across the local distortion between the two
sections.

![Figure 1b](figures/fig1b_worked_example.png)

## A.3 Geometry-aware matching enforces global consistency

Appearance similarity alone cannot disambiguate thousands of similar fibres, so F2FMatcher combines it with
geometry. For every cross-section pair the matching cost is the product of the **classifier score** and a
**spatial-signature similarity** — each fibre's geometric niche is described by multiscale k-nearest-neighbour
distance vectors (k = 3, 5, 7) compared across scales by Wasserstein distance, so a fibre is only a good match
if both its appearance *and* its local neighbourhood agree. Matching proceeds by (i) selecting high-scoring
seed pairs and keeping only those that are mutually consistent under a **triangle-geometry test** (side
lengths and angles of centroid triangles must agree within tolerance), (ii) **iterative local propagation**
that grows the match set outward from validated seeds, each new pair re-checked against its three nearest
already-matched neighbours, until fewer than 0.25% of fibres are added per step, (iii) an **affine fill** that
transforms unmatched centroids into the partner section and accepts geometry-consistent candidates, and (iv)
**post-hoc validation** that removes duplicates and re-checks every match against its local neighbourhood.
Geometry thereby imposes the global spatial consistency that appearance cannot, while tolerating the local
non-rigid distortion between serial sections.

## A.4 From a stained stack to a per-fibre multiplex table

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
fibre are encoded by a VAE into stain-invariant 256-d embeddings; a pairwise classifier scores same-fibre
probability; a geometry-aware algorithm (cost = classifier × spatial signature; triangle-geometry seeds →
iterative propagation → affine fill → validation) assigns fibre correspondences; per-fibre staining is
quantified in four compartments, giving an 807-feature vector per fibre. (b) Worked example: two serial
sections stained for laminin (immunofluorescence, left) and NADH (brightfield, right); Cellpose segmentation
(coloured overlays) and F2FMatcher correspondence link fibres across the modality gap. (c) Each fibre is
partitioned into whole/mem/cyto1/cyto2 compartments (left); the per-fibre feature vector concatenates 15
morphological features with 9 statistics × 4 compartments for each of 22 channels = 807 features (right).

*Optional additional panel (to add): a UMAP of the VAE latent space coloured by fibre, illustrating that the
embedding groups the same fibre across stains — requires exporting the (transient) VAE embeddings.*

---

## Methods (§A — mapping and feature extraction)

**Imaging and segmentation.** Whole-slide CZI images were exported to PNG and resized to a common pixel
resolution (`f2fmatcher.io.czi_reader`). Fibres were segmented with fine-tuned Cellpose 2 models
(`CellPose2_finetuned`; cellprob threshold 0, flow threshold 0.4), producing per-fibre masks and flow fields;
objects below 100 px² were discarded.

**VAE embedding.** For each fibre a 256×256 crop of (flow_x, flow_y, mask) was converted to (magnitude, angle,
mask), resized to 128×128, and encoded by `SharedMultiHeadVAE` (six-conv encoder → 256-d latent via the
reparameterization trick; a shared transposed-conv decoder with three heads reconstructing flow_x, flow_y and
mask). Training minimised reconstruction MSE + β_KL·KL (β_KL = 0.001) + β_consistency·latent-consistency
between paired augmented views (β_consistency = 1.0); 50 augmentations per fibre (rotation ≤ 90°, shear ≤ 5°,
scale ± 10%); Adam (lr 1e-3), batch 32, early stopping (patience 10).

**Pairwise classifier.** Two 256-d embeddings were concatenated and passed through an MLP
(512→128→64→1, ReLU, dropout 0.5, sigmoid). Training used BCE loss, Adam (lr 1e-4), batch 256, negatives at
4× positives, early stopping on validation F1 (patience 10); the selected checkpoint reached validation
F1 ≈ 0.944.

**Matching.** For each image pair the cost was `classifier_score × spatial_signature`, the spatial signature
being the multiscale (k = 3, 5, 7) kNN-distance-vector similarity (per-scale Wasserstein distance,
geometrically averaged). Seeds: the 80 top-scoring pairs, filtered by triangle geometry over combinations of
4 (side cost < 30, angle cost < 0.15). Iterative propagation used a 200-px neighbourhood in each section, each
candidate validated by triangle geometry against its 3 nearest matched neighbours, to convergence (< 0.25%
new pairs/step). Unmatched fibres were filled by an affine transform of matched centroids (max distance 150
px, min classifier score 0.5) with geometry validation; duplicates were removed and all matches re-validated.
Key defaults are in `configs/default.yaml`.

**Feature extraction.** Each matched fibre's mask was eroded/dilated to define whole, mem (dilate 4 → erode
4), cyto1 (erode 4 → erode 12) and cyto2 (erode 12) compartments. Per channel per compartment, nine intensity
statistics (mean, s.d., p10, p25, p50, p75, p90, skewness, kurtosis) were computed; with 15 morphological
descriptors this gave 15 + 36 × 22 = 807 features per QUA fibre (663 for TA, which lacks the WGA/Myh panel).
Unmatched channels were left as NaN and excluded from statistics without imputation.
