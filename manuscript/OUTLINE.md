# F2FMatcher — publication outline

*Working draft of the argument structure, figures, and suggested additional analyses.
Numbers below come from the QUA analysis notebooks (`notebooks/Q1–Q7`); TA is available as a
second cohort but not yet analysed. Verified citations are collected at the end.*

---

## 0. Strategy: one paper or two?

The material supports **either** a single method-forward paper **or** a two-paper program. Recommendation:

- **Option A (recommended first submission) — one integrated paper**, methods-forward with a strong DMD
  application: *"F2FMatcher enables multiplex single-fibre phenotyping across serial muscle sections and
  reveals fibre-state–resolved responses to micro-dystrophin gene therapy."* Method + evaluation as Fig 1–2,
  biology as Fig 3–6. Fits a methods/tech journal with biological payoff (e.g. *Nat Commun / Nat Methods /
  Cell Rep Methods / eLife*).
- **Option B — two papers**: (1) a **methods** paper (development + evaluation/benchmark) in a
  methods venue; (2) an **application** paper (dystrophic phenotypes + gene-therapy evaluation) in a
  muscle/gene-therapy venue. Choose this if the benchmark grows large enough to stand alone, or the
  biology needs the TA cohort + more animals to be definitive.

The four themes the PI listed map onto sections as: **(1) development → §A**, **(2) evaluation/benchmark →
§B**, **(3) dystrophic phenotypes → §C**, **(4) gene-therapy evaluation → §D**.

---

## §A. Development of fibre-to-fibre mapping (F2FMatcher)

**Core claim.** Corresponding muscle fibres can be matched across serial sections stained with *different*
markers/modalities, at single-fibre resolution and whole-section scale, by combining a shape-based deep
embedding with a geometry-aware assignment — turning a stack of differently-stained sections into one
**multiplex per-fibre table** (here 807 features × 111,403 QUA fibres across 7 slides).

**Arguments / sub-claims**
- A1. The problem is *instance correspondence across modality gaps*, unlike registration (returns a warp,
  not identities) or same-modality keypoint matching. Motivate why naïve approaches fail (stain appearance
  gap; inter-section tissue distortion/tears; thousands of near-identical objects).
- A2. Design: (i) Cellpose segmentation → flow fields; (ii) **VAE** encodes flow fields (shape, not stain
  intensity) into a 256-d latent → stain-invariant by construction; (iii) **pairwise classifier**
  (same-fibre vs not, val **F1 ≈ 0.944**); (iv) **geometry-aware matching**: cost = classifier × spatial
  signature (multiscale kNN Wasserstein), triangle-geometry seed filtering, iterative local propagation,
  affine fill.
- A3. Design rationale: encoding **shape/flow rather than intensity** is what makes it cross-stain;
  geometry enforces global consistency that appearance alone cannot.
- A4. Scales to whole slides and to a 7-panel serial stack anchored on one reference section.

**Figures**
- **Fig 1a** Schematic pipeline (CZI → Cellpose → VAE → classifier → geometry matching → per-fibre table).
- **Fig 1b** VAE latent space (UMAP of embeddings; reconstruction examples of flow fields).
- **Fig 1c** Worked example: two serial sections (different stains), matched contours with shared IDs
  (the repo's `video/comparison.gif` → static panel).
- **Fig 1d** From stack to table: the 22-channel × 4-compartment feature vector per fibre (schematic +
  the compartment definitions whole/mem/cyto1/cyto2).
- *Supp:* architecture details, training curves, hyper-parameters.

**Status:** method exists & documented (`F2FMatcher/README.md`); figures 1a/1c/1d are assembly of existing
assets; 1b needs a short embedding-plot script.

---

## §B. Evaluation of the method (accuracy, robustness, benchmark)

**Core claim.** F2FMatcher matches fibres **accurately, at high coverage, and more robustly across the
stain gap than registration- or appearance-only baselines**, with each algorithmic component contributing.

**Arguments / sub-claims**
- B1. **Coverage** scales with stain-pair difficulty (measured, QUA): same-segmentation IF (Laminin→IgG/CD11b)
  **83%**, NADH **77%**, HE **67%**, COX **61%**, WGA/Myh **48%**, lysosomal slide-8 **44%** (mapping failures,
  see §D caveat); overall **64%** of anchor fibres matched per pair. Coverage is a function of segmentation
  concordance and section quality, not a ceiling of the matcher.
- B2. **Accuracy vs ground truth.** Needed: a labelled correspondence set. Two complementary GT sources:
  (i) **expert manual annotation** of K fibre pairs on N section-pairs spanning easy→hard stain gaps;
  (ii) **self-supervised consistency** — *cycle/transitivity closure* (all slides map to the slide-1 anchor;
  also map slides pairwise and measure A→B→anchor→A closure error) and the **10X-HE ↔ 5X-HE** same-stain,
  different-magnification pair as a near-GT sanity check.
- B3. **Benchmark against baselines** (all producing a per-fibre assignment; report precision/recall/F1,
  coverage, runtime):
  1. Global **affine/rigid centroid registration + nearest-neighbour**.
  2. **Elastic registration** (e.g. ANTs/bUnwarpJ/SIFT-flow) + NN.
  3. **DINOv2** ROI embeddings + mutual-NN / Hungarian (appearance-only, no geometry).
  4. **General keypoint matchers** (SuperPoint+SuperGlue, LoFTR, DISK, RoMa — via
     `alexstoken/image-matching-models`) → dense correspondence → propagate to cells.
  5. **F2FMatcher ablations**: VAE-only (no spatial signature); classifier-only; no triangle geometry; no
     iterative propagation; full model.
  Hypothesis: appearance-only baselines (DINOv2, keypoint matchers) degrade sharply across the IF↔brightfield
  gap and among thousands of similar fibres; registration baselines fail under local tears/distortion;
  F2FMatcher's geometry + shape-embedding wins, and each ablated component costs accuracy.
- B4. **Robustness**: accuracy vs inter-section distance/distortion; vs fibre density; vs segmentation errors.

**Figures**
- **Fig 2a** Coverage per stain-pair (bar; the measured numbers above).
- **Fig 2b** Accuracy (P/R/F1) vs ground truth, F2FMatcher vs baselines 1–4, split by stain-gap difficulty.
- **Fig 2c** Ablation of the four components.
- **Fig 2d** Cycle-consistency closure error (F2FMatcher vs baselines).
- **Fig 2e** Robustness curves (accuracy vs distortion/density) + runtime/scale.
- *Supp:* DINOv2 / image-matching-models configs; example failure cases of each baseline.

**Status:** B1 is **ready** (numbers computed from `mapping_summary.csv`). B2–B4 are the **main new work
for the methods story** — see "Suggested additional analyses". This directly answers the PI's benchmark ask
(DINOv2, `image-matching-models`).

---

## §C. F2FMatcher reveals dystrophic phenotypes of DMD muscle (mdx vs WT)

**Core claim.** Multiplex single-fibre phenotyping resolves DMD pathology into **discrete fibre states** and
recovers **marker–marker relationships invisible to bulk or single-stain analysis** — because every fibre
carries all markers at once.

**Arguments / sub-claims** (all from `Q1`, `Q5`, `Q6`, `Q7`)
- C1. **Unbiased clustering** (morphology + HE only, no fluorescence) yields **6 reproducible fibre states**
  (deterministic GMM, seed=42); staining read out *post hoc* is therefore discovery, not circular. Name them:
  Cls4/Cls6 healthy (oxidative-I / large-fast), Cls2 necrotic-atrophic, Cls3 regenerating-inflamed, Cls5
  large IIb dystrophin-low, Cls1 neutral large-pale.
- C2. **The core DMD signature, fibre-intrinsic and universal:** sarcolemmal dystrophin ↓ (Cohen d≈−3.5,
  p=0.008, *every* cluster); IgG↑ (necrosis, d≈+5.2); CD11b↑ (macrophages, d≈+3.7); central nucleation ↑;
  fibre-size heterogeneity CV 0.40→0.54.
- C3. **Correlation structure is lost in disease** — only measurable with per-fibre multiplexing:
  NADH~MyHC-I 0.88→0.67; a new negative Dystrophin~area coupling; Laminin~WGA and LAMP2~SQSTM1 couplings
  loosen.
- C4. **Fibre states are distinct, sequential lesions** (Q5/Q6): **regenerated (centronucleated)** fibres are
  a separate **oxidative, autophagy-active** state, **not** the same as **necrotic** fibres (IgG-influx,
  small, round, **macrophage-encircled** IgG~CD11b r≈0.75, lysosome-high). Necrosis hits small-fibre clusters
  (Cls2 73%, Cls3 60%), spares large Cls5/6.
- C5. **Lysosomal/autophagy axis** (Q7): LAMP2/LGALS3/SQSTM1 ↑ in mdx; **SQSTM1/p62 is fibre-intrinsic**
  while **LAMP2/LGALS3 are peri-fibre (macrophage-associated)** — separated by the compartment split. Marks
  the necrotic/regenerating niche.

**Figures**
- **Fig 3a** UMAP by group + by 6 clusters (from `fig_q1_umap`).
- **Fig 3b** Cluster-identity heatmap (marker z per cluster; `fig_q1_identity`).
- **Fig 3c** Global mdx–WT effect sizes (`fig_q1_effects`) + key violins.
- **Fig 3d** Within-cluster mdx–WT (fibre-intrinsic dystrophin loss; `fig_q1_wdelta`).
- **Fig 3e** Lost correlations (`fig_q1_corr`).
- **Fig 4** Fibre-state biology: regeneration niche (`fig_q5_niche`), necrosis (`fig_q6_necrosis`),
  lysosomal axis (`fig_q7_lyso`).

**Status:** **fully analysed & figured** (`Q1`, `Q5`, `Q6`, `Q7`). Strengthen with **TA replication** and
**cluster-robustness** (see suggestions).

---

## §D. F2FMatcher reveals a fibre-state–resolved gene-therapy evaluation (AAV9 vs LICA1)

**Core claim.** Per-fibre multiplexing provides a **quantitative, mechanism-resolved efficacy readout** that
bulk assays cannot: it shows *which fibres* are rescued, *how completely*, and *why some resist* — and it
distinguishes two serotypes at the same dose.

**Arguments / sub-claims** (from `Q2`, `Q3`, `Q4`, `Q6`, `Q7`)
- D1. **LICA1 ≫ AAV9 dystrophin restoration**: fraction of fibres reaching WT-level dystrophin **3% (mdx) →
  28% (AAV9) → 65% (LICA1)**; LICA1 ≥ WT in every cluster; AAV9 partial and *highly variable between animals*.
- D2. **Compositional rescue toward healthy**: healthy fraction (Cls4+6) 12%→31%(AAV9)→38%(LICA1); LICA1
  significantly clears the regenerating cluster Cls3, AAV9 does not.
- D3. **Two determinants of restoration** (Q3/Q4): **(i) laminin/basal-lamina scaffold** — predicts ~50% of
  per-fibre membrane-dystrophin variance (β≈0.75, R²≈0.5), *specific* vs collagen-IV → the permissive anchor;
  **(ii) regeneration/vector-retention** — regenerated (CN+) fibres carry **less** transgene (β≈−0.27,
  independent of laminin), the in-situ signature of **AAV episome dilution during regeneration**. LICA1's
  stronger restoration lowers regeneration burden (68→49%), partly breaking the loop.
- D4. **Laminin is dystrophin-independent** (present in dystrophin-null mdx; not raised by restoring
  dystrophin) → the correct causal statement is *laminin (scaffold) + dystrophin → function*, not
  "dystrophin restores laminin".
- D5. **Resistance is architectural, not delivery failure** (Q3): in Cls5, LICA1 loads dystrophin ≥WT yet the
  fibre stays a remodelled fast-glycolytic state and pathology is not resolved → **molecular ≠ morphological
  rescue**.
- D6. **Safety/efficacy signal**: **AAV9 leaves *more* necrosis than untreated mdx** (60% vs 46%), while
  **LICA1 reduces it (38%)**; residual **galectin-3 lysosomal damage persists despite gene therapy** —
  independently reproducing Jaber et al. 2025 (same group) and motivating lysosome-directed combinations.

**Figures**
- **Fig 5a** Cluster-composition shift + healthy/diseased/resistant composites (`fig_q2_comp`,
  `fig_q2_composite`).
- **Fig 5b** Dystrophin restoration: per-animal, %-reaching-WT (3/28/65), per-cluster (`fig_q2_dystr`) +
  dystrophin UMAP maps (`fig_q2_dumap`).
- **Fig 6a** Two-factor model: laminin scaffold (`fig_q3_laminin`) + regeneration penalty
  (`fig_q4_cn_transduction`).
- **Fig 6b** Laminin causality (`fig_q4_laminin_causality`).
- **Fig 6c** Cls5 paradox — molecular vs morphological rescue (`fig_q3_cls5paradox`).
- **Fig 6d** Necrosis & lysosomal residue (`fig_q6_necrosis`, `fig_q7_lyso`).

**Status:** **fully analysed & figured** (`Q2`, `Q3`, `Q4`, `Q6`, `Q7`). Strengthen with more animals /
TA replication / a durability time-course (see suggestions).

---

## Suggested additional analyses (prioritised)

**High priority (make the story publishable & rebuttal-proof)**
1. **Method ground truth + benchmark (§B2–B4).** Build the GT (manual + cycle-consistency + HE-5X/10X) and
   run baselines: registration (affine/elastic), **DINOv2** ROI embeddings, **`image-matching-models`**
   (LoFTR/RoMa/SuperGlue), and F2FMatcher **ablations**. This is the single biggest gap and directly answers
   the PI's benchmark question. *Deliverable: Fig 2.*
2. **TA cross-muscle replication.** Re-run the entire Q1–Q7 pipeline on the TA cohort (data present,
   `results/TA/`, 19 samples, 663 features). Showing the 6-state structure, the LICA1≫AAV9 result, the
   laminin scaffold and the regeneration penalty **replicate in a second muscle** would greatly strengthen
   every biological claim. *Deliverable: a TA panel in Fig 3–6 or a dedicated supp figure.*
3. **Clustering robustness.** BIC/silhouette vs k; bootstrap/seed stability of the 6 states; robustness to the
   feature preset (ch3 vs ch9 vs ch18); show cluster identity is not an artefact of k=6. *Supp.*
4. **Proper mixed-effects statistics.** Replace/complement per-animal medians with linear mixed models
   (fibre nested in animal; group fixed effect) for the headline claims, to satisfy statistical review.

**Medium priority (mechanistic depth)**
5. **Spatial / niche analysis.** Fibre centroids exist (in `paired_labels`/masks): test whether necrotic,
   regenerating, and low-transgene fibres are **spatially clustered** (focal lesions), and whether
   macrophage-high fibres neighbour necrotic ones. Turns "states" into "microenvironments".
6. **Dose/transduction heterogeneity.** Quantify per-animal AAV9 variability (already visible) and relate
   fibre-level dystrophin to a delivered-dose proxy; model the regeneration↔episome-loss feedback explicitly.
7. **Fibre-type–specific rescue.** Formal Myh7/Myh2/Myh4 typing per fibre (not just intensity) and rescue by
   true type; test the IIb-vulnerability/oxidative-protection axis directly.
8. **"Virtual co-staining" demonstrations.** Showcase relationships only obtainable by matching (e.g.
   dystrophin×laminin×IgG per fibre; NADH×MyHC per fibre) as the unique value proposition of F2FMatcher.

**Lower priority / future**
9. **Durability time-course** (new experiment) to convert the inferred regeneration→episome-loss feedback and
   the "no-reversal" resistance into longitudinal evidence.
10. **Combination-therapy pilot** (µDys + laminin-111/α7β1 or + lysosome protection/trehalose, per Jaber 2025)
    predicted to widen the rescuable pool and clear residual galectin-3 damage.
11. **Human-sample application** (mouse is IIb-rich vs human) to establish translational relevance.

---

## Novelty summary (what is new here)

- A **cross-stain, whole-section, single-fibre instance-matching** method that yields a true multiplex
  per-fibre table (method + benchmark).
- **Laminin/basal lamina quantitatively predicts micro-dystrophin membrane restoration in situ** (β≈0.75,
  specific vs collagen-IV) — a single-fibre quantification of the DGC–laminin scaffold requirement.
- **Regeneration imposes a fibre-intrinsic transduction penalty** (regenerated fibres carry less transgene) —
  an in-situ, single-fibre embodiment of AAV episome dilution, independent of the scaffold (two-factor model).
- **Resistance = molecular≠morphological rescue** in a discrete, remodelled fast-glycolytic fibre population
  that nonetheless receives the transgene.
- **AAV9 leaves more necrosis than untreated mdx; LICA1 reduces it**, and **gene therapy leaves residual
  galectin-3 lysosomal damage** (reproducing Jaber 2025).

---

## Verified references (PubMed; carry DOIs into the manuscript)

1. Jaber A, *et al.* Lysosomal damage is a therapeutic target in Duchenne muscular dystrophy. *Sci Adv* 2025;11(43):eadv6805. doi:10.1126/sciadv.adv6805 (PMID 41124255).
2. Le Hir M, *et al.* AAV genome loss from dystrophic mouse muscles during AAV-U7 snRNA-mediated exon-skipping therapy. *Mol Ther* 2013;21(8):1551-8. doi:10.1038/mt.2013.121 (PMID 23752313).
3. Weinmann J, *et al.* Identification of a myotropic AAV by massively parallel in vivo evaluation of barcoded capsid variants. *Nat Commun* 2020;11:5432. doi:10.1038/s41467-020-19230-w (PMID 33116134).
4. Burke MJ, *et al.* AAV8 and AAV9 myofibre type/size tropism profiling reveals therapeutic effect of microdystrophin in canines. *J Cachexia Sarcopenia Muscle* 2025;16(1):e13681. doi:10.1002/jcsm.13681 (PMID 39790021).
5. Bostick B, *et al.* AAV micro-dystrophin gene therapy alleviates stress-induced cardiac death but not myocardial fibrosis in >21-m-old mdx mice. *J Mol Cell Cardiol* 2012;53(2):217-22. doi:10.1016/j.yjmcc.2012.05.002 (PMID 22587991).
6. Bonato A, *et al.* Cyclin D3 deficiency promotes a slower, more oxidative skeletal muscle phenotype and ameliorates pathophysiology in the mdx mouse. *FASEB J* 2023;37(7):e23025. doi:10.1096/fj.202201769R (PMID 37309599).
7. Gumerson JD, Michele DE. The dystrophin-glycoprotein complex in the prevention of muscle damage. *J Biomed Biotechnol* 2011;2011:210797. doi:10.1155/2011/210797 (PMID 22007139).
8. Spaulding HR, *et al.* Autophagic dysfunction and autophagosome escape in the mdx model of DMD. *Acta Physiol* 2018;222(2):e12944. doi:10.1111/apha.12944 (PMID 28834378).

*Tool baselines to cite in §B:* Cellpose; DINOv2 (Oquab et al. 2023); `alexstoken/image-matching-models`
(SuperGlue/LoFTR/RoMa/DISK) — add formal citations when the benchmark is run.
