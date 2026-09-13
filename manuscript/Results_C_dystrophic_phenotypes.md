# Results §C — F2FMatcher reveals the single-fibre pathology of dystrophic muscle

*Draft Results text for theme 3 ("F2FMatcher reveals dystrophic phenotypes of DMD muscle"). Figures are
numbered from 3 onward, reserving Figures 1–2 for the method and its evaluation (§A, §B). All statistics are
computed at the animal level (per-animal medians; Mann–Whitney U and Cohen's d; n = 5 WT, 5 mdx, 4 AAV9, 5
LICA1) unless stated; effect sizes are quoted because with n = 5 vs 5 the smallest attainable U-test p is
0.008 (complete separation). Only observed — i.e. actually matched — marker values enter the statistics;
unmatched channels are never imputed.*

By matching every fibre across the seven serial-section panels, F2FMatcher assembled a single multiplex table
of 807 features (15 morphological + 36 intensity statistics × 22 staining channels, in four sub-cellular
compartments) for each of 111,403 quadriceps fibres from 19 animals (WT, mdx, and two gene-therapy groups).
This per-fibre multiplexing is what enables the analyses below, which no single stain or bulk assay can
provide.

## C.1 An unbiased six-state fibre landscape

To describe fibre heterogeneity without circularity, we clustered fibres using **only morphology and
haematoxylin–eosin brightfield texture** (123 features; no immunofluorescence), then read out every other
marker *post hoc*. A z-score → PCA(30) → UMAP → Gaussian-mixture pipeline (k = 6, fixed seed) partitioned the
fibres into six reproducible states (Fig 3a,b). The same six states, in identical proportions, were recovered
by an independent re-analysis, and their group composition already separated the genotypes on the unbiased
map: WT fibres occupied the upper "healthy" territory, mdx fibres spread into lower regions with little WT
support, and treated fibres shifted back toward the WT core (Fig 3a).

![Figure 3a,b](figures/figC3_landscape.png)

Reading out the full marker panel per state defined their biology (Fig 3c): **Cls4** (WT-dominated, 62%) —
small, oxidative type-I fibres with the highest sarcolemmal dystrophin, high laminin and MyHC-I; **Cls6** (WT
55%) — large fast fibres, dystrophin-high, lowest inflammation; **Cls3** (mdx 38%) — small, regenerating
fibres with the highest central nucleation and high oxidative, immune and lysosomal signals; **Cls2** (mdx
40%) — the smallest, most necrotic/fibrotic fibres (highest IgG, high CD11b/collagen-IV/LAMP2); **Cls5**
(virtually no WT, 1%) — large, fast-glycolytic (highest MyHC-IIb) fibres with the *lowest* dystrophin and
laminin; and **Cls1** — very large, pale, MyHC-IIb fibres distributed evenly across groups.

![Figure 3c](figures/figC3_identity.png)

Thus two healthy poles (Cls4, Cls6) oppose two diseased poles (Cls2 necrotic, Cls3 regenerating), with Cls5 a
distinct large dystrophin-low IIb population — a decomposition obtained without using any disease marker.

## C.2 The mdx disease signature is fibre-intrinsic and universal

Comparing mdx with WT across the marker panel recovered the textbook DMD signature with large, fully
separated effects (Fig 4a): loss of sarcolemmal **dystrophin** (49.8 → 30.6 a.u.; Cohen's d ≈ −3.5, p =
0.008), a rise in intracellular **IgG** marking sarcolemmal leak/necrosis (d ≈ +5.2), increased **CD11b⁺**
myeloid infiltration (d ≈ +3.7), higher central nucleation and roundness (regeneration), and elevated
lysosomal LAMP2. Fibre-size heterogeneity — a histological hallmark of dystrophy — increased sharply: the
per-animal coefficient of variation of fibre area rose from 0.40 (WT) to 0.54 (mdx), with both atrophic
(<2000 px²: 1%→7%) and hypertrophic (>10,000 px²: 7%→19%) fibres expanding (Fig 4b).

![Figure 4a,b](figures/figC4_effects.png)

![Figure 4b](figures/figC4_hetero.png)

Because fibres within a cluster are morphologically matched, comparing WT and mdx *inside each cluster*
separates molecular from morphological change. Dystrophin was reduced in **every** cluster, including the
healthy Cls4/Cls6, and IgG/CD11b rose in every cluster (Fig 4c) — establishing that dystrophin loss and
membrane leak are **fibre-intrinsic**, not merely a consequence of the shift in cluster membership. The
diseased clusters (Cls2/Cls3) additionally gained central nucleation and lysosomal LAMP2/LGALS3.

![Figure 4c](figures/figC4_within.png)

## C.3 Disease dissolves the physiological coupling between markers

Per-fibre multiplexing further revealed that *relationships* between markers — not only their levels — are
lost in disease (Fig 4d). The tight coupling between oxidative metabolism and slow fibre type
(NADH ~ MyHC-I) fell from r = 0.88 (WT) to 0.67 (mdx); a negative dystrophin–area coupling appeared in mdx
(larger fibres losing more dystrophin); and the basal-lamina/membrane (Laminin ~ WGA) and lysosome–autophagy
(LAMP2 ~ SQSTM1) couplings loosened. The dystrophin ~ laminin correlation, by contrast, remained high in mdx
(0.83 → 0.91): residual/revertant dystrophin still tracked the intact basal lamina. Loss of coordinated
marker scaling is a signature obtainable only when all markers are measured in the same fibre.

![Figure 4d](figures/figC4_corr.png)

## C.4 Regeneration, necrosis and lysosomal stress are distinct fibre states

We next defined two fibre states directly from single-fibre signals: **regenerated (centronucleated, CN⁺)**
fibres, scored from DAPI intensity variance in the deep fibre core (a central nucleus injects bright DAPI
into an otherwise-dark core; the peripheral-nucleus compartment was deliberately avoided), and **necrotic
(IgG⁺)** fibres, defined by sarcolemmal IgG influx above the WT baseline.

Centronucleation rose from 10% of fibres in WT to 68% in mdx and marked the disease clusters (Cls5 87%, Cls3
75%; healthy Cls4/Cls6 8–16%). Within mdx, CN⁺ (regenerating) fibres were distinctly **oxidative** (NADH
+0.49, COX +0.54 SD) and **autophagy-active** (LAMP2/LGALS3/SQSTM1 up), but were **no more necrotic or
macrophage-encircled than CN⁻ fibres** (42% vs 40% necrotic; 34% vs 34% macrophage-high; Fig 5a). Necrotic
fibres, in turn, were small, round, **macrophage-encircled** (peri-fibre CD11b +0.86 SD; IgG ~ CD11b Spearman
r ≈ 0.75), lysosome-loaded, and concentrated in the small-fibre clusters (Cls2 73%, Cls3 60%, Cls4 61%),
sparing the large clusters (Cls5 20%, Cls6 13%) (Fig 5b). Regeneration and necrosis are therefore **distinct,
sequential fibre states** — a necrotic fibre is cleared and a *new* fibre regenerates — that co-elevate in
disease but rarely coincide in the same fibre.

![Figure 5a](figures/figC5_regen.png)

![Figure 5b](figures/figC5_necrosis.png)

Finally, the autophagy–lysosome axis was perturbed in dystrophic fibres: LAMP2, galectin-3 (LGALS3) and
p62/SQSTM1 all increased in mdx (Fig 5c), consistent with the known block in autophagic flux and p62
accumulation in dystrophic muscle [Spaulding 2018]. Compartmental read-out separated their origin: **p62/
SQSTM1 was intra-fibre** (rising toward the core), whereas **LAMP2 and galectin-3 were concentrated
peri-fibre** — a large macrophage contribution — and all three tracked the necrotic/inflammatory/regenerating
niche (Fig 5c). Galectin-3, a sensor of lysosomal-membrane permeabilization, thus reports an active lysosomal
damage compartment in dystrophic muscle [Jaber 2025].

![Figure 5c](figures/figC5_lyso.png)

Together, §C shows that F2FMatcher converts a stack of differently-stained sections into a single-fibre atlas
that recovers the complete DMD phenotype, demonstrates its fibre-intrinsic nature, exposes the loss of
physiological marker coupling, and resolves the degeneration–regeneration cycle into distinct fibre states —
the baseline against which therapy is evaluated in §D.

---

### Figure legends (§C)

**Figure 3. An unbiased six-state fibre landscape of quadriceps muscle.**
(a) UMAP of 111,403 fibres embedded from morphology + HE features only (no immunofluorescence), coloured by
group. (b) The same map coloured by the six Gaussian-mixture states (Cls1–Cls6; deterministic, seed = 42).
(c) Cluster identity: median level of each marker per state, z-scored across states (red = high, blue = low),
naming the states biologically. Markers span all seven panels and four compartments.

**Figure 4. The mdx single-fibre disease signature.**
(a) Standardized effect size (Cohen's d, per-animal medians) of mdx vs WT for each marker; red = higher in
mdx. (b) Fibre-size distribution and per-animal coefficient of variation of fibre area by group. (c)
Within-cluster mdx − WT standardized median difference per marker × cluster: dystrophin falls and IgG/CD11b
rise in every state (fibre-intrinsic). (d) Marker–marker Spearman couplings by group and NADH ~ MyHC-I
hexbin (WT vs mdx), showing loss of physiological coupling in disease.

**Figure 5. Regeneration, necrosis and lysosomal stress are distinct fibre states.**
(a) Centronucleated (regenerated) vs non-regenerated fibres within mdx: profile across metabolism, matrix,
immune and lysosomal markers (left), NADH by group in CN⁺ vs CN⁻ (middle), and near-identical overlap with
necrosis/macrophage (right). (b) Necrotic (IgG⁺) fibres within mdx: profile (left), per-animal necrotic
fraction by group (middle), and necrotic prevalence by cluster with the IgG ~ CD11b coupling (right). (c)
Lysosomal markers (LAMP2/LGALS3/SQSTM1; AAV9 excluded for slide-8 mapping): per-animal levels (left),
compartmental localization showing p62 intra-fibre vs LAMP2/galectin-3 peri-fibre (middle), and correlations
with the necrotic/inflammatory/regenerating niche (right).

*Note: for slide-8 markers (LAMP2/LGALS3/SQSTM1), AAV9 is excluded throughout because its slide-8 mapping
coverage is only ~16%; brightfield NADH/COX intensity can be inflated by staining uptake in leaky fibres, so
their metabolic direction in necrotic fibres is read cautiously.*
