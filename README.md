# F2FMatcher_DDC

Fiber-to-fiber matching and multi-modal feature analysis pipeline for muscle histology study **22-082**. Matches individual muscle fibers across multiple staining protocols on serial sections, then extracts morphological and intensity features per fiber compartment for downstream analysis.

## Study Design

- **Muscles**: Tibialis Anterior (TA), Quadriceps (QUA)
- **Groups**: WT, mdx, AAV9, LICA1
- **Samples**: 18 per muscle (60 total fibers across groups)
- **Source data**: CZI whole-slide images on `LAB_ROOT/DB_DDC/serverGPU/AJ/22-082_{TA,QUA}/`

## Staining Panels

| Slide | Type | Stainings | Segmentation | Muscles |
|-------|------|-----------|-------------|---------|
| 1 | IHF, 10X | DAPI, Laminin, Dystrophin, Collagen4 | Laminin (ch 1) | TA, QUA |
| 2 | IHF, 10X | IgG, CD11b | Laminin (ch 1) | TA, QUA |
| 3 | Brightfield, 10X | NADH (R, G, B) | NADH | TA, QUA |
| 4 | IHF, 10X | WGA, Myh7, Myh2, Myh4 | WGA (ch 0) | QUA only |
| 6 | Brightfield, 10X | HE (R, G, B) | HE | TA, QUA |
| 7 | Brightfield, 10X | COX (R, G, B) | HE | TA, QUA |
| 8 | IHF, 10X | LAMP2, LGALS3, SQSTM1 | LAMP2 (ch 1) | TA, QUA |

## Fiber Compartments

Each fiber is partitioned into four sub-regions via morphological operations:

- **whole** — full fiber mask
- **mem** — membrane ring (dilate 4 px, erode 4 px)
- **cyto1** — peripheral cytoplasm (erode 4 px, erode 12 px)
- **cyto2** — deep cytoplasm (erode 12 px)

## Pipeline

### 0. Fiber Mapping (F2FMatcher)

```bash
python scripts/run_mapping.py            # all pairs in config/QUA_Myh4_mapping_files.csv
python scripts/evaluate_mapping.py       # per-pair/per-panel coverage -> results/QUA/mapping_summary.csv
```

Producer/consumer pipeline: `n_seg_workers` threads precompute PNG export + CellPose
segmentation (persistent masks/flows) while the main thread runs the matching stage
of the previous pair (npz crops + VAE embeddings are transient, deleted after each pair).
Outputs under `results/QUA/{images_segmentation,out_CP_masks,prediction_output}`.

### 1. Feature Extraction

```bash
python scripts/extract_features.py                       # all muscles, all slides
python scripts/extract_features.py --muscle QUA --slide 4
```

For each CZI image, each staining channel:
1. Resizes CZI to PNG via `f2fmatcher.io.czi_reader`
2. Segments fibers with CellPose (finetuned models in `/DATA/F2FMatcher/models/CellPose2_finetuned`)
3. Filters fibers by area threshold (≥100 px²)
4. Computes per-fiber features:
   - **15 morphological features**: area, perimeter, bbox, eccentricity, solidity, roundness, etc.
   - **36 intensity features**: 9 statistics (mean, std, p10-p90, skew, kurt) × 4 compartments
5. Saves as pickle: `results/{muscle}/features/slide_{N}/{image}_c{channel}.pkl`

Resumes automatically — skips images with existing `.pkl` files. Runs with 32 parallel workers.

### 2. Feature Combination

```bash
python scripts/combine_features.py
```

For each sample:
1. Loads per-slide feature pickles
2. Uses F2FMatcher paired labels (`prediction_output/*___vs__*/paired_labels.pkl`) to match fibers across slides
3. Concatenates features into a single vector per fiber:
    - 15 mask features (from slide 1, ch 1)
    - 36 intensity features × N channels (N = channels of the slides available for that muscle)
    - **Total: TA = 15 + 36×18 = 663, QUA = 15 + 36×22 = 807 features per fiber**
4. Unmatched fibers get NaN padding for missing channels
5. Saves (replaces) `results/{muscle}/features_combined/{sample}.pkl`

### 3. Visualization

Open `scripts/visualization.ipynb` (kernel: `fibermatcher`):
- Loads combined features, adds sample/group metadata
- Standardizes and runs PCA
- Plots PC1 vs PC2 colored by group (WT, mdx, AAV9, LICA1)

### 4. Clustering (unbiased) + Q1/Q2 analysis notebooks

**Clustering setup (precomputed, do not re-run):** unbiased on morphology (15 features) +
HE_10x RGB (108) = 123 features, global z-score, PCA(30), UMAP(nn=15, min_dist=0),
GMM(k=6) on the PCA scores (full covariance, n_init=1, max_iter=200, random_state=42).
Cache: `results/QUA/clustering_cache/ch3_zscore_pca30_mindist0.0_nn15.npz` (111,403 fibers,
filter = detected on slide 6). Plots: `scripts/plot_umap.py`.

**Analysis bundle:** `scripts/build_analysis_bundle.py` derives
`results/QUA/analysis_bundle.npz` from the precomputed cache + `features_combined`
(raw 807-dim features, UMAP, groups, samples, GMM k=6 labels) — no re-clustering.

**Notebooks** (kernel: `fibermatcher`, run top-to-bottom; figures in `visualizations/Q1|Q2/`):
- `notebooks/Q1_mdx_vs_WT.ipynb` — the mdx disease signature: cluster identities,
  mdx vs WT per-cluster marker effects (Cliff's delta), dystrophin loss, correlations lost in mdx.
- `notebooks/Q2_treatment_effect_AAV9_vs_LICA1.ipynb` — treatment effect: cluster
  composition shift (mdx -> AAV9 -> LICA1), dystrophin restoration (LICA1 > AAV9 in all 6 clusters),
  correlations restored, slide-8 lysosomal panel (AAV9 excluded: mapping failed for QUAG27/28/29),
  C5 deep dive (treatment-resistant cluster).

**Clusters (k=6):** C1 large quiet mixed; C2 small fibrotic/inflammatory; C3 small active + immune
(mdx-enriched); C4 small healthy fast (WT-enriched); C5 large Myh4/BM-defect (mdx-specific,
treatment-resistant); C6 large quiet WT-like.

## Project Structure

```
F2FMatcher_DDC/
├── config/
│   └── ddc_config.py          # Paths, samples, slides, compartments, analysis params
├── scripts/
│   ├── run_mapping.py         # Step 0: F2FMatcher mapping (segmentation pipelined with matching)
│   ├── evaluate_mapping.py    # Step 0: mapping quality (per-pair/per-panel coverage)
│   ├── extract_features.py    # Step 1: per-image feature extraction
│   ├── combine_features.py    # Step 2: cross-slide fiber matching + concatenation
│   ├── visualization.ipynb    # Step 3: PCA visualization
│   ├── cluster_fibers.py      # Step 4: UMAP cache builder (filter->impute->zscore->PCA->UMAP)
│   ├── plot_umap.py           # Step 4: UMAP plots from cache (group / cluster / feature)
│   ├── plot_stain_correlation.py  # Step 4: density scatter of 2 stainings per group
│   └── build_analysis_bundle.py   # Step 5: analysis_bundle.npz for the Q1/Q2 notebooks
├── notebooks/
│   ├── Q1_mdx_vs_WT.ipynb
│   └── Q2_treatment_effect_AAV9_vs_LICA1.ipynb
├── results/
│   ├── TA/
│   │   ├── features/slide_{N}/   # Per-image .pkl files
│   │   └── features_combined/    # Per-sample .pkl files (663 features)
│   └── QUA/
│       ├── features/slide_{N}/
│       └── features_combined/
├── visualizations/
└── README.md
```

## Dependencies

- `f2fmatcher` (internal library for CZI reading, CellPose segmentation, fiber matching)
- CellPose 2 (finetuned models)
- scikit-image, scipy, opencv-python, numpy, pandas
- matplotlib, seaborn, scikit-learn
- tqdm, Pillow

## Configuration

All paths, sample lists, slide definitions, and analysis parameters are in `config/ddc_config.py`. Edit this file to:
- Change source/output directories
- Add/remove samples or slides
- Adjust fiber area threshold (`threshold_fiber_area`)
- Modify compartment erosion values (`list_erosion`)
- Tune analysis parameters (`ANALYSIS_CONFIG`)
