# F2FMatcher_DDC

Fiber-to-fiber matching and multi-modal feature analysis pipeline for muscle histology study **22-082**. Matches individual muscle fibers across multiple staining protocols on serial sections, then extracts morphological and intensity features per fiber compartment for downstream analysis.

## Study Design

- **Muscles**: Tibialis Anterior (TA), Quadriceps (QUA)
- **Groups**: WT, mdx, AAV9, LICA1
- **Samples**: 18 per muscle (60 total fibers across groups)
- **Source data**: CZI whole-slide images on `LAB_ROOT/DB_DDC/serverGPU/AJ/22-082_{TA,QUA}/`

## Staining Panels

| Slide | Type | Stainings | Segmentation |
|-------|------|-----------|-------------|
| 1 | IHF, 10X | DAPI, Laminin, Dystrophin, Collagen4 | Laminin (ch 1) |
| 2 | IHF, 10X | IgG, CD11b | Laminin (ch 1) |
| 3 | Brightfield, 10X | NADH (R, G, B) | NADH |
| 6 | Brightfield, 10X | HE (R, G, B) | HE |
| 7 | Brightfield, 10X | COX (R, G, B) | HE |
| 8 | IHF, 10X | LAMP2, LGALS3, SQSTM1 | LAMP2 (ch 1) |

## Fiber Compartments

Each fiber is partitioned into four sub-regions via morphological operations:

- **whole** — full fiber mask
- **mem** — membrane ring (dilate 4 px, erode 4 px)
- **cyto1** — peripheral cytoplasm (erode 4 px, erode 12 px)
- **cyto2** — deep cytoplasm (erode 12 px)

## Pipeline

### 1. Feature Extraction

```bash
python scripts/extract_features.py
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
   - 648 intensity features (36 × 18 channels across all slides)
   - **Total: 663 features per fiber**
4. Unmatched fibers get NaN padding for missing channels
5. Saves as `results/{muscle}/features_combined/{sample}.pkl`

### 3. Visualization

Open `scripts/visualization.ipynb` (kernel: `fibermatcher`):
- Loads combined features, adds sample/group metadata
- Standardizes and runs PCA
- Plots PC1 vs PC2 colored by group (WT, mdx, AAV9, LICA1)

## Project Structure

```
F2FMatcher_DDC/
├── config/
│   └── ddc_config.py          # Paths, samples, slides, compartments, analysis params
├── scripts/
│   ├── extract_features.py    # Step 1: per-image feature extraction
│   ├── combine_features.py    # Step 2: cross-slide fiber matching + concatenation
│   └── visualization.ipynb    # Step 3: PCA visualization
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
