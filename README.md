# F2FMatcher_DDC

F2FMatcher_DDC is a workflow for analyzing fiber data from histological slides with different staining protocols. This project implements the complete pipeline for:

1. **Mapping CSV Generation**: Creates a mapping CSV with image_name (original image name without extension), slide, sample_id, muscle_type, and group columns
2. **Feature Extraction**: Extracts meaningful features from masks across different staining types and slides
3. **Data Compilation**: Compiles features into a single CSV where each fiber contains all information across slides/stainings
4. **Clustering Analysis**: Performs PCA, tSNE, and UMAP clustering to group fibers across samples

## Project Structure

```
/DATA/F2FMatcher_DDC/
├── config/
│   └── ddc_config.py          # Configuration file with paths and sample definitions
├── scripts/
├── results/
│   └── TA/
│   └── QUA/  
└── README.md                  # This file
```