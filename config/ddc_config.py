"""
Minimal working configuration for testing
"""
import os
from pathlib import Path

# Mock configuration for testing
CZI_BASE_DIR_TA = "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_TA/"
CZI_BASE_DIR_QUA = "/media/DATABRUT/DB_DDC/serverGPU/AJ/22-082_QUA/"
CP_MASKS_DIR_TA = "/DATA/fiber_mapping_DDC/AJ/analysis_TA_run1/run/out_CP_masks"
CP_MASKS_DIR_QUA = "/DATA/fiber_mapping_DDC/AJ/analysis_Qua_run1/run/out_CP_masks"
PAIR_DIRS_BASE_TA = "/DATA/fiber_mapping_DDC/AJ/analysis_TA_run1/run/prediction_output"
PAIR_DIRS_BASE_QUA = "/DATA/fiber_mapping_DDC/AJ/analysis_Qua_run1/run/prediction_output"
OUTPUT_DIR = "/DATA/F2FMatcher_DDC/results"

# Sample definitions (mock data)
TA_SAMPLES = [
    "TAG01", "TAG02", "TAG03", "TAG04", "TAG05",
    "TAG21", "TAG22", "TAG23", "TAG24", "TAG25",
    "TAG26", "TAG27", "TAG28", "TAG29",
    "TAG31", "TAG32", "TAG33", "TAG34", "TAG35",
]
QUA_SAMPLES = [
    "QUAG01", "QUAG02", "QUAG03", "QUAG04", "QUAG05",
    "QUAG21", "QUAG22", "QUAG23", "QUAG24", "QUAG25",
    "QUAG26", "QUAG27", "QUAG28", "QUAG29",
    "QUAG31", "QUAG32", "QUAG33", "QUAG34", "QUAG35"
]

# Mock staining configuration
"""
'slide_id': {
    'czi_dir': 'directory_name_in_czi_base_dir',
    'IHF': True/False, 
        IHF=True: multiple stainings corresponding to different channels (0-3), 
        IHF=False: brightfield BGR image (e.g. HE), channel_index set to 0
    'stainings': {
        'channel_index': 'staining_name',
        ...    }
    }
"""
SLIDES = {
    1: {
        "czi_dir": "IHF_Lam-Dys-Col4",
        "IHF": True,
        "stainings": {
            0: "DAPI", 1: "Laminin", 2: "Dystrophin", 3: "Collagen4"
        }
    },
    2: {
        "czi_dir": "IHF_Lam-IgG-CD11b",
        "IHF": True,
        "stainings": {
            2: "IgG", 3: "CD11b"
        }
    },
    3: {
        "czi_dir": "NADH",
        "IHF": False,
        "stainings": {
            0: "NADH"
        }
    },
    6: {
        "czi_dir": "HE_10x",
        "IHF": False,
        "stainings": {
            0: "HE_10x"
        }
    },
    7: {
        "czi_dir": "COX",
        "IHF": False,
        "stainings": {
            0: "COX"
        }
    },
    8: {
        "czi_dir": "IHF_LAMP2-LGALS3-SQSTM1",
        "IHF": True,
        "stainings": {
            1: "LAMP2", 2: "LGALS3", 3: "SQSTM1"
        }
    }, 
}


# Compartments to analyze
"""
4 compartments will be analyzed:
- whole: no dilation or erosion, corresponds to the original mask
- mem: dilation of 4 pixels followed by erosion of 4 pixels, corresponds to the muscle fiber membrane
- cyto1: erosion of 4 pixels followed by erosion of 8 pixels, corresponds to the inner cytoplasm close to the membrane, set to 0 if the resulting mask is empty
- cyto2: erosion of 8 pixels followed by no further erosion, corresponds to the inner cytoplasm close to the nucleus, set to 0 if the resulting mask is empty
"""
COMPARTMENTS = {
    "whole": {"dilation": None, "erosion": None}, 
    "mem": {"dilation": 4, "erosion": 4}, 
    "cyto1": {"erosion_1": 4, "erosion_2": 8},
    "cyto2": {"erosion_1": 8, "erosion_2": None},
}

# Feature statistics to compute
FEATURE_STATISTICS = ["mean", "std", "p10", "p25", "p50", "p75", "p90", "skew", "kurt"]

# group mapping
GROUP_MAP = {
    "WT": ["TAG01", "TAG02", "TAG03", "TAG04", "TAG05",
           "QUAG01", "QUAG02", "QUAG03", "QUAG04", "QUAG05"],
    "mdx": ["TAG21", "TAG22", "TAG23", "TAG24", "TAG25",
            "QUAG21", "QUAG22", "QUAG23", "QUAG24", "QUAG25"],
    "AAV9": ["TAG26", "TAG27", "TAG28", "TAG29",
             "QUAG26", "QUAG27", "QUAG28", "QUAG29"],
    "LICA1": ["TAG31", "TAG32", "TAG33", "TAG34", "TAG35",
              "QUAG31", "QUAG32", "QUAG33", "QUAG34", "QUAG35"],
}

GROUP_OF_SAMPLE = {}
for group, samples in GROUP_MAP.items():
    for s in samples:
        GROUP_OF_SAMPLE[s] = group

def get_group(sample_name):
    return GROUP_OF_SAMPLE.get(sample_name, "unknown")

# Analysis parameters
ANALYSIS_CONFIG = {
    "pca_n_components": 50,
    "tsne_perplexity": 30,
    "umap_n_neighbors": 15,
    "umap_min_dist": 0.1,
    "missingness_threshold": 0.2,
    "random_state": 42,
}

# Mask feature descriptions (simplified)
MASK_FEATURES = [
    "area", "perimeter",
    "bbox_height", "bbox_width", "bbox_aspect_ratio", 
    "eccentricity", "solidity",
    "extent", "major_axis_length", "minor_axis_length", "orientation", 
    "roundness", "edge_pixels", "edge_density", "circularity", "compactness"
]