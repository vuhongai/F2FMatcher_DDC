from f2fmatcher.io import czi_reader
from f2fmatcher.segmentation import cellpose_seg
import os, sys, multiprocessing, pickle
from tqdm import tqdm
import pandas as pd
import numpy as np
from PIL import Image
from scipy.ndimage import binary_dilation, binary_erosion
from skimage.morphology import disk
from matplotlib import pyplot as plt
from scipy import ndimage
from skimage import measure, morphology
import cv2
from scipy import stats
Image.MAX_IMAGE_PIXELS = None

sys.path.append("/DATA/F2FMatcher_DDC")
from config.ddc_config import *

def dilate_mask(mask_i, n_px):
    return binary_dilation(mask_i, structure=disk(n_px)).astype(np.uint8)

def erode_mask(mask_i, n_px):
    return binary_erosion(mask_i, structure=disk(n_px)).astype(np.uint8)

def crop_image_with_bbox(image, bbox_i):
    y_min, x_min, y_max, x_max = bbox_i
    image_i = image[y_min:y_max, x_min:x_max]
    return image_i

def compute_features_channel(mask_i, image_i):
    """
    Calculate intensity statistics for each compartment.
    """
    if mask_i.sum() == 0:
        stats_dict = {
            'mean': np.nan,
            'std': np.nan,
            'p10': np.nan,
            'p25': np.nan,
            'p50': np.nan,
            'p75': np.nan,
            'p90': np.nan,
            'skew': np.nan,
            'kurt': np.nan
        }
    else:
        pixels = image_i[mask_i.astype(bool)]
        stats_dict = {
            'mean': np.mean(pixels),
            'std': np.std(pixels),
            'p10': np.percentile(pixels, 10),
            'p25': np.percentile(pixels, 25),
            'p50': np.percentile(pixels, 50),
            'p75': np.percentile(pixels, 75),
            'p90': np.percentile(pixels, 90),
            'skew': stats.skew(pixels),
            'kurt': stats.kurtosis(pixels)
        }
    return stats_dict

def calculate_mask_features(mask):
    """
    Calculate various morphological features from a binary mask.
    
    Parameters:
    -----------
    mask : numpy.ndarray
        Binary mask (2D array with values 0 or 1)
    
    Returns:
    --------
    dict
        Dictionary containing calculated features
    """
    # Ensure the mask is binary
    mask = mask.astype(bool)
    
    # Calculate basic properties using skimage.measure.regionprops
    props = measure.regionprops(mask.astype(int))[0]
    
    # Initialize result dictionary
    features = {}
    
    # Basic features
    features['area'] = props.area
    features['perimeter'] = props.perimeter
    
    # Bounding box features
    min_row, min_col, max_row, max_col = props.bbox
    features['bbox_height'] = max_row - min_row
    features['bbox_width'] = max_col - min_col
    features['bbox_aspect_ratio'] = (max_col - min_col) / max((max_row - min_row), 1)
    
    # Shape features
    features['eccentricity'] = props.eccentricity
    features['solidity'] = props.solidity
    features['extent'] = props.extent
    features['major_axis_length'] = props.major_axis_length
    features['minor_axis_length'] = props.minor_axis_length
    features['orientation'] = props.orientation
    
    # Additional derived features
    # Roundness: ratio of area to area of a circle with same perimeter
    if props.perimeter > 0:
        features['roundness'] = (4 * np.pi * props.area) / (props.perimeter ** 2)
    else:
        features['roundness'] = 0
    
    # Edge pixels: count of pixels on the boundary
    # Find contours and count edge pixels
    contours, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    edge_pixels = sum(len(contour) for contour in contours)
    features['edge_pixels'] = edge_pixels
    
    # Edge density: ratio of edge pixels to area
    if props.area > 0:
        features['edge_density'] = edge_pixels / props.area
    else:
        features['edge_density'] = 0
    
    # Compactness: ratio of area to cube of perimeter
    if props.perimeter > 0:
        features['compactness'] = props.area / (props.perimeter ** 3)
    else:
        features['compactness'] = 0
    
    return features


dir_save_png = Path("./")
n_process = 32


for muscle in ["TA", "QUA"]:
    if muscle == "TA":
        dir_czi_source = CZI_BASE_DIR_TA
        dir_CP_MASKS = CP_MASKS_DIR_TA
        dir_pair_output = PAIR_DIRS_BASE_TA
    elif muscle == "QUA":
        dir_czi_source = CZI_BASE_DIR_QUA
        dir_CP_MASKS = CP_MASKS_DIR_QUA
        dir_pair_output = PAIR_DIRS_BASE_QUA

    for slide in SLIDES.keys():
        list_images = [f.split(".czi")[0] for f in os.listdir(dir_czi_source / f'{SLIDES[slide]["czi_dir"]}') \
                        if f.endswith(".czi")]
        for img in list_images:
            for channel in SLIDES[slide]["stainings"].keys():
    
                IHF = SLIDES[slide]["IHF"]
                staining = SLIDES[slide]["stainings"][channel]
                if IHF:
                    param_img = ("fluorescence", SLIDES[slide]["scanning_objective"], 1.0)
                else:
                    param_img = ("brightfield", SLIDES[slide]["scanning_objective"], 1.0)

                dir_czi = dir_czi_source / f'{SLIDES[slide]["czi_dir"]}' / f'{img}.czi'
                dir_png = dir_save_png / f'{img}.png'

                CP_model_name = SLIDES[slide]["segmentation_model"]
                CP_model_path = Path("/DATA/F2FMatcher/models/CellPose2_finetuned")

                dir_save_features = Path(f"/DATA/F2FMatcher_DDC/results/{muscle}/features/slide_{slide}/")
                os.makedirs(dir_save_features, exist_ok=True)

                path_save_features_dict = dir_save_features / f"{img}_c{channel}.pkl"
                # only process if the feature dictionary does not exist yet
                if path_save_features_dict.exists():
                    print(f"Feature dictionary for {img}_c{channel} already exists, skipping...")
                    continue
                else:
                    # import czi, resize and export as png
                    czi_reader.import_resize_export_czi(
                        czi_path=dir_czi,
                        IHF=IHF,
                        channel_index=channel,
                        dir_save_png=dir_save_png,
                        param_img=param_img,
                    )

                    # get cellpose segmentation masks and properties
                    masks, props, _ = cellpose_seg.get_CP_masks(
                        img_path=dir_png, 
                        CP_model_name=CP_model_name, 
                        CP_model_path=CP_model_path,
                        savedir=dir_CP_MASKS, 
                        channels = [0, 0],
                    )

                    # filter small fibers
                    region_labels = {
                        int(r.label):[r.bbox, r.area] for r in props if r.area >= threshold_fiber_area
                    }

                    # analysis per ROI
                    def calc_intensity_staining(label_id):
                        bbox_i = region_labels[label_id][0]

                        # update bbox after dilation, used it for cropping the image later
                        n_px_dilation = -min(list_erosion)
                        y_min, x_min, y_max, x_max = bbox_i
                        y_min = max(0, y_min-n_px_dilation)
                        y_max = min(y_full, y_max+n_px_dilation)
                        x_min = max(0, x_min-n_px_dilation)
                        x_max = min(x_full, x_max+n_px_dilation)

                        # crop around mask (prevent from calculation in large matrix --> speed up calculation)
                        bbox_i = [y_min, x_min, y_max, x_max]
                        image_i = crop_image_with_bbox(image, bbox_i)
                        mask_i = (masks == label_id).astype(np.uint8)
                        mask_i = mask_i[y_min:y_max, x_min:x_max]
                        
                        # quantification
                        cyto2 = erode_mask(mask_i, list_erosion[2])
                        cyto1 = erode_mask(mask_i, list_erosion[1])-cyto2
                        mem = dilate_mask(mask_i, n_px_dilation)-(cyto1+cyto2)

                        # mask features
                        features = calculate_mask_features(mask_i)

                        # intensity features
                        stats_whole = compute_features_channel(mask_i, image_i)
                        stats_cyto1 = compute_features_channel(cyto1, image_i)
                        stats_cyto2 = compute_features_channel(cyto2, image_i)
                        stats_mem = compute_features_channel(mem, image_i)

                        return {
                            "label_id": label_id,
                            "mask_features": features,
                            "intensity": {
                                "whole": stats_whole,
                                "cyto1": stats_cyto1,
                                "cyto2": stats_cyto2,
                                "mem": stats_mem
                            }
                        }

                    # open resized PNG image
                    image = np.array(Image.open(dir_png))
                    if not IHF:
                        # select the channel
                        image = image[:,:,channel]
                        # convert to grayscale with black background if brightfield
                        image = 255 - image

                    assert masks.shape == image.shape
                    y_full, x_full = masks.shape

                    args_list = list(region_labels.keys())
                    with multiprocessing.Pool(processes=n_process) as pool:
                        feature_dicts = list(
                            tqdm(
                                pool.imap_unordered(calc_intensity_staining, args_list),
                                total=len(args_list),
                                desc=f"Processing image {img} channel {channel}"
                            )
                        )

                    with open(path_save_features_dict, "wb") as f:
                        pickle.dump(feature_dicts, f)

                    os.system(f"rm {dir_png}")