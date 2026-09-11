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
from collections import abc

def get_filename(name_sample, list_names):
    names = [n for n in list_names if name_sample.upper() in n.upper()]
    if len(names)==0:
        print(f'No image of sample {name_sample} found.')
        return None
    else:
        name = names[0].split(".")[0]
        return name


def convert_to_float_array(features):
    """
    Convert dict_values or similar iterable to numpy array with float64 type.
    Handles NaN values properly.
    """
    # Convert to list first
    if isinstance(features, abc.ValuesView):
        features_list = list(features)
    else:
        features_list = list(features)
    # Convert to numpy array with float64 type
    return np.array(features_list, dtype=np.float64)

# get all images
dict_all_images = {"TA": {}, "QUA": {}}
dict_img2metadata = {"TA": {}, "QUA": {}}

# slides that actually have CZI data for each muscle (e.g. slide 4 = QUA only)
valid_slides = {
    "TA": [s for s in SLIDES if (CZI_BASE_DIR_TA / SLIDES[s]["czi_dir"]).is_dir()],
    "QUA": [s for s in SLIDES if (CZI_BASE_DIR_QUA / SLIDES[s]["czi_dir"]).is_dir()],
}
n_channels = {m: sum(len(SLIDES[s]["stainings"]) for s in valid_slides[m]) for m in valid_slides}
print(f"Valid slides: TA={valid_slides['TA']} ({n_channels['TA']} channels), "
      f"QUA={valid_slides['QUA']} ({n_channels['QUA']} channels)")

for muscle in ["TA", "QUA"]:
    dir_czi_source = CZI_BASE_DIR_TA if muscle == "TA" else CZI_BASE_DIR_QUA
    dir_CP_MASKS = CP_MASKS_DIR_TA if muscle == "TA" else CP_MASKS_DIR_QUA
    dir_pair_output = PAIR_DIRS_BASE_TA if muscle == "TA" else PAIR_DIRS_BASE_QUA
    list_samples = TA_SAMPLES if muscle == "TA" else QUA_SAMPLES

    for sample in list_samples:
        dict_all_images[muscle][sample] = {}
        for slide in valid_slides[muscle]:
            list_images = [f.split(".czi")[0] for f in os.listdir(dir_czi_source / f'{SLIDES[slide]["czi_dir"]}') \
                            if f.endswith(".czi")]
            img = get_filename(sample, list_images)
            dict_all_images[muscle][sample][slide] = img
            dict_img2metadata[muscle][img] = {"sample": sample, "slide": slide}

# list all pairs
dict_mapping = {
    "TA": {sample: [] for sample in TA_SAMPLES}, 
    "QUA": {sample: [] for sample in QUA_SAMPLES}
}
for muscle in ["TA", "QUA"]:
    dir_pair_output = PAIR_DIRS_BASE_TA if muscle == "TA" else PAIR_DIRS_BASE_QUA
    for pairs in os.listdir(dir_pair_output):
        if ("___vs___" in pairs):
            img1, img2 = pairs.split("___vs___")
            meta1 = dict_img2metadata[muscle].get(img1)
            meta2 = dict_img2metadata[muscle].get(img2)
            if meta1 is not None and meta2 is not None:
                assert meta1["sample"] == meta2["sample"]
                p_paired_labels = dir_pair_output / pairs / "paired_labels.pkl"
                if os.path.exists(p_paired_labels):
                    with open(p_paired_labels, "rb") as f:
                        paired_labels = pickle.load(f)
                else:
                    paired_labels = None
                dict_mapping[muscle][meta1["sample"]].append((meta1["slide"], meta2["slide"], pairs, paired_labels))

dict_features_intensity = {
    "TA": {sample: {} for sample in TA_SAMPLES}, 
    "QUA": {sample: {} for sample in QUA_SAMPLES}
}

for muscle in ["TA", "QUA"]:
    dir_pair_output = PAIR_DIRS_BASE_TA if muscle == "TA" else PAIR_DIRS_BASE_QUA
    list_samples = TA_SAMPLES if muscle == "TA" else QUA_SAMPLES

    # dir save all combined features for each sample
    dir_save_features_muscle = Path(f"/DATA/F2FMatcher_DDC/results/{muscle}/features_combined")
    os.makedirs(dir_save_features_muscle, exist_ok=True)
    
    for sample in list_samples:
        for slide in valid_slides[muscle]:
            dir_save_features = Path(f"/DATA/F2FMatcher_DDC/results/{muscle}/features/slide_{slide}/")
            img = dict_all_images[muscle][sample][slide]
            if img is None:
                print(f"[skip] {muscle} {sample} slide {slide}: no image found")
                continue

            dict_features_intensity[muscle][sample][slide] = {}

            list_channels = SLIDES[slide]["stainings"].keys()

            # concatenate features intensity across channel
            for channel in list_channels:
                path_save_features_dict = dir_save_features / f"{img}_c{channel}.pkl"
                with open(path_save_features_dict, "rb") as f:
                    features_dict = pickle.load(f)
                    for f in features_dict:
                        label_id = f["label_id"]
                        features_intensity = np.concatenate([
                            convert_to_float_array(f["intensity"]["whole"].values()),
                            convert_to_float_array(f["intensity"]["mem"].values()),
                            convert_to_float_array(f["intensity"]["cyto1"].values()),
                            convert_to_float_array(f["intensity"]["cyto2"].values())
                            ])
                        assert features_intensity.shape[0] == 36
                        
                        if dict_features_intensity[muscle][sample][slide].get(label_id) is None:
                            dict_features_intensity[muscle][sample][slide][label_id] = []

                        dict_features_intensity[muscle][sample][slide][label_id].append(features_intensity)

        ### COMBINATION OF ALL FEATURES
        # label_id in slide 1 as key, list of features_mask + features_intensity from all channels in slide 2 as value
        features_sample = {}
        list_labels = list(dict_features_intensity[muscle][sample][1].keys())
        for slide in valid_slides[muscle]:
            if slide==1: # reference slide
                for label_id in list_labels:
                    features_sample[label_id] = dict_features_intensity[muscle][sample][slide][label_id]

            else: # to be matched across slides

                mapping_list = [d for d in dict_mapping[muscle][sample] if d[0]==1 and d[1]==slide]
                if mapping_list:
                    dict_map_slide = mapping_list[0][-1]
                    dict_map_slide = {k:v for k,v in dict_map_slide} if dict_map_slide is not None else {}
                else:
                    dict_map_slide = {}

                slide_features = dict_features_intensity[muscle][sample].get(slide, {})
                for label_id in list_labels:
                    label_id_matched = dict_map_slide.get(label_id)
                    if label_id_matched is not None and label_id_matched in slide_features:
                        features_sample[label_id].extend(slide_features[label_id_matched])
                    else:
                        features_sample[label_id].extend([np.full(36, np.nan)] * len(SLIDES[slide]["stainings"].keys()))

        for label_id in list_labels:
            features_sample[label_id] = np.concatenate(features_sample[label_id])
            assert features_sample[label_id].shape[0] == 36*n_channels[muscle]

        # add mask features (slide 1, channel 1)
        dir_save_features = Path(f"/DATA/F2FMatcher_DDC/results/{muscle}/features/slide_1/")
        img = dict_all_images[muscle][sample][1]

        path_save_features_dict = dir_save_features / f"{img}_c1.pkl"
        with open(path_save_features_dict, "rb") as f:
            features_dict = pickle.load(f)

        for f in features_dict:
            label_id = f["label_id"]
            f_mask = convert_to_float_array(f["mask_features"].values())
            assert f_mask.shape[0] == 15
            features_sample[label_id] = np.concatenate([f_mask, features_sample[label_id]])
            
        for label_id in list_labels:
            assert features_sample[label_id].shape[0] == 15+36*n_channels[muscle]

        with open(dir_save_features_muscle / f"{sample}.pkl", "wb") as f:
            pickle.dump(features_sample, f)    