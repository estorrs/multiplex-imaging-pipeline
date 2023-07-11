import logging
import math
from collections import Counter

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import tifffile
import torchvision.transforms.functional as TF
import torch
from scipy.ndimage import binary_fill_holes
from skimage.filters import gaussian
from skimage.measure import regionprops_table, regionprops
from skimage.morphology import label, binary_erosion
from skimage.segmentation import expand_labels
from skimage.exposure import rescale_intensity
from einops import rearrange

import multiplex_imaging_pipeline.utils as utils


def fast_expand_labels(labeled, dist, props=None):
    if props is None:
        props = regionprops(labeled)
    new = np.zeros_like(labeled)
    for prop in props:
        r1, c1, r2, c2 = prop.bbox
        r1, c1 = max(r1 - dist * 2, 0), max(c1 - dist * 2, 0)
        r2, c2 = r2 + dist * 2, c2 + dist * 2
        new[r1:r2, c1:c2][expand_labels(labeled[r1:r2, c1:c2], dist)==prop.label] = prop.label
    return new

def fast_shrink_labels(labeled, dist, props=None):
    if props is None:
        props = regionprops(labeled)
    new = np.zeros_like(labeled)
    for prop in props:
        r1, c1, r2, c2 = prop.bbox
        r1, c1 = max(r1 - dist * 2, 0), max(c1 - dist * 2, 0)
        r2, c2 = r2 + dist * 2, c2 + dist * 2
        
        tile = labeled[r1:r2, c1:c2]==prop.label
        for i in range(dist):
            tile = binary_erosion(tile)
        new[r1:r2, c1:c2][labeled[r1:r2, c1:c2]==prop.label] = 0
        new[r1:r2, c1:c2][tile] = prop.label
    return new

def compute_polsby_popper(area, perimeter):
    try:
        return 4 * math.pi * area / perimeter ** 2
    except ZeroDivisionError:
        return math.nan

def get_regionprops_df(labeled_img,
                       feats=('label', 'bbox', 'area', 'perimeter', 'centroid')):
    regions = regionprops_table(labeled_img, properties=feats)
    region_df = pd.DataFrame(regions)
    region_df = region_df.set_index('label')
    region_df.columns = ['r1', 'c1', 'r2', 'c2', 'area', 'perimeter', 'row', 'col']
    region_df['compactness'] = compute_polsby_popper(region_df['area'], region_df['perimeter'])
    return region_df

def get_morphology_masks(labeled_img, myoepi_dist=(20, 40), boundary_dist=150):
    props = regionprops(labeled_img)

    mask = labeled_img > 0
    
    ductal = fast_shrink_labels(labeled_img, myoepi_dist[0], props=props)
    ductal_mask = ductal > 0

    # myoepithelial
    outer_labeled = fast_expand_labels(labeled_img, myoepi_dist[1], props=props)
    myoepi_labeled = np.zeros_like(labeled_img)
    outer_mask = outer_labeled>0
    myoepi_mask = outer_mask ^ ductal_mask
    myoepi_labeled[myoepi_mask] = outer_labeled[myoepi_mask]

    # boundary
    outer_labeled = fast_expand_labels(labeled_img, boundary_dist, props=props)
    boundary_labeled = np.zeros_like(labeled_img)
    outer_mask = outer_labeled>0
    boundary_mask = outer_mask ^ mask
    boundary_labeled[boundary_mask] = outer_labeled[boundary_mask]

    labeled_dict = {
        'region': labeled_img.astype(np.int16),
        'ductal': ductal.astype(np.int16),
        'myoepi': myoepi_labeled.astype(np.int16),
        'boundary': boundary_labeled.astype(np.int16),
        'expanded': outer_labeled.astype(np.int16),
    }
    mask_dict = {
        'region': mask,
        'ductal': ductal_mask,
        'myoepi': myoepi_mask,
        'boundary': boundary_mask,
        'expanded': outer_mask
    }

    return labeled_dict, mask_dict

def calculate_annotation_fractions(adata, labeled, key):
    pool = sorted(set(adata.obs[key]))
    regions = []
    for row, col in zip(adata.obs['row'], adata.obs['col']):
        row, col = int(row), int(col)
        regions.append(labeled[row, col])
    regions = np.asarray(regions)
    
    data = []
    for region in range(1, labeled.max() + 1, 1):
        f = adata.obs[regions==region]
        counts = Counter(f[key])
        if f.shape[0]:
            data.append([counts.get(cat, 0) / f.shape[0] for cat in pool])
        else:
            data.append([np.nan] * len(pool))
    df = pd.DataFrame(data=data, columns=pool, index=list(range(1, labeled.max() + 1, 1)))
    return df
        

def calculate_marker_intensities(channel_to_img, labeled):
    data, labels = [], []
    channels = sorted(set(channel_to_img.keys()))
    stacked = np.concatenate([np.expand_dims(channel_to_img[c], -1) for c in channels], axis=-1)
    props = regionprops(labeled, intensity_image=stacked)
    for prop in props:
        labels.append(prop.label)
        data.append(prop.intensity_mean.tolist())
   
    df = pd.DataFrame(data=data, columns=channels, index=labels)
    return df

def calculate_marker_fractions(channel_to_img, adata, labeled):
    data, labels = [], []
    channels = [c.replace('_fraction', '') for c in adata.var.index.to_list()]
    channel_to_thresh = {c:thresh for c, thresh in zip(channels, adata.uns['thresholds'])
                if thresh > 0 and thresh not in [254., 255., 65534., 65535.]}
    logging.info(f'Using the following threhsolds to calculate positive fraction: {channel_to_thresh}')
    channels = sorted(set(channel_to_thresh.keys()))
    stacked = np.concatenate([np.expand_dims(channel_to_img[c] >= channel_to_thresh[c], -1)
                              for c in channels], axis=-1)
    props = regionprops(labeled, intensity_image=stacked)
    for p in props:
        labels.append(p.label)
        data.append(
            p.image_intensity.sum(axis=(0, 1)) / (p.image_intensity.shape[0] * p.image_intensity.shape[1]))

    df = pd.DataFrame(data=data, columns=channels, index=labels)
    return df

def generate_mask(mask, sigma=2., min_area=10000):
    mask = gaussian(mask, sigma=sigma)
    mask = binary_fill_holes(mask)
    
    labeled = label(mask)
    props = regionprops(labeled)
    
    mask = labeled > 0
    for prop in props:
        if prop.area < min_area:
            r1, c1, r2, c2 = prop.bbox
            r1, c1 = max(r1 - 10, 0), max(c1 - 10, 0)
            r2, c2, = r2 + 10, c2 + 10
            mask[r1:r2, c1:c2][labeled[r1:r2, c1:c2]==prop.label] = 0
 
    return mask

def generate_region_mask(channel_to_img, channel_to_thresholds, min_area=10000, sigma=1.):
    channels = sorted(set(channel_to_thresholds.keys()))
    mask = np.zeros_like(next(iter(channel_to_img.values())), dtype=bool)
    for c in channels:
        if c in channel_to_thresholds:
            val = channel_to_thresholds[c]
            m = channel_to_img[c] >= val
            mask |= m
        else:
            logging.info(f'Warning: channel {c} does not have a threshold.')

    mask = generate_mask(mask, sigma=sigma, min_area=min_area)

    return mask

def get_region_features(
        adata_fp, ome_fp, mask_fp=None,
        mask_markers=None, min_area=10000, sigma=1., default_threshold=2.):
    a = sc.read_h5ad(adata_fp)
    channel_to_img = utils.extract_ome_tiff(ome_fp)
    channel_to_img = {utils.R_CHANNEL_MAPPING.get(k, k):v for k, v in channel_to_img.items()}
    # channel_to_img_scaled = {c:img / img.std() for c, img in channel_to_img.items()}

    if mask_fp is not None:
        mask = tifffile.imread(mask_fp)
    else:
        logging.info(f'no mask detected, generating mask from {mask_markers}')
        channels = [c.replace('_fraction', '') for c in a.var.index.to_list()]
        channel_to_thresh = {c:thresh for c, thresh in zip(channels, a.uns['thresholds'])
                            if not pd.isnull(thresh)
                            if thresh > 0
                            if thresh not in [254., 255., 65534., 65535.]}
        mask_channel_to_img = {c:img for c, img in channel_to_img.items()
                               if c in mask_markers}
        channel_to_thresholds = {c:v for c, v in channel_to_thresh.items()
                                 if c in mask_markers}
        # mask_channel_to_img = {c:img for c, img in channel_to_img_scaled.items()
        #                        if c in mask_markers}
        # channel_to_thresholds = {c:default_threshold for c in mask_channel_to_img.keys()}
        mask = generate_region_mask(mask_channel_to_img, channel_to_thresholds,
                                    min_area=min_area, sigma=sigma)

    labeled = label(mask)
    labeled_dict, _ = get_morphology_masks(labeled)

    type_to_table = {}
    for k, img in labeled_dict.items():
        logging.info(f'generating region features for {k}')
        df_meta = get_regionprops_df(img)

        logging.info(f'generating cell fractions')
        df_fracs = calculate_annotation_fractions(a, img, 'cell_type')
        df_fracs.columns = [f'cell_type_fraction_{c}' for c in df_fracs.columns]

        logging.info(f'generating marker intensities')
        df_marker_intensities = calculate_marker_intensities(channel_to_img, img)
        df_marker_intensities.columns = [f'marker_intensity_{c}'
                                         for c in df_marker_intensities.columns]
        
        logging.info(f'generating marker intensities scaled')
        df_marker_intensities_scaled = calculate_marker_intensities(channel_to_img_scaled, img)
        df_marker_intensities_scaled.columns = [f'marker_intensity_scaled_{c}'
                                         for c in df_marker_intensities_scaled.columns]

        combined = pd.merge(
            df_meta, df_fracs, left_index=True, right_index=True, how='left')
        combined = pd.merge(
            combined, df_marker_intensities, left_index=True, right_index=True, how='left')
        combined = pd.merge(
            combined, df_marker_intensities_scaled, left_index=True, right_index=True, how='left')
        
        logging.info(f'generating marker fractions')
        if 'thresholds' in a.uns and a.uns['thresholds'] is not None:
            df_marker_fracs = calculate_marker_fractions(channel_to_img, a, img)
            df_marker_fracs.columns = [f'marker_fraction_{c}'
                                    for c in df_marker_fracs.columns]
            combined = pd.merge(
                combined, df_marker_fracs, left_index=True, right_index=True, how='left')
        
        type_to_table[k] = combined

    combined = None
    for name, df in type_to_table.items():
        df.columns = [f'{name}_{c}' for c in df.columns]
        if combined is None:
            combined = df.copy()
        else:
            combined = pd.concat((combined, df), axis=1)

    return combined, labeled_dict
