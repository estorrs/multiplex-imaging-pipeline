import logging

import anndata
import pandas as pd
import numpy as np
import tifffile
from einops import rearrange
from skimage.measure import regionprops

import multiplex_imaging_pipeline.utils as utils

def generate_feature_table(ome_fp, seg_fp, thresholds=None):
    logging.info(f'extracting {ome_fp}')
    channels, imgs = utils.extract_ome_tiff(ome_fp, as_dict=False)
    assert len(channels) == len(thresholds)
    
    logging.info(f'extracting {seg_fp}')
    seg = tifffile.imread(seg_fp)
    if thresholds is not None:
        thresholds = np.asarray(thresholds) # make sure numpy
        logging.info(f'thresholds detected: {thresholds}')
        masks = rearrange(
            rearrange(imgs, 'c h w -> h w c') > thresholds,
            'h w c -> c h w')
    
    props = regionprops(seg)
    logging.info(f'num cells: {len(props)}')
    
    data = []
    for i, prop in enumerate(props):
        label = i + 1
        row = []
        r1, c1, r2, c2 = prop['bbox']

        area = prop['area']
        seg_tile = seg[r1:r2, c1:c2]
        imgs_tile = imgs[..., r1:r2, c1:c2]
        if thresholds is not None:
            masks_tile = masks[..., r1:r2, c1:c2]

        cell_mask = seg_tile==label

        row = [label, prop['centroid'][0], prop['centroid'][1], r1, c1, r2, c2, area]
        for j in range(imgs_tile.shape[0]):
            img = imgs_tile[j]

            if thresholds is not None:
                mask = masks_tile[j]
                counts = (cell_mask & mask).sum()
                row.append(counts / area)
            
            intensity = img[cell_mask].mean()
            row.append(intensity)

        data.append(row)

    cols = ['label', 'row', 'col', 'bbox-r1', 'bbox-c1', 'bbox-r2', 'bbox-c2', 'area']
    for c in channels:
        converted = utils.R_CHANNEL_MAPPING.get(c, c)
        roots = ['fraction', 'intensity'] if thresholds is not None else ['intensity']
        for identifier in roots:
            cols.append(f'{converted}_{identifier}')
    df = pd.DataFrame(data=data, columns=cols)
    return df
