import logging
import json
from collections import Counter

import anndata
import pandas as pd
import numpy as np
import tifffile
from einops import rearrange
from skimage.measure import regionprops

import multiplex_imaging_pipeline.utils as utils

DEFAULT_PIXEL_FRACTION = .05
DEFAULT_GATING_STRATEGY = [
    {
        'name': 'Epithelial',
        'strategy': [
            {'channel': 'Pan-Cytokeratin', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'}
        ],
    },
    {
        'name': 'Epithelial',
        'strategy': [
            {'channel': 'E-cadherin', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'}
        ],
    },
    {
        'name': 'Treg',
        'strategy': [
            {'channel': 'CD3e', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
            {'channel': 'CD4', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
            {'channel': 'FOXP3', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'}
        ],
    },
    {
        'name': 'CD4 T cell',
        'strategy': [
            {'channel': 'CD3e', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
            {'channel': 'CD4', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'CD8 T cell',
        'strategy': [
            {'channel': 'CD3e', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
            {'channel': 'CD8', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'}
        ],
    },
    {
        'name': 'B cell',
        'strategy': [
            {'channel': 'CD20', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
            {'channel': 'Pan-Cytokeratin', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'neg'},
            {'channel': 'E-cadherin', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'neg'},
        ],
    },
    {
        'name': 'Macrophage - M1',
        'strategy': [
            {'channel': 'CD68', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'}
        ],
    },
    {
        'name': 'Macrophage - M2',
        'strategy': [
            {'channel': 'CD68', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
            {'channel': 'CD163', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'}
        ],
    },
    {
        'name': 'Endothelial',
        'strategy': [
            {'channel': 'CD31', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'Immune',
        'strategy': [
            {'channel': 'CD45', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'Immune',
        'strategy': [
            {'channel': 'CD3e', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'Immune',
        'strategy': [
            {'channel': 'CD4', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'Immune',
        'strategy': [
            {'channel': 'CD68', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'Immune',
        'strategy': [
            {'channel': 'CD163', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'Stroma',
        'strategy': [
            {'channel': 'SMA', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
    {
        'name': 'Stroma',
        'strategy': [
            {'channel': 'Podoplanin', 'value': DEFAULT_PIXEL_FRACTION, 'direction': 'pos'},
        ],
    },
]


def auto_calculate_threshold(img, frac=.003, size=10000, n=5):
    std = img.std()
    x = img / std
    rows, cols = np.random.choice(x.shape[0], size=size), np.random.choice(x.shape[1], size=size)
    vals = x[rows, cols]
    _, counts = np.unique(vals, return_counts=True)
    peak = counts[1:n + 1].mean()
    thresh = peak * frac
    thresh *= std
    return thresh


def generate_feature_table(ome_fp, seg_fp, thresholds=None):
    logging.info(f'extracting {ome_fp}')
    channels, imgs = utils.extract_ome_tiff(ome_fp, as_dict=False)

    logging.info(f'extracting {seg_fp}')
    seg = tifffile.imread(seg_fp)
    if thresholds is not None:
        thresholds = np.asarray([thresholds.get(c, 0.) for c in channels])
        logging.info(f'thresholds detected: {list(zip(channels, thresholds))}')

    else:
        logging.info('thresholds not detected. auto-calculating.')
        thresholds = np.asarray([auto_calculate_threshold(img) for img in imgs])
        logging.info(f'thresholds: {list(zip(channels, thresholds))}')

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
        masks_tile = masks[..., r1:r2, c1:c2]

        cell_mask = seg_tile==label

        row = [label, prop['centroid'][0], prop['centroid'][1], r1, c1, r2, c2, area]
        for j in range(imgs_tile.shape[0]):
            img = imgs_tile[j]

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

    # scaled
    cols = [c for c in df.columns if '_intensity' in c]
    for c in cols:
        df[f'{c}_scaled'] = df[c] / df[c].std()
    
    order = ['label', 'row', 'col', 'bbox-r1', 'bbox-c1', 'bbox-r2', 'bbox-c2', 'area']
    rest = [c for c in df.columns if c not in order]
    order += sorted(rest)
    df = df[order]

    return df, thresholds

def gate_cells(df, key=None, gating_strategy=None, default_threshold=None):
    if key is None:
        key = 'fraction' if any(['_fraction' in c for c in df.columns]) else 'intensity_scaled'
    logging.info(f'cell typing key is: {key}')
    markers = [c.replace(f'_{key}', '') for c in df.columns if f'_{key}' in c]
    logging.info('gating cells with the following markers: ' + str(list(markers)))
    cell_types = np.asarray(['Unlabeled'] * df.shape[0], dtype=object)
    gating_strategy = DEFAULT_GATING_STRATEGY if gating_strategy is None else gating_strategy
    for d in gating_strategy:
        is_valid = cell_types=='Unlabeled'
        mask = np.ones_like(cell_types, dtype=np.bool)
        for strategy in d['strategy']:
            channel, val = strategy['channel'], strategy['value']
            if default_threshold is None:
                if key == 'fraction':
                    pass
                elif key == 'intensity_scaled':
                    val = 1.
            else:
                val = default_threshold

            x = df[f'{channel}_{key}'].to_list()[0] if f'{channel}_{key}' in df.columns else 0
            if pd.isnull(x) or f'{channel}_{key}' not in df.columns:
                if len(d) == 1:
                    m = np.ones_like(mask, dtype=np.bool)
                else:
                    m = np.zeros_like(mask, dtype=np.bool)
            elif strategy['direction'] == 'pos': 
                m = df[f'{channel}_{key}'] >= val
            else:
                m = df[f'{channel}_{key}'] < val
            mask &= m

        if ~mask.sum() == 0: # if all true then pass bc that means marker isn't present or something is wrong with gate
            pass
        else:
            mask &= is_valid
            cell_types[mask] = d['name']
    df['cell_type'] = cell_types

    counts = Counter(cell_types).most_common()
    logging.info(f'cells gated: {counts}')

    return df


def get_spatial_features(labeled_fp, ome_fp, thresholds=None, output_prefix='output', key=None, gating_strategy=None):
    df, thresholds = generate_feature_table(ome_fp, labeled_fp, thresholds=thresholds)

    if key is None:
        key = 'fraction' if any(['_fraction' in c for c in df.columns]) else 'intensity_scaled'

    df = gate_cells(df, key=key, gating_strategy=gating_strategy)

    val_cols = [c for c in df.columns if f'_{key}' in c]
    meta = df[[c for c in df.columns if c not in val_cols]]
    meta = meta.set_index('label')
    a = anndata.AnnData(X=df[val_cols].values.astype(np.float32), obs=meta)
    a.var.index = [c.replace(f'_{key}', '') for c in val_cols]

    a.uns['thresholds'] = thresholds
    a.uns['gating_strategy'] = json.dumps(gating_strategy)

    a.obsm['spatial'] = meta[['col', 'row']].values

    return df, a
    