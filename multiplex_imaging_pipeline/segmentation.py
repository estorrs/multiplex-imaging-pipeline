import logging

import numpy as np
from skimage.morphology import label
from skimage.measure import regionprops
from deepcell.applications import Mesmer

from multiplex_imaging_pipeline.utils import extract_ome_tiff, merge_channels, R_CHANNEL_MAPPING


APP = Mesmer()


def remove_labeled_edges(labeled, criteria=['top', 'bottom', 'left', 'right']):
    new = labeled.copy()
    pool = set()
    for c in criteria:
        if c == 'top':
            pool.update(np.unique(labeled[0]).tolist())
        if c == 'bottom':
            pool.update(np.unique(labeled[-1]).tolist())
        if c == 'left':
            pool.update(np.unique(labeled[:, 0]).tolist())
        if c == 'right':
            pool.update(np.unique(labeled[:, -1]).tolist())
    
    props = regionprops(labeled)
    for p in props:
        if p.label in pool:
            r1, c1, r2, c2 = p.bbox
            new[r1:r2, c1:c2][labeled[r1:r2, c1:c2]==p.label] = 0
    return new


def keep_first(labeled1, labeled2, threshold=.1, start_idx=None):
    props1 = regionprops(labeled1)
    props2 = regionprops(labeled2)
    new = labeled1.copy()
    
    idx = new.max() + 1 if start_idx is None else start_idx
    for p in props2:
        r1, c1, r2, c2 = p.bbox
        r1, c1 = max(0, r1 - 1), max(0, c1 - 1)
        r2, c2 = r2 + 1, c2 + 1
        tile1 = new[r1:r2, c1:c2]
        tile2 = labeled2[r1:r2, c1:c2].copy()
        
        mask1 = tile1>0
        mask2 = tile2==p.label

        m = mask1 & mask2
        
        overlap_frac = m.sum() / mask2.sum()
                
        if overlap_frac < threshold:
            tile1[mask2] = idx
            new[r1:r2, c1:c2] = tile1
            idx += 1
    return new


def segment_cells(
        fp, split_size=25000, overlap_padding=200, nuclei_channels=['DAPI'],
        membrane_channels=['Pan-Cytokeratin', 'E-cadherin', 'CD45', 'CD8', 'CD3e', 'Vimentin', 'SMA', 'CD31', 'CD20']):
    channel_to_img = extract_ome_tiff(fp)
    
    channel_to_img = {R_CHANNEL_MAPPING.get(k, k):v for k, v in channel_to_img.items()}
    
    nuclei = merge_channels(channel_to_img, nuclei_channels)
    membrane = merge_channels(
        channel_to_img,
        membrane_channels)
    
    # if bigger than ~25k x 25k then they fail sometimes, so spliting them apart
    if nuclei.shape[-2] > split_size or nuclei.shape[-1] > split_size:
        cell_labels, nuclei_labels = (np.zeros_like(nuclei, dtype=np.int32),
                                      np.zeros_like(nuclei, dtype=np.int32))
        nrows = nuclei.shape[-2] // (split_size - overlap_padding) + 1
        ncols = nuclei.shape[-1] // (split_size - overlap_padding) + 1
        logging.info(f'spliting into nrows: {nrows}, ncols: {ncols}')
        for r in range(nrows):
            r1 = r * (split_size - overlap_padding)
            r2 = r1 + split_size
            for c in range(ncols):
                logging.info(f'{r} {c}')
                c1 = c * (split_size - overlap_padding)
                c2 = c1 + split_size
                if r1 < nuclei.shape[-2] and c1 < nuclei.shape[-1]:
                    logging.info(f'window: {r1}, {r2}, {c1}, {c2}')
                    input_img = np.concatenate(
                        (np.expand_dims(nuclei[r1:r2, c1:c2], axis=-1),
                         np.expand_dims(membrane[r1:r2, c1:c2], axis=-1)), axis=-1)
                    input_img = np.expand_dims(input_img, axis=0)
                    if input_img.shape[1] > 256 and input_img.shape[2] > 256:

                        segmentation_predictions = APP.predict(input_img, compartment='both')

                        c_labels = segmentation_predictions[0, ..., 0]
                        n_labels = segmentation_predictions[0, ..., 1]

                        # remove correct edges
                        criteria = ['top', 'bottom', 'left', 'right']
                        if r == 0: criteria.remove('top')
                        if c == 0: criteria.remove('left')
                        if r == nrows - 1: criteria.remove('bottom')
                        if c == ncols - 1: criteria.remove('right')
                        c_labels = remove_labeled_edges(c_labels, criteria=criteria)
                        n_labels = remove_labeled_edges(n_labels, criteria=criteria)

                        # if cells overlap keep original cell
                        cell_labels[r1:r2, c1:c2] = keep_first(
                            cell_labels[r1:r2, c1:c2], c_labels, threshold=.1,
                            start_idx=cell_labels.max() + 1)
                        nuclei_labels[r1:r2, c1:c2] = keep_first(
                            nuclei_labels[r1:r2, c1:c2], n_labels, threshold=.05,
                            start_idx=nuclei_labels.max() + 1)

        labeled_cells, labeled_nuclei = label(cell_labels), label(nuclei_labels)

    else:
        input_img = np.concatenate((np.expand_dims(nuclei, axis=-1),
                                np.expand_dims(membrane, axis=-1)), axis=-1)
        input_img = np.expand_dims(input_img, axis=0)

        segmentation_predictions = APP.predict(input_img, compartment='both')

        labeled_cells, labeled_nuclei = (segmentation_predictions[0, ..., 0],
                                         segmentation_predictions[0, ..., 1])
    
    return labeled_cells, labeled_nuclei