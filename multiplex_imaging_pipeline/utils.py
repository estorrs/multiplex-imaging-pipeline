import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import tifffile
from tifffile import TiffFile
from ome_types import from_tiff, from_xml
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from skimage.exposure import rescale_intensity
from skimage.segmentation import find_boundaries



CHANNEL_MAPPING = {
    'Pan-Cytokeratin': ['Pan-Cytokeratin', 'Pan-CK', 'Pan-CK (D)', 'PanCK (D)'],
    'E-cadherin': ['E-cadherin'],
    'CD45': ['CD45 (D)', 'CD45', 'CD45-(D)', 'CD45RO'],
    'CD8': ['CD8', 'CD8a'],
    'DAPI': ['DAPI'],
    'CD4': ['CD4', 'CD4 (D)'],
    'CD3e': ['CD3e', 'CD3'],
    'Vimentin': ['Vimentin-(D)', 'Vimentin', 'Vimentin (D)'],
    'SMA': ['SMA-(D)', 'SMA', 'SMA (D)', 'a-SMA (D)'],
    'CD31': ['CD31', 'CD31-(D)'],
    'CD20': ['CD20-(D)', 'CD20', 'CD20 (D)', 'CD20-Akoya'],
    'CD68': ['CD68', 'CD68 (D)', 'CD68-(D)'],
    'CD163': ['CD163'],
    'FOXP3': ['FoxP3', 'FOXP3', 'Foxp3'],
    'cKit': ['cKIT (D)', 'cKit (D)', 'cKit', 'cKIT', 'cKIT-(D)'],
    'MGP': ['MGP', 'MGP (D)', 'MGP-(D)'],
    'CD36': ['CD36-(D)', 'CD36'],
    'PR': ['PR', 'PR-(D)'],
    'ER': ['ER', 'ER-(D)'],
    'P21': ['P21', 'P21-(D)', 'P21-(D)', 'P21 (D)'],
    'P16': ['P16-(D)'],
    'CK5': ['Keratin 5', 'KRT5'],
    'TFF1': ['TFF1', 'TFF1-(D)'],
    'beta-integrin': ['beta-integrin', 'beta3-integrin'],
    'CK14': ['CK14', 'Keratin 14', 'CK14 (D)'],
    'CK17': ['CK17', 'CK17 (D)', 'CK17-(D)', 'Keratin 17'],
    'CK19': ['CK19', 'CK19 (D)', 'CK19-(D)'],
    'CD11b': ['CD11b', 'CD11b (D)', 'CD11b-(D)'],
    'GATA3': ['GATA3', 'GATA3-(D)', 'GATA3 (D)'],
    'PLAT/tPA': ['PLAT/tPA', 'PLAT/tPA (D)'],
    'COX6c': ['COX6c (D)'],
    'Her2': ['Her2', 'Her2 (D)'],
    'Bap1': ['Bap1 (D)'],
    'GLUT1': ['Glut1', 'GLUT1', 'GLUT1-(D)'],
    'CD11c': ['CD11c'],
    'HLA-DR': ['HLA-DR'],
    'Ki67': ['Ki67', 'KI67'],
    'Podoplanin': ['Podoplanin', 'PDPN', 'Podoplanin (D)'],
    'CTLA4': ['CTLA4', 'CTLA4 (D)', 'CTLA4-(D)'],
    'SLC39A6': ['SLC39A6'],
    'BCA1': ['BCA1'],
    'BCAL': ['BCAL'],
    'TUBB3': ['TUBB3', 'TUBB3 (D)'],
    'PTPRZ1': ['PTPRZ1', 'PTPRZ1 (D)'],
    'HIF1a': ['HIF1a', 'HIF1a (D)'],
    'PAI1': ['PAI1', 'PAI1-(D)'],
    'GFAP': ['GFAP', 'GFAP (D)'],
    'VEGFA': ['VEGFA', 'VEGFA (D)'],
    'IBA1': ['IBA1', 'IBA1 (D)'],
    'OLIG2': ['OLIG2', 'OLIG2 (D)'],
    'FN1': ['FN1', 'FN1 (D)'],
}
R_CHANNEL_MAPPING = {v:k for k, vs in CHANNEL_MAPPING.items() for v in vs}


def listfiles(folder, regex=None):
    """Return all files with the given regex in the given folder structure"""
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            if regex is None:
                yield os.path.join(root, filename)
            elif re.findall(regex, os.path.join(root, filename)):
                yield os.path.join(root, filename)


def extract_ome_tiff(fp, channels=None, as_dict=True):   
    tif = TiffFile(fp)
    ome = from_xml(tif.ome_metadata)
    im = ome.images[0]
    d = {}
    img_channels, img_channels, imgs = [], [], []
    for c, p in zip(im.pixels.channels, tif.pages):
        if channels is None or c.name in channels:
            img = p.asarray()
            d[c.name] = img
            imgs.append(img)
            img_channels.append(c.name)
    if channels is not None and len(set(channels).intersection(set(img_channels))) != len(channels):
        raise RuntimeError(f'Not all channels were found in ome tiff: {channels} | {img_channels}')
    
    if as_dict:
        return d
    return img_channels, np.asarray(imgs)



def get_ome_tiff_channels(fp):   
    tif = TiffFile(fp)
    ome = from_xml(tif.ome_metadata)
    im = ome.images[0]
    return [c.name for c in im.pixels.channels]


def create_circular_mask(h, w, center=None, radius=None):
    """
    https://stackoverflow.com/questions/44865023/how-can-i-create-a-circular-mask-for-a-numpy-array
    """

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask


def display_slide_from_adata(a, x='centroid_col', y='centroid_row_inverted', color=None, scale=1000, size=10):
    fig, ax = plt.subplots(
        figsize=(int(max(a.obs[x]) / scale), int(max(a.obs[y]) / scale)))
    sc.pl.scatter(a, x=x, y=y, color=color, ax=ax, size=size)


def display_region(labeled_img, region_to_value, region_to_bbox,
                   cmap=None, add_alpha=True, save_fp=None, vmin=None, vmax=None, cats=None):
    regions, vals = zip(*region_to_value.items())
    vals_original = np.asarray(vals)


    try:
        float(vals[0])
        is_numeric = True
        cmap = cmap if cmap is not None else 'viridis'
        pool = [v for v in vals if not pd.isnull(v)] + [0]
        max_val = np.max(pool) if vmax is None else vmax
        min_val = np.min(pool) if vmin is None else vmin
        bins = np.linspace(min_val, max_val, 100)
        val_to_color = {b:c for b, c in zip(np.arange(100), sns.color_palette(cmap, n_colors=100))}
    except:
        is_numeric = False
        cats = sorted(set(vals)) if cats is None else cats
        if cmap is None:
            if len(cats) <= 10:
                cmap = 'tab10'
            else:
                cmap = 'tab20'
        if isinstance(cmap, str):
            val_to_color = {b:c for b, c in zip(cats, sns.color_palette(cmap))}
        else:
            val_to_color = cmap
    for val in [np.nan, -1, -1., None]:
        val_to_color[val] = (.6, .6, .6)

    rgb = np.full((labeled_img.shape[0], labeled_img.shape[1], 3), 1., dtype=np.float32)

    if is_numeric:
        vals = np.digitize(vals, bins) - 1
        vals[pd.isnull(vals_original)] = -1

    region_mask = np.zeros((labeled_img.shape[0], labeled_img.shape[1]), dtype=bool)
    for region_id, val in zip(regions, vals):
        if region_to_bbox is not None:
            r1, c1, r2, c2 = region_to_bbox[region_id]
            cropped_labeled = labeled_img[r1:r2, c1:c2]
            cropped_rgb = rgb[r1:r2, c1:c2]
            cropped_rgb[cropped_labeled==int(region_id)] = val_to_color[val] # faster on cropped
            rgb[r1:r2, c1:c2] = cropped_rgb
            region_mask[r1:r2, c1:c2] = cropped_labeled==int(region_id)
        else:
            rgb[labeled_img==int(region_id)] = val_to_color[val]
            region_mask[labeled_img==int(region_id)] = True

    # if add_alpha:
    #     rgb = np.concatenate((rgb, np.ones((rgb.shape[0], rgb.shape[1], 1))), axis=-1)
    #     # rgb[labeled_img==0] = [1., 1., 1., 0.]
    #     rgb[~region_mask] = [1., 1., 1., 0.]

    return rgb


def make_pseudo(channel_to_img, cmap=None, contrast_pct=20.):
    cmap = sns.color_palette('tab10') if cmap is None else cmap

    new = np.zeros_like(next(iter(channel_to_img.values())))
    img_stack = []
    for i, (channel, img) in enumerate(channel_to_img.items()):
        color = cmap[i] if not isinstance(cmap, dict) else cmap[channel]
        new = img.copy().astype(np.float32)
        new -= new.min()
        new /= new.max()

        try:
            vmax = np.percentile(new[new>0], (contrast_pct)) if np.count_nonzero(new) else 1.
            new = rescale_intensity(new, in_range=(0., vmax))
        except IndexError:
            pass

        new = np.repeat(np.expand_dims(new, -1), 3, axis=-1)
        new *= color
        img_stack.append(new)
    stack = np.mean(np.asarray(img_stack), axis=0)
    stack -= stack.min()
    stack /= stack.max()
    return stack


def merge_channels(channel_to_img, channels, thresh=.01, contrast_pct=95.):
    img = None
    for c in channels:
        if c not in channel_to_img:
            print(f'warning: {c} not in image')
            pass
        else:
            X = channel_to_img[c].copy().astype(np.float32)
            X -= X.min()
            X /= X.max()
            vmax = np.percentile(X[X>0], (contrast_pct)) if np.count_nonzero(X) else 1.
            X = rescale_intensity(X, in_range=(0., vmax))

            X = np.expand_dims(X, 0)

            if img is None:
                img = X
            else:
                img = np.concatenate((img, X), axis=0)

    img = np.mean(img, axis=0)
    return img


def color_regions(labeled, label_to_annotation, label_to_bbox, mapping):
    rgba = np.zeros((labeled.shape[0], labeled.shape[1], 4))
    for label, annot in label_to_annotation.items():
        color = mapping[annot]
        r1, r2, c1, c2 = label_to_bbox[label]
        
        rgba[r1:r2, c1:c2][labeled[r1:r2, c1:c2]==label] = color
        
        bounds = find_boundaries(labeled[r1:r2, c1:c2]==label, mode='thick')
        
        rgba[r1:r2, c1:c2][bounds] = [.85, .85, .85, 1.]
    
    return rgba


def save_annotated_rgba(labeled_img, adata, filepath, figsize=(10, 10), dpi=300):
    label_to_annotation = {int(l):annot
                           for l, annot in zip(adata.obs.index, adata.obs['cell_type'])}
    label_to_bbox = {int(l):(r1, r2, c1, c2) for l, r1, r2, c1, c2 in zip(
        adata.obs.index, adata.obs['bbox-r1'], adata.obs['bbox-r2'],
        adata.obs['bbox-c1'], adata.obs['bbox-c2'])}
    mapping = {annot:c + (1.,)
            for annot, c in zip(
                sorted(set(label_to_annotation.values())), sns.color_palette('tab20'))}

    rgba = color_regions(labeled_img, label_to_annotation, label_to_bbox, mapping)
    mosaic = [
        ['a', 'a', 'a', 'a', 'b'],
        ['a', 'a', 'a', 'a', 'b'],
        ['a', 'a', 'a', 'a', 'b'],
        ['a', 'a', 'a', 'a', 'b'],
    ]
    fig, axs = plt.subplot_mosaic(mosaic, layout='constrained', figsize=figsize)
    for name, ax in axs.items():
        ax.axis('off')

    ax = axs['a']
    ax.imshow(rgba)

    ax = axs['b']
    custom = [Line2D([0], [0], marker='s', color=(.85, .85, .85), label=label,
                    markerfacecolor=color, markersize=10, linestyle='')
            for label, color in mapping.items()]
    ax.legend(handles=custom, loc='center')

    plt.savefig(filepath, dpi=dpi)