import os
import re

import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import tifffile
from tifffile import TiffFile
from ome_types import from_tiff, from_xml
from skimage.exposure import rescale_intensity


def listfiles(folder, regex=None):
    """Return all files with the given regex in the given folder structure"""
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            if regex is None:
                yield os.path.join(root, filename)
            elif re.findall(regex, os.path.join(root, filename)):
                yield os.path.join(root, filename)


def extract_ome_tiff(fp, channels=None):   
    tif = TiffFile(fp)
    ome = from_xml(tif.ome_metadata)
    im = ome.images[0]
    d = {}
    img_channels = []
    for c, p in zip(im.pixels.channels, tif.pages):
        img_channels.append(c.name)

        if channels is None:
            img = p.asarray()
            d[c.name] = img
        elif c.name in channels:
            img = p.asarray()
            d[c.name] = img

    if channels is not None and len(set(channels).intersection(set(img_channels))) != len(channels):
        raise RuntimeError(f'Not all channels were found in ome tiff: {channels} | {img_channels}')

    return d


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
                   cmap=None, add_alpha=True, save_fp=None, vmin=None, vmax=None):
    regions, vals = zip(*region_to_value.items())
    try:
        float(vals[0])
        is_numeric = True
        cmap = cmap if cmap is not None else 'viridis'
        max_val = max(vals) if vmax is None else vmax
        min_val = min(vals) if vmin is None else vmin
        bins = np.linspace(min_val, max_val, 100)
        val_to_color = {b:c for b, c in zip(np.arange(100), sns.color_palette(cmap, n_colors=100))}
    except:
        is_numeric = False
        cats = sorted(set(vals))
        if cmap is None:
            if len(cats) <= 10:
                cmap = 'tab10'
            else:
                cmap = 'tab20'
        val_to_color = {b:c for b, c in zip(cats, sns.color_palette(cmap))}

    rgb = np.full((labeled_img.shape[0], labeled_img.shape[1], 3), 1., dtype=np.float32)

    if is_numeric:
        vals = np.digitize(vals, bins) - 1

    for region_id, val in zip(regions, vals):
        if region_to_bbox is not None:
            r1, c1, r2, c2 = region_to_bbox[region_id]
            cropped_labeled = labeled_img[r1:r2, c1:c2]
            cropped_rgb = rgb[r1:r2, c1:c2]
            cropped_rgb[cropped_labeled==int(region_id)] = val_to_color[val] # faster on cropped
            rgb[r1:r2, c1:c2] = cropped_rgb
        else:
            rgb[labeled_img==int(region_id)] = val_to_color[val]

    if add_alpha:
        rgb = np.concatenate((rgb, np.ones((rgb.shape[0], rgb.shape[1], 1))), axis=-1)
        rgb[labeled_img==0] = [1., 1., 1., 0.]

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
