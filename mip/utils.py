import os
import re

import numpy as np
import tifffile
from tifffile import TiffFile
from ome_types import from_tiff, from_xml


def listfiles(folder, regex=None):
    """Return all files with the given regex in the given folder structure"""
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            if regex is None:
                yield os.path.join(root, filename)
            elif re.findall(regex, os.path.join(root, filename)):
                yield os.path.join(root, filename)


def extract_ome_tiff(fp):   
    tif = TiffFile(fp)
    ome = from_xml(tif.ome_metadata)
    im = ome.images[0]
    d = {}
    for c, p in zip(im.pixels.channels, tif.pages):
        img = p.asarray()
        d[c.name] = img
    return d


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
