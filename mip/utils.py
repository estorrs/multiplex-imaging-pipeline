import os
import re

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
