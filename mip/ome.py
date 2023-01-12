import logging
import re
import os

import numpy as np
import tifffile
from ome_types import from_tiff, from_xml, to_xml, model
from ome_types.model.simple_types import UnitsLength

from mip.utils import listfiles

# an effort to have consistent marker names in all level 2 ome.tiffs
# getting rid of it for now because it messes with some of the qc
d = {
}
CHANNEL_MAP = {v:k for k, vs in d.items() for v in vs}

# calculated from andrew shinkle measurement on HT206 test slide
PHENOCYCLER_PIXELS_PER_MICRON = 1.9604911906033102


def parse_codex_channel_name_from_raw(c):
    pieces = c.split('_')
    c = pieces[-1]
    if 'DAPI1' in c or 'DAPI-01' in c or 'DAPI-1' in c:
        return 'DAPI'
    elif 'DAPI' in c:
        return None
    elif 'blank' in c.lower():
        return None
    elif 'empty' in c.lower():
        return None
    elif 'Experiment' in pieces[0]:
        return None
    return c

def identity(c):
    return c


def generate_ome_from_tifs(fps, output_fp, platform='codex', bbox=None):
    """
    Generate an HTAN compatible ome tiff from a list of filepaths, where each filepath is a tiff representing a different channel.

    Assumes files have .tif extension and channel name is parsable from filename root.
    """
    names = [fp.split('/')[-1].replace('.tif', '') for fp in fps]

    if platform == 'codex':
        mapping = parse_codex_channel_name_from_raw
    elif platform == 'raw':
        mapping = identity
    else:
        raise RuntimeError(f'The platform {platform} is not a valid platform')

    keep_idxs, keep = zip(*[(i, n) for i, n in enumerate(names)
                            if mapping(n) is not None])
    new = [CHANNEL_MAP.get(mapping(n), mapping(n))
           for n in keep]
    name_to_identifier = {k:n for k, n in zip(keep, new)}
    
    x, y = None, None
    logging.info(f'bbox detected, cropping to {bbox}')
    with tifffile.TiffWriter(output_fp, ome=True, bigtiff=True) as out_tif:
        for i, fp in enumerate(fps):
            if i in keep_idxs:
                img = tifffile.imread(fp)
                x, y = img.shape[1], img.shape[0]
                if bbox is not None:
                    r1, r2, c1, c2 = bbox
                    y = r2 - r1
                    x = c2 - c1
                    img = img[r1:r2, c1:c2]
                out_tif.write(img.astype(np.uint16))
        o = model.OME()
        o.images.append(model.Image(id='Image:0', pixels=model.Pixels(dimension_order='XYCZT',
              size_c=len(keep_idxs),
              size_t=1,
              size_x=x,
              size_y=y,
              size_z=1,
              type='float',
              big_endian=False,
              channels=[model.Channel(id=f'Channel:{i}', name=f'{name_to_identifier[c]}') for i, c in enumerate(keep)],
              physical_size_x=float(x),
              physical_size_y=float(y),)))
        
        im = o.images[0]
        im.pixels.physical_size_x_unit = 'µm'
        im.pixels.physical_size_y_unit = 'µm'
        for i in range(len(im.pixels.channels)):
            im.pixels.planes.append(model.Plane(the_c=i, the_t=0, the_z=0))
        im.pixels.tiff_data_blocks.append(model.TiffData(plane_count=len(im.pixels.channels)))
        xml_str = to_xml(o)
        out_tif.overwrite_description(xml_str.encode())


def generate_ome_from_qptiff(qptiff_fp, output_fp, bbox=None):
    """
    Generate an HTAN compatible ome tiff from qptiff output by codex phenocycler
    """
    tf = tifffile.TiffFile(qptiff_fp)

    s = tf.series[0] # full res tiffs are in the first series
    n_channels = s.get_shape()[0]
    logging.info(f'image has {n_channels} total biomarkers')
    logging.info(f'bbox detected, cropping to {bbox}')
    x, y = None, None
    with tifffile.TiffWriter(output_fp, ome=True, bigtiff=True) as out_tif:
        biomarkers = []
        for i, p in enumerate(s.pages):
            img = p.asarray()
            x, y = img.shape[1], img.shape[0]
            d = tifffile.xml2dict(p.description)['PerkinElmer-QPI-ImageDescription']
            biomarker = d['Biomarker']
            biomarkers.append(biomarker)
            logging.info(f'writing {biomarker}')
            if bbox is not None:
                r1, r2, c1, c2 = bbox
                y = r2 - r1
                x = c2 - c1
                img = img[r1:r2, c1:c2]

            out_tif.write(img.astype(np.uint8)) # phenocycler outputs are uint8
        o = model.OME()
        o.images.append(
            model.Image(
                id='Image:0',
                pixels=model.Pixels(
                    dimension_order='XYCZT',
                    size_c=n_channels,
                    size_t=1,
                    size_x=x,
                    size_y=y,
                    size_z=1,
                    type='float',
                    big_endian=False,
                    channels=[model.Channel(id=f'Channel:{i}', name=c) for i, c in enumerate(biomarkers)],
                    physical_size_x=x / PHENOCYCLER_PIXELS_PER_MICRON,
                    physical_size_y=x / PHENOCYCLER_PIXELS_PER_MICRON,
                )
            )
        )

        im = o.images[0]
        im.pixels.physical_size_x_unit = 'µm'
        im.pixels.physical_size_y_unit = 'µm'
        for i in range(len(im.pixels.channels)):
            im.pixels.planes.append(model.Plane(the_c=i, the_t=0, the_z=0))
        im.pixels.tiff_data_blocks.append(model.TiffData(plane_count=len(im.pixels.channels)))
        xml_str = to_xml(o)
        out_tif.overwrite_description(xml_str.encode())
