import logging
import re
import os

import numpy as np
import tifffile
import torchvision.transforms.functional as TF
import torch
from ome_types import from_tiff, from_xml, to_xml, model
from ome_types.model.simple_types import UnitsLength
from einops import rearrange
from skimage.transform import rescale

from multiplex_imaging_pipeline.utils import listfiles

# an effort to have consistent marker names in all level 2 ome.tiffs
# getting rid of it for now because it messes with some of the qc
d = {
}
CHANNEL_MAP = {v:k for k, vs in d.items() for v in vs}

# calculated from andrew shinkle measurement on HT206 test slide
PHENOCYCLER_PIXELS_PER_MICRON = 1.9604911906033102

# in 20220106_HT206B1_H1_A1 and A4_reg001 metadata
CODEX_OLD_PIXELS_PER_MICRON = 2.6495274


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


def generate_ome_from_tifs(fps, output_fp, platform='codex', bbox=None, pixel_type='uint8', subresolutions=4):
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
    n_channels = len(keep)

    biomarkers = []
    data = None
    for i, fp in enumerate(fps):
        if i in keep_idxs:
            img = tifffile.imread(fp)
            biomarker = mapping(names[i])
            biomarkers.append(biomarker)
            if bbox is not None:
                r1, r2, c1, c2 = bbox
                y = r2 - r1
                x = c2 - c1
                img = img[r1:r2, c1:c2]
            img = img.astype(np.float32)
            img /= img.max()
            if pixel_type == 'uint8':
                img *= 255
                img = img.astype(np.uint8)
            else:
                img *= 65535
                img = img.astype(np.uint16) 

            if data is None:
                if bbox is None:
                    data = np.zeros((img.shape[1], img.shape[0], n_channels, 1, 1),
                                    dtype=np.uint8 if pixel_type=='uint8' else np.uint16)
                else:
                    logging.info(f'bbox detected, cropping to {bbox}')
                    data = np.zeros((bbox[3] - bbox[2], bbox[1] - bbox[0], n_channels, 1, 1),
                                    dtype=np.uint8 if pixel_type=='uint8' else np.uint16)
            else:
                data[..., i, 0, 0] = np.swapaxes(img, 0, 1)
    
    o = model.OME()
    o.images.append(
        model.Image(
            id='Image:0',
            pixels=model.Pixels(
                dimension_order='XYCZT',
                size_c=n_channels,
                size_t=1,
                size_x=data.shape[0],
                size_y=data.shape[1],
                size_z=data.shape[3],
                type=pixel_type,
                big_endian=False,
                channels=[model.Channel(id=f'Channel:{i}', name=c) for i, c in enumerate(biomarkers)],
                physical_size_x=1 / PHENOCYCLER_PIXELS_PER_MICRON,
                physical_size_y=1 / PHENOCYCLER_PIXELS_PER_MICRON,
                physical_size_x_unit='µm',
                physical_size_y_unit='µm'
            )
        )
    )

    im = o.images[0]
    for i in range(len(im.pixels.channels)):
        im.pixels.planes.append(model.Plane(the_c=i, the_t=0, the_z=0))
    im.pixels.tiff_data_blocks.append(model.TiffData(plane_count=len(im.pixels.channels)))

    with tifffile.TiffWriter(output_fp, ome=True, bigtiff=True) as out_tif:
        opts = {
            'compression': 'LZW',
        }
        out_tif.write(
            rearrange(data, 'x y c z t -> t c y x z'),
            subifds=subresolutions,
            **opts
        )
        for level in range(subresolutions):
            mag = 2**(level + 1)
            x = torch.tensor(rearrange(data[..., 0, 0], 'w h c -> c h w'))
            sampled = rearrange(
                TF.resize(x, (int(x.shape[-2] / mag), int(x.shape[-1] / mag)), antialias=True),
                'c h w -> w h c 1 1'
            )
            sampled = sampled.numpy().astype(np.uint8)
            out_tif.write(
                rearrange(sampled, 'x y c z t -> t c y x z'),
                subfiletype=1,
                **opts
            )
        xml_str = to_xml(o)
        out_tif.overwrite_description(xml_str.encode())


def generate_ome_from_qptiff(qptiff_fp, output_fp, bbox=None, pixel_type='uint8', subresolutions=4):
    """
    Generate an HTAN compatible ome tiff from qptiff output by codex phenocycler
    """
    tf = tifffile.TiffFile(qptiff_fp)
    s = tf.series[0] # full res tiffs are in the first series
    n_channels = s.get_shape()[0]
    shape = s.pages[0].shape

    logging.info(f'{n_channels} total markers')
    logging.info(f'writing as {pixel_type}')

    if bbox is None:
        data = np.zeros((shape[1], shape[0], n_channels, 1, 1),
                        dtype=np.uint8 if pixel_type=='uint8' else np.uint16)
    else:
        logging.info(f'bbox detected: {bbox}')
        data = np.zeros((bbox[3] - bbox[2], bbox[1] - bbox[0], n_channels, 1, 1),
                        dtype=np.uint8 if pixel_type=='uint8' else np.uint16)
    biomarkers = []
    for i, p in enumerate(s.pages):
        img = p.asarray()
        d = tifffile.xml2dict(p.description)['PerkinElmer-QPI-ImageDescription']
        biomarker = d['Biomarker']
        logging.info(f'loading {biomarker}')
        biomarkers.append(biomarker)
        if bbox is not None:
            r1, r2, c1, c2 = bbox
            img = img[r1:r2, c1:c2]

        if pixel_type=='uint8':
            img = img.astype(np.uint8)
        else:
            img = img.astype(np.float32)
            img /= img.max()
            img *= 65535
            img = img.astype(np.uint16) 
            logging.info('writing as uint16')

        data[..., i, 0, 0] = np.swapaxes(img, 0, 1)

    o = model.OME()
    o.images.append(
        model.Image(
            id='Image:0',
            pixels=model.Pixels(
                dimension_order='XYCZT',
                size_c=n_channels,
                size_t=1,
                size_x=data.shape[0],
                size_y=data.shape[1],
                size_z=data.shape[3],
                type=pixel_type,
                big_endian=False,
                channels=[model.Channel(id=f'Channel:{i}', name=c) for i, c in enumerate(biomarkers)],
                physical_size_x=1 / PHENOCYCLER_PIXELS_PER_MICRON,
                physical_size_y=1 / PHENOCYCLER_PIXELS_PER_MICRON,
                physical_size_x_unit='µm',
                physical_size_y_unit='µm'
            )
        )
    )

    im = o.images[0]
    for i in range(len(im.pixels.channels)):
        im.pixels.planes.append(model.Plane(the_c=i, the_t=0, the_z=0))
    im.pixels.tiff_data_blocks.append(model.TiffData(plane_count=len(im.pixels.channels)))

    with tifffile.TiffWriter(output_fp, ome=True, bigtiff=True) as out_tif:
        opts = {
            'compression': 'LZW',
        }
        logging.info(f'writting full res image')
        out_tif.write(
            rearrange(data, 'x y c z t -> t c y x z'),
            subifds=subresolutions,
            **opts
        )
        for level in range(subresolutions):
            mag = 2**(level + 1)
            logging.info(f'writting subres {mag}')
            x = torch.tensor(rearrange(data[..., 0, 0], 'w h c -> c h w'))
            sampled = rearrange(
                TF.resize(x, (int(x.shape[-2] / mag), int(x.shape[-1] / mag)), antialias=True),
                'c h w -> w h c 1 1'
            )
            sampled = sampled.numpy().astype(np.uint8)
            out_tif.write(
                rearrange(sampled, 'x y c z t -> t c y x z'),
                subfiletype=1,
                **opts
            )
        xml_str = to_xml(o)
        out_tif.overwrite_description(xml_str.encode())


def generate_ome_from_codex_imagej_tif(tif_fp, output_fp, bbox=None, pixel_type='uint8', subresolutions=4):
    """
    Generate an HTAN compatible ome tiff from qptiff output by codex phenocycler
    """
    tf = tifffile.TiffFile(tif_fp)
    channels = tf.imagej_metadata['Labels']
    keep_idxs, keep_channels = zip(*[(i, c) for i, c in enumerate(channels)
                     if c.lower()!='blank'
                     if c.lower()!='empty'
                     if c[:4]!='DAPI'])
    keep_idxs = list(keep_idxs)
    keep_channels = list(keep_channels)

    n_channels = len(keep_channels) + 1
    logging.info(f'{n_channels} total markers')
    logging.info(f'writing as {pixel_type}')

    img = tifffile.imread(tif_fp)

    if bbox is not None:
        logging.info(f'bbox detected: {bbox}')
        r1, r2, c1, c2 = bbox
        img = img[..., r1:r2, c1:c2]

    if pixel_type=='uint8':
        img = img.astype(np.float32)
        img /= img.max()
        img *= 255.
        img = img.astype(np.uint8)
    else:
        img = img.astype(np.float32)
        img /= img.max()
        img *= 65535
        img = img.astype(np.uint16)

    data = rearrange(img, 'n_cycles n_per_cycle h w -> (n_cycles n_per_cycle) h w')
    data = data[keep_idxs]

    # add DAPI
    data = np.concatenate((np.expand_dims(img[0, 0], 0), data))
    keep_channels.insert(0, 'DAPI')
    logging.info(f'image channels: {keep_channels}')

    # add z and t dims
    data = rearrange(data, 'c h w -> w h c 1 1')

    o = model.OME()
    o.images.append(
        model.Image(
            id='Image:0',
            pixels=model.Pixels(
                dimension_order='XYCZT',
                size_c=n_channels,
                size_t=1,
                size_x=data.shape[0],
                size_y=data.shape[1],
                size_z=data.shape[3],
                type=pixel_type,
                big_endian=False,
                channels=[model.Channel(id=f'Channel:{i}', name=c) for i, c in enumerate(keep_channels)],
                physical_size_x=1 / CODEX_OLD_PIXELS_PER_MICRON,
                physical_size_y=1 / CODEX_OLD_PIXELS_PER_MICRON,
                physical_size_x_unit='µm',
                physical_size_y_unit='µm'
            )
        )
    )

    im = o.images[0]
    for i in range(len(im.pixels.channels)):
        im.pixels.planes.append(model.Plane(the_c=i, the_t=0, the_z=0))
    im.pixels.tiff_data_blocks.append(model.TiffData(plane_count=len(im.pixels.channels)))

    with tifffile.TiffWriter(output_fp) as out_tif:
        opts = {
            'compression': 'LZW',
        }
        logging.info(f'writting full res image')
        out_tif.write(
            rearrange(data, 'x y c z t -> t c y x z'),
            subifds=subresolutions,
            **opts
        )
        for level in range(subresolutions):
            mag = 2**(level + 1)
            logging.info(f'writting subres {mag}')
            x = torch.tensor(rearrange(data[..., 0, 0], 'w h c -> c h w'))
            sampled = rearrange(
                TF.resize(x, (int(x.shape[-2] / mag), int(x.shape[-1] / mag)), antialias=True),
                'c h w -> w h c 1 1'
            )
            sampled = sampled.numpy().astype(np.uint8)
            out_tif.write(
                rearrange(sampled, 'x y c z t -> t c y x z'),
                subfiletype=1,
                **opts
            )
        xml_str = to_xml(o)
        out_tif.overwrite_description(xml_str.encode())