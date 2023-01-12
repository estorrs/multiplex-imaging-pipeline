import os
import logging
import math
import re
import json
from pathlib import Path
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import shapely
import rasterio
import tifffile

from skimage.morphology import label, remove_small_objects, erosion, binary_dilation, binary_erosion
from skimage.measure import regionprops_table, regionprops
from skimage.segmentation import expand_labels
from shapely import Polygon, Point, STRtree
from shapely.wkt import dumps
from rasterio import features

from mip.utils import extract_ome_tiff

DEBUG_DIR = '/diskmnt/Projects/Users/estorrs/sandbox'


def compute_polsby_popper(area, perimeter):
    try:
        return 4 * math.pi * area / perimeter ** 2
    except ZeroDivisionError:
        return math.nan


def get_regionprops_df(labeled_img,
                       feats=('label', 'bbox', 'area', 'centroid')):
    regions = regionprops_table(labeled_img, properties=feats)
    region_df = pd.DataFrame(regions)
    region_df = region_df.set_index('label')
    return region_df


def get_morphology_masks(labeled_img, myoepi_dist=40, boundary_dist=150):
    #expand by ~1 micron to avoid errors
    labeled_img = expand_labels(labeled_img, distance=2)
    mask = labeled_img > 0

    # myoepithelial
    outer_labeled = expand_labels(labeled_img, distance=myoepi_dist)
    myoepi_labeled = np.zeros_like(labeled_img)
    outer_mask = outer_labeled>0
    myoepi_mask = outer_mask ^ mask
    myoepi_labeled[myoepi_mask] = outer_labeled[myoepi_mask]

    # boundary
    outer_labeled = expand_labels(labeled_img, distance=boundary_dist)
    boundary_labeled = np.zeros_like(labeled_img)
    outer_mask = outer_labeled>0
    boundary_mask = outer_mask ^ mask
    boundary_labeled[boundary_mask] = outer_labeled[boundary_mask]

    labeled_dict = {
        'region': labeled_img,
        'myoepi': myoepi_labeled,
        'boundary': boundary_labeled,
        'expanded': outer_labeled,
    }
    mask_dict = {
        'region': mask,
        'myoepi': myoepi_mask,
        'boundary': boundary_mask,
        'expanded': outer_mask
    }

    return labeled_dict, mask_dict


def mask_to_polygon(mask, use_biggest=True):
    for i, (shape, _) in enumerate(features.shapes(mask.astype(np.int32), connectivity=4, mask=(mask>0),
                                    transform=rasterio.Affine(1.0, 0, 0, 0, 1.0, 0))):
        ring = shapely.geometry.shape(shape)
        print('ring', i, ring.exterior.length)
    for i, (shape, _) in enumerate(features.shapes(mask.astype(np.int32), connectivity=4, mask=(mask>0),
                                    transform=rasterio.Affine(1.0, 0, 0, 0, 1.0, 0))):
        return shapely.geometry.shape(shape)


def generate_labeled_expansion_grid(
        mask, parallel_step=50, perp_steps=5, expansion=80, grouping_dist=10):
    regions = regionprops(label(mask))
    assert len(regions) == 1
    print('perimeter', regions[0].perimeter)
    tifffile.imsave(os.path.join(DEBUG_DIR, str(regions[0].perimeter) + 'region_mask.npy'), mask)
    # np.save(os.path.join(DEBUG_DIR, 'region_mask.npy'), mask)
    ring = mask_to_polygon(mask).exterior
    print('first ring length', ring.length)
    
    # test for left and right
    left_ring = ring.parallel_offset(5, side='left')
    right_ring = ring.parallel_offset(5, side='right')
    if left_ring.length > right_ring.length:
        side = 'left'
    else:
        side = 'right'

    # logging.info(f'using {side} side to expand')

    perp_step = expansion // perp_steps
    rings = []
    for i in range(perp_steps):
        new = ring.parallel_offset(perp_step * i, side=side)
        rings.append(new)
    rings[0] = ring # fails on offset of zero so replacing
    
    # generate grid
    polys = []
    for i in range(len(rings) - 1):
        inner_ring = rings[i]
        outer_ring = rings[i + 1]
        parallel_steps = int(inner_ring.length // parallel_step)
        first_inner_pt, last_inner_pt, first_outer_pt, last_outer_pt = None, None, None, None
        for j in range(parallel_steps + 1):
            inner_pt1 = inner_ring.interpolate(j * parallel_step)
            inner_pt2 = inner_ring.interpolate((j + 1) * parallel_step)

            inner_pts = [inner_ring.interpolate(x)
                         for x in range(j * parallel_step, (j + 1) * parallel_step, 1)]
            p1 = shapely.shortest_line(inner_pt1, outer_ring)
            p2 = shapely.shortest_line(inner_pt2, outer_ring)
            if p1 is not None and p2 is not None:
                outer_pt1 = Point(p1.coords[1]) if j else outer_ring.interpolate(0.)
                outer_pt2 = Point(p2.coords[1])
                d1, d2 = outer_ring.project(outer_pt1), outer_ring.project(outer_pt2)
                outer_pts = [outer_ring.interpolate(x)
                            for x in range(int(d2), int(d1) - 1, -1)]

                poly_pts = inner_pts
                poly_pts += outer_pts
                poly_pts = [pt.coords[0] for pt in poly_pts]
                poly = Polygon(poly_pts)
                polys.append(poly)

                if not j:
                    first_inner_pt = inner_pt1
                    first_outer_pt = outer_pt1
                elif j == parallel_steps:
                    last_inner_pt = inner_pt2
                    last_outer_pt = outer_pt2

        # close outer rings
        poly_pts = [last_inner_pt, last_outer_pt, first_outer_pt, first_inner_pt]
        if any([x is None for x in poly_pts]):
            pass
            # logging.info(f'arc {j} closing failed: last inner - {last_inner_pt}, last outer - {last_outer_pt}, first outer - {first_outer_pt}, first inner - {first_inner_pt}')
        else:
            poly_pts = [pt.coords[0] for pt in poly_pts]
            poly = Polygon(poly_pts)
            polys.append(poly)

    # generate tile groupings
    tree = STRtree(polys)
    inner_ring = rings[0]
    outer_ring = rings[-1]
    parallel_steps = int(inner_ring.length // grouping_dist)
    groups = []
    group_lines = []
    for j in range(parallel_steps):
        dist = (j * grouping_dist) + (grouping_dist // 2)
        inner_pt = inner_ring.interpolate(dist)
        line = shapely.shortest_line(inner_pt, outer_ring)
        idxs = tree.query(line, predicate="intersects")
        groups.append(idxs + 1) # adding 1 to match label mask
        group_lines.append(line)

    labeled_grid = np.zeros(mask.shape, dtype=np.int32)
    for i, poly in enumerate(polys):
        m = rasterio.features.rasterize([poly], out_shape=mask.shape)
        labeled_grid[m!=0] = i + 1
    return labeled_grid, polys, groups, group_lines, rings


def generate_grid_props(region_id, labeled_dict, region_to_bbox,
                        parallel_step=50, perp_steps=5,
                        expansion=80, grouping_dist=10):
    print('region id', region_id)
    r1, c1, r2, c2 = region_to_bbox[region_id]
    print('bbox', (r1, r2, c1, c2))
    mask = (labeled_dict['region'][r1:r2, c1:c2]==region_id).astype(np.int32)
    print('mask count', np.count_nonzero(mask))
    labeled_grid, polys, groups, group_lines, rings = generate_labeled_expansion_grid(
        mask, parallel_step=parallel_step, perp_steps=perp_steps,
        expansion=expansion, grouping_dist=grouping_dist
    )
    print('ring length', rings[0].length)

    return labeled_grid, polys, groups, group_lines, rings


def generate_grid_dict(region_to_bbox, labeled_dict,
                       parallel_step=50, perp_steps=5,
                       expansion=80, grouping_dist=10):
    grid_dict = {}
    for region_id in region_to_bbox.keys():
        logging.info(f'Generating grid properties for region {region_id}')
        print('region id', region_id)
        labeled_grid, polys, groups, group_lines, rings = generate_grid_props(
            region_id, labeled_dict, region_to_bbox,
            parallel_step=parallel_step, perp_steps=perp_steps,
            expansion=expansion, grouping_dist=grouping_dist
        )
        print('region id', region_id)
        grid_dict[region_id] = {
            'labeled_grid': labeled_grid,
            'polys': polys,
            'groups': groups,
            'group_lines': group_lines,
            'rings': rings
        }
    return grid_dict


def has_valid_coords(poly, x_max, y_max):
    xs, ys = poly.exterior.xy
    xs, ys = np.asarray(xs), np.asarray(ys)

    if not xs.max() > 0:
        return False
    if not ys.max() > 0:
        return False
    if not xs.min() < x_max:
        return False
    if not ys.min() < y_max:
        return False

    return True


def filter_polys(grid_dict, region_to_bbox, area_thresh=200, group_line_thresh=100):
    filtered_grid_dict = {}
    for region_id, d in grid_dict.items():
        _, _, r2, c2 = region_to_bbox[region_id]
        remove_idxs = [i for i, p in enumerate(d['polys'])
                       if p.area > area_thresh]
        pool = set(remove_idxs)
        keep_idxs = [i for i in range(len(d['polys'])) if i not in pool]
        filtered_polys = [d['polys'][i] for i in keep_idxs]

        filtered_labeled_grid = d['labeled_grid'].copy()
        for i in remove_idxs:
            filtered_labeled_grid[d['labeled_grid']==i] = 0

        idxs = [i for i, p in enumerate(d['group_lines']) if p is not None and p.length <= group_line_thresh]
        filtered_group_lines = [d['group_lines'][i] for i in idxs]
        filtered_groups = [d['groups'][i] for i in idxs]

        filtered_grid_dict[region_id] = {
            'polys': filtered_polys,
            'labeled_grid': filtered_labeled_grid,
            'group_lines': filtered_group_lines,
            'groups': filtered_groups,
            'rings': d['rings']
        }
    return filtered_grid_dict


def get_positive_grid(labeled_grid, df):
    new = np.zeros_like(labeled_grid)
    for l in df[df['is_positive']].index:
        new[labeled_grid==l] = l
    return new


def get_num_pieces(pos_grid):
    new = pos_grid!=0
    new = binary_dilation(new)

    labeled = label(new)
    return np.unique(labeled).shape[0] - 1


def get_thicknesses(labeled_grid, df, groups):
    pool = set(np.unique(labeled_grid))
    valid = []
    for group in groups:
        # check to make sure it's not out of bounds
        keep = True
        for idx in group:
            if idx not in pool:
                keep = False
                break
        valid.append(keep)

    groups = [g for g, keep in zip(groups, valid) if keep]

    thicknesses = []
    for group in groups:
        thicknesses.append(df.loc[group]['is_positive'].sum())

    return np.asarray(thicknesses)


def generate_grid_metrics(region_id, grid_dict, region_to_bbox, channel_to_img,
                       channel_to_thresh={'SMA': 1500, 'Podoplanin': 5000}):
    r1, c1, r2, c2 = region_to_bbox[region_id]
    d = grid_dict[region_id]
    labeled_grid = d['labeled_grid']
    groups = d['groups']

    channel_to_props = {}
    channel_to_df = {}
    channel_to_pos_grid = {}
    for channel, thresh in channel_to_thresh.items():
        img = channel_to_img[channel][r1:r2, c1:c2]
        df = pd.DataFrame(regionprops_table(labeled_grid, intensity_image=img,
                                            properties=('label', 'area', 'intensity_mean')))
        df = df.set_index('label')
        df['is_positive'] = df['intensity_mean']>=thresh
        channel_to_df[channel] = df

        thicknesses = get_thicknesses(labeled_grid, df, groups)
        pos_grid = get_positive_grid(labeled_grid, df)
        channel_to_pos_grid[channel] = pos_grid
        channel_to_props[channel] = {
            'num_grids': df.shape[0],
            'total_area': df['area'].sum(),
            'mean_intensity': df['intensity_mean'].mean(),
            'num_pieces': get_num_pieces(pos_grid),
            'fraction_positive': df['is_positive'].sum() / df.shape[0],
            'integrity': 1 - (np.count_nonzero(thicknesses==0) / len(thicknesses)) if len(thicknesses) else np.nan,
            'mean_thickness': thicknesses.mean(),
            'mean_thickness_nonbreaks': thicknesses[thicknesses!=0].mean() if np.count_nonzero(thicknesses!=0) else 0.
        }

    combined_props = {}
    pos_grids = np.asarray(list(channel_to_pos_grid.values())) # (num_channels, H, W)
    combined_props['num_pieces'] = get_num_pieces(pos_grids.max(axis=0))

    ls = np.asarray([df['is_positive'].to_list() for df in channel_to_df.values()]) # (num_channels, num_grids)
    is_positive = ls.sum(axis=0)
    df = next(iter(channel_to_df.values())).copy()
    df['is_positive'] = is_positive
    thicknesses = get_thicknesses(labeled_grid, df, groups)
    combined_props.update({
        'overlap_fraction': np.count_nonzero(ls.sum(axis=0)==ls.shape[0]) / ls.shape[1],
        'fraction_positive': df['is_positive'].sum() / df.shape[0],
        'integrity': 1 - (np.count_nonzero(thicknesses==0) / len(thicknesses)) if len(thicknesses) else np.nan,
        'mean_thickness': thicknesses.mean(),
        'mean_thickness_nonbreaks': thicknesses[thicknesses!=0].mean() if np.count_nonzero(thicknesses!=0) else 0.
    })

    channel_to_pos_grid['combined'] = pos_grids.max(axis=0)
    channel_to_props['combined'] = combined_props

    return channel_to_df, channel_to_props, channel_to_pos_grid


def add_polsby_popper(props_dict):
    for name, df in props_dict.items():
        df['compactness'] = [compute_polsby_popper(area, per)
                             for area, per in zip(df['area'], df['perimeter'])]                    


def add_cell_fractions(spatial_feats_df, props_dict, cell_type_col='cell_type'):
    cell_types = sorted(set(spatial_feats_df[cell_type_col]))
    for name, df in props_dict.items():
        data = []
        for region in df.index:
            ls = spatial_feats_df[spatial_feats_df[f'{name}_region_id']==region][cell_type_col]
            counts = Counter(ls)
            data.append([counts.get(ct, 0) / (len(ls) + 1) for ct in cell_types])
        frac_df = pd.DataFrame(
            data=data, columns=[f'{cell_type_col}_{ct}_cell_fraction'
                                for ct in cell_types],
            index=df.index.to_list())
        props_dict[name] = pd.concat((df, frac_df), axis=1)


def add_pixel_metrics(labeled_dict, props_dict, region_to_bbox, channel_to_img, channel_to_thresh):
    channels = list(channel_to_thresh.keys())
    grid = np.meshgrid(channels, channels)
    tups = [tuple(sorted([x, y])) for x, y in zip(grid[0].flatten(), grid[1].flatten()) if x!=y]
    combos = sorted(set(tups))

    for name, df in props_dict.items():
        img = labeled_dict[name]
        data, cols = [], []
        for region_id in df.index:
            ls = []
            metric_dict = {}
            r1, c1, r2, c2 = region_to_bbox[region_id]
            cropped = img[r1:r2, c1:c2]
            m = cropped==region_id
            channel_to_thresh_mask = {c:channel_to_img[c][r1:r2, c1:c2] >= t for c, t in channel_to_thresh.items()}
            for k, v in channel_to_thresh_mask.items():
                v[~m] = 0
                channel_to_thresh_mask[k] = v
            for channel in channels:
                m1 = channel_to_thresh_mask[channel]
                pos_fraction = np.count_nonzero(m1) / np.count_nonzero(m)
                metric_dict[f'intensity_positive_fraction_{channel}'] = pos_fraction

            for channel_1, channel_2 in combos:
                m1 = channel_to_thresh_mask[channel_1]
                m2 = channel_to_thresh_mask[channel_2]
                both = (m1>0) & (m2>0)
                pos_fraction = np.count_nonzero(both) / np.count_nonzero(m)
                metric_dict[f'intensity_overlap_{channel_1}_{channel_2}'] = pos_fraction
            
            cols, ls = zip(*metric_dict.items())
            data.append(ls)
        new = pd.DataFrame(data=data, columns=cols, index=df.index)
        df = pd.concat((df, new), axis=1)
        props_dict[name] = df



def get_grid_metrics(labeled_dict, region_to_bbox, channel_to_img,
                     parallel_step=50, perp_steps=5,
                     expansion=80, grouping_dist=10,
                     area_thresh=200, group_line_thresh=100,
                     channel_to_thresh={'SMA': 1500, 'Podoplanin': 5000}):
    logging.info('Genreating grid properties')
    grid_dict = generate_grid_dict(
        region_to_bbox, labeled_dict,
        parallel_step=parallel_step, perp_steps=perp_steps, expansion=expansion,
        grouping_dist=grouping_dist)
    filtered_grid_dict = filter_polys(
        grid_dict, region_to_bbox,
        area_thresh=area_thresh, group_line_thresh=group_line_thresh)

    logging.info('Calculating grid metrics')
    region_to_results = {}
    for region_id in filtered_grid_dict:
        try:
            channel_to_df, channel_to_props, channel_to_pos_grid = generate_grid_metrics(
                region_id, filtered_grid_dict, region_to_bbox, channel_to_img,
                channel_to_thresh=channel_to_thresh
            )
            region_to_results[region_id] = {
                'channel_to_df': channel_to_df,
                'channel_to_props': channel_to_props,
                'channel_to_pos_grid': channel_to_pos_grid
            }
        except ZeroDivisionError:
            logging.info(f'zero division error in region {region_id}')
            logging.info(f'channel all dataframe shape: {channel_to_df.shape}')

    return region_to_results, grid_dict, filtered_grid_dict


def combine_metrics(props_dict, grid_metric_results):
    combined_df = None
    for k, df in props_dict.items():
        temp = df.copy()
        temp.columns = [f'{k}_{c}' for c in temp.columns]
        if combined_df is None:
            combined_df = temp
        else:
            combined_df = pd.concat((combined_df, temp), axis=1)

    if grid_metric_results is None:
        return combined_df

    cols = [f'{channel}_{k}'
            for channel, meta in next(iter(grid_metric_results.values()))['channel_to_props'].items()
            for k, v in meta.items()]
    idxs = []
    data = []
    for region_id in combined_df.index:
        d = grid_metric_results[region_id]
        idxs.append(region_id)
        row = []
        for col in cols:
            channel = re.sub(r'^([^_]+)_.*$', r'\1', col)
            metric = re.sub(r'^[^_]+_(.*)$', r'\1', col)
            val = d['channel_to_props'][channel][metric]
            row.append(val)
        data.append(row)

    grid_df = pd.DataFrame(data=data, columns=cols, index=idxs)
    grid_df.columns = [f'region_grid_metrics_{c}' for c in grid_df.columns]
    combined_df = pd.concat((combined_df, grid_df), axis=1)

    return combined_df


def save_results(output_dir, metric_df, labeled_dict, mask_dict,
                 grid_metric_results, grid_dict, filtered_grid_dict):
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # metrics
    metric_df.to_csv(os.path.join(output_dir, 'region_metrics.txt'), sep='\t')

    # masks
    to_save = {k:{
                   'labeled': labeled_dict[k],
                   'mask': mask_dict[k],
              } for k in labeled_dict.keys()}
    np.save(os.path.join(output_dir, 'regions.npy'), to_save)

    if grid_metric_results is not None: # only save if we calculated grid stuff
        # non shapely grid regions
        region_ids = metric_df.index.to_list()
        to_save = {k:{
                    'grid_binary_positive_regions': grid_metric_results[k]['channel_to_pos_grid'],
                    'grid_labeled_regions_filtered': filtered_grid_dict[k]['labeled_grid'],
                    'grid_labeled_regions_unfiltered': grid_dict[k]['labeled_grid'],
                } for k in region_ids}
        np.save(os.path.join(output_dir, 'grid_regions.npy'), to_save)

        # grid polys and lines
        to_save = {k:{
                        'polygons_filtered': {k:list(dumps(v))
                                            for k, v in filtered_grid_dict[k].items()
                                            if k in ['polys', 'group_lines', 'rings']},
                        'polygons_unfiltered': {k:list(dumps(v))
                                                for k, v in grid_dict[k].items()
                                                if k in ['polys', 'group_lines', 'rings']},
                    } for k in region_ids}
        json.dump(to_save, open(os.path.join(output_dir, 'grid_polygons.json'), 'w'))

        # the leftovers
        to_save = {k:{
                    'grid_channel_to_df': {
                            c:df.to_dict()
                            for c, df in grid_metric_results[k]['channel_to_df'].items()},
                } for k in region_ids}

        json.dump(to_save, open(os.path.join(output_dir, 'other_metadata.json'), 'w'))


def generate_region_metrics(
        spatial_features_df, ome_fp, region_mask_fp, output_dir,
        y_col='centroid_row', x_col='centroid_col', cell_metadata_cols=('cell_type',),
        boundary_dist=150, parallel_step=50, perp_steps=5,
        expansion=80, grouping_dist=10, area_thresh=200, group_line_thresh=100,
        channel_to_thresh_grid={'SMA': 1500}, channel_to_thresh_pixel={'SMA': 1500},
        min_region_size=None, max_region_size=None, calculate_grid_metrics=True):
    
    initial_region_mask = tifffile.imread(region_mask_fp) > 0
    logging.info(f'mask loaded. Size: {initial_region_mask.shape}')

    logging.info('making morphology masks')
    initial_labeled = label(initial_region_mask)
    labeled_dict, mask_dict = get_morphology_masks(
        initial_labeled, boundary_dist=boundary_dist, myoepi_dist=expansion)
    region_labeled = labeled_dict['region'].copy()
    # region_mask = mask_dict['region'].copy()
    # initial_region_df = get_regionprops_df(region_labeled)

    # if min_region_size is not None:
    #     logging.info('removing small regions')
    #     region_mask = remove_small_objects(region_mask, min_size=min_region_size)
    #     region_labeled[~region_mask] = 0

    region_df = get_regionprops_df(region_labeled)
    n_regions = region_df.shape[0]
    logging.info(f'Labeled region props generated. {n_regions} total regions.')

    # # # allows for max size threshold, speeds up runs for debugging purposes
    # if max_region_size is not None:
    #     logging.info(f'filtering regions with area greater than {max_region_size}')
    #     remove = [l for l, area in zip(region_df.index, region_df['area'])
    #               if area > max_region_size]

    #     for l in remove:
    #         r1, c1, r2, c2 = (region_df.loc[l, 'bbox-0'], region_df.loc[l, 'bbox-1'],
    #                           region_df.loc[l, 'bbox-2'], region_df.loc[l, 'bbox-3'])
    #         cropped = region_labeled[r1:r2, c1:c2] # cropping is much faster
    #         new = cropped.copy()
    #         new[cropped==l] = 0
    #         region_labeled[r1:r2, c1:c2] = new
    #     region_mask = region_labeled > 0

    #     # regenerate region_df
    #     region_df = get_regionprops_df(region_labeled)
    #     n_regions = region_df.shape[0]
    #     logging.info(f'{n_regions} regions remaining after area filtering')

    # # resychronize masks and labeled imgs
    # logging.info('resynchronizing labeled masks')
    # region_ids = region_df.index.to_list()
    # previous_ids = initial_region_df.index.to_list()
    # to_remove = set(previous_ids) - set(region_ids)
    # for name, img in labeled_dict.items():
    #     for l in to_remove:
    #         r1, c1, r2, c2 = (initial_region_df.loc[l, 'bbox-0'], initial_region_df.loc[l, 'bbox-1'],
    #                         initial_region_df.loc[l, 'bbox-2'], initial_region_df.loc[l, 'bbox-3'])
    #         cropped = img[r1:r2, c1:c2] # cropping is much faster
    #         new = cropped.copy()
    #         new[cropped==l] = 0
    #         img[r1:r2, c1:c2] = new
    #     labeled_dict[name] = img
    #     mask_dict[name] = img > 0

    props_dict = {}
    for name, img in labeled_dict.items():
        props_dict[name] = get_regionprops_df(
            img, ('label', 'bbox', 'area', 'perimeter', 'centroid', 'eccentricity', 'extent'))

    region_to_bbox = {l:(r1, c1, r2, c2)
                      for l, r1, c1, r2, c2 in zip(
                          props_dict['expanded'].index,
                          props_dict['expanded']['bbox-0'],
                          props_dict['expanded']['bbox-1'],
                          props_dict['expanded']['bbox-2'],
                          props_dict['expanded']['bbox-3'])}

    logging.info('extracting ome.tiff')
    total_channels = set(channel_to_thresh_grid.keys())
    total_channels.update(set(channel_to_thresh_pixel.keys()))
    channel_to_img = extract_ome_tiff(ome_fp, channels=list(total_channels))
    channels = list(channel_to_img.keys())
    logging.info(f'channels in ome.tiff: {channels}')
    assert next(iter(channel_to_img.values())).shape == labeled_dict['region'].shape

    for name, img in labeled_dict.items():
        spatial_features_df[f'{name}_region_id'] = [
            img[int(r), int(c)]
            for r, c in zip(spatial_features_df[y_col], spatial_features_df[x_col])]

    logging.info('calculating polsby popper compactness')
    # compactness
    add_polsby_popper(props_dict)

    logging.info('calculating border region cellular metadata fractions')
    # cell fractions
    for col in cell_metadata_cols:
        add_cell_fractions(spatial_features_df, props_dict, cell_type_col=col)

    logging.info('adding pixel level metrics')
    add_pixel_metrics(labeled_dict, props_dict, region_to_bbox, channel_to_img, channel_to_thresh_pixel)

    if not calculate_grid_metrics: # stop early if we are not calculating grid metrics
        logging.info('skipping grid-based features')
        metric_df = combine_metrics(props_dict, None)
        logging.info('saving results')
        save_results(output_dir, metric_df, labeled_dict, mask_dict,
                     None, None, None)
        logging.info('finished saving')
        return None

    # grid metrics
    logging.info('adding grid metrics')
    grid_metric_results, grid_dict, filtered_grid_dict = get_grid_metrics(
        labeled_dict, region_to_bbox, channel_to_img,
        parallel_step=parallel_step, perp_steps=perp_steps,
        expansion=expansion, grouping_dist=grouping_dist, area_thresh=area_thresh,
        group_line_thresh=group_line_thresh, channel_to_thresh=channel_to_thresh_grid)

    logging.info('formatting metrics')
    metric_df = combine_metrics(props_dict, grid_metric_results)

    logging.info('saving results')
    save_results(output_dir, metric_df, labeled_dict, mask_dict,
                 grid_metric_results, grid_dict, filtered_grid_dict)
    logging.info('finished saving')
    