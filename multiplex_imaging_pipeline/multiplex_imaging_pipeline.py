import argparse
import logging

import numpy as np
import pandas as pd
import tifffile

import multiplex_imaging_pipeline.utils as utils
from multiplex_imaging_pipeline.ome import generate_ome_from_tifs, generate_ome_from_qptiff, generate_ome_from_codex_imagej_tif
from multiplex_imaging_pipeline.spatial_features import get_spatial_features
# from multiplex_imaging_pipeline.region_analysis import generate_region_metrics
from multiplex_imaging_pipeline.segmentation import segment_cells

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('mode', type=str,
    choices=['make-ome', 'segment-ome', 'generate-spatial-features', 'show-channels'],
    help='Which task mip is to execute.')


###################
## show-channels ##
###################

parser.add_argument('--sep', type=str, default='\n',
    help='Seperator between channel names. Defaults to newline character (i.e. channels are displayed on seperate lines)')

##############
## make-ome ##
##############
parser.add_argument('--input-tif', type=str,
    help='Used in make-ome mode. Either directory of stitched tif files that will be combined into a single ome.tiff file, a multichannel .tif (for original codex platform), or a .qptiff (phenocycler platform).')

parser.add_argument('--output-filepath', type=str,
    help='Location to write ome.tiff file')

parser.add_argument('--platform', type=str,
    choices=['codex', 'phenocycler', 'raw'], default='phenocycler',
    help='Which platform produced the input images. phenocycler assumes a .qptiff from the akoya phenocycler platform. codex assumes a multichannel .tif output by the original akoya codex platform. raw will save a directory of .tifs together into a multiplex image.')

parser.add_argument('--bbox', type=str,
    help='If desired, bbox in to crop image with. Must be the following format: "top,bottom,left,right"')

# #################
# ## segment-ome ##
# #################
# parser.add_argument('--input-tif', type=str,
#     help='ome.tiff file to segment')

parser.add_argument('--output-prefix', type=str, default='output',
    help='Output prefix to use when writing cell segmentation files. Two files will be written: *_cell_segmentation.tif and *_nuclei_segmentation.tif. Default is "output". For example, if --output-prefix is "path/to/out/directory/sample" then the two output files will be named path/to/out/directory/sample_cell_segmentation.tif and path/to/out/directory/sample_nuclei_segmentation.tif.')

parser.add_argument('--split-size', type=int, default=25000,
    help='If image width or height is larger than --split-size, then image will be split into multiple pieces for segmentation and stitched back together. Decrease if running into memory issues.')

parser.add_argument('--nuclei-channels', type=str, default='DAPI',
    help='List of nuclei markers to use during segmentation. Must be the following format: "MARKER_NAME,MARKER_NAME,....". Default is "DAPI".')

parser.add_argument('--membrane-channels', type=str, default='Pan-Cytokeratin,E-cadherin,CD45,CD8,CD3e,Vimentin,SMA,CD31,C20',
    help='List of markers to use during membrane segmentation. Must be the following format: "MARKER_NAME,MARKER_NAME,....". Default is "Pan-Cytokeratin,E-cadherin,CD45,CD8,CD3e,Vimentin,SMA,CD31,CD20". Note that image marker names are automatically converted using mip.utils.R_CHANNEL_MAPPING')

###############################
## generate-spatial-features ##
###############################
# parser.add_argument('--input-tif', type=str,
#     help='ome.tiff file to generate spatial features for')

# parser.add_argument('--output-prefix', type=str, default='output',
#     help='Output prefix to use when writing output files. Two files will be written: *_spatial_features.h5ad and *_spatial_features.txt. Default is "output". For example, if --output-prefix is "path/to/out/directory/sample" then the two output files will be named path/to/out/directory/sample_spatial_features.h5ad and path/to/out/directory/sample_spatial_features.txt.')

parser.add_argument('--labeled-image', type=str,
    help='Filepath of labeled cell segmentation tif.')

parser.add_argument('--thresholds', type=str,
    help='Filepath to tab-seperated .txt file containing marker threshold values. The file should have two columns in the following order: <marker>\t<theshold>. Do not include column names/headers in the file.')

# ###############################
# ## generate-region-features ##
# ###############################
# parser.add_argument('--spatial-features', type=str,
#     help='Filepath of a tab-seperated .txt file with columns specifying coordinates cell annotations in slide. First column is cell ID, second and third columns are treated as "x" and "y" coordinates respectively. All following columns are treated as cell metadata features and seperate fractions/metrics will be generated for each feature.')

# parser.add_argument('--regions-mask', type=str,
#     help='Filepath of region mask that will be used when calculating metrics. Mask should be a .tif file and the same height and width as --ome-tiff.')

# parser.add_argument('--channel-thresholds-grid', type=str,
#     help='Filepath of tab-seperated .txt file where the first column is a channel name in the --ome-tiff and the second column is threshold values to use when determining positive grid polygons.')

# parser.add_argument('--channel-thresholds-pixel', type=str,
#     help='Filepath of tab-seperated .txt file where the first column is a channel name in the --ome-tiff and the second column is threshold values to use when determining positive pixels when generating region metrics.')

# parser.add_argument('--output-dir', type=str, default='output',
#     help='Location to write generate-region-features output files.')

# parser.add_argument('--boundary-dist', type=int, default=150,
#     help='Distance (in pixels) to draw boundary around each region.')

# parser.add_argument('--perp-steps', type=int, default=10,
#     help='Number of arcs to generate for region grid when drawing grid polygons.')

# parser.add_argument('--expansion', type=int, default=40,
#     help='Distance (in pixels) of innermost arc to outermost arc.')

# parser.add_argument('--parallel-step', type=int, default=50,
#     help='Step size (in pixels) of steps along arcs to use when drawing grid polygons.')

# parser.add_argument('--breakage-dist', type=int, default=10,
#     help='Distance (in pixels) along innermost arc to use when drawing breakage lines.')

# parser.add_argument('--area_thresh', type=int, default=2000,
#     help='Filter out grid polygons with area greater than area-thresh.')

# parser.add_argument('--breakage-line-thresh', type=int, default=100,
#     help='Filter out grid polygons with area greater than area-thresh.')

# parser.add_argument('--min-region-size', type=int,
#     help='Skip regions below --min-region-size when calculating metrics.')

# parser.add_argument('--max-region-size', type=int,
#     help='Skip regions over --max-region-size when calculating metrics. Helps speed up runs for debugging purposes since there are non-linear increases in runtime with region size.')

# parser.add_argument('--skip-grid-metrics', action='store_true',
#     help='If --slip-grid-metrics, then do not calculate polygon-mesh related metrics and only calculate basic region metrics.')


args = parser.parse_args()



def run_show_channels(ome_tiff_fp, sep):
    channels = utils.get_ome_tiff_channels(ome_tiff_fp)
    print(sep.join(channels))
    

def make_ome(input_tif, output_fp, platform='phenocycler', bbox=None):
    ext = input_tif.split('.')[-1]
    if platform == 'phenocycler' and ext != 'qptiff':
        raise RuntimeError('phenocycler platform option must use .qptiff as input file')
    if platform == 'codex' and ext != 'tif':
        raise RuntimeError('codex platform option must use multichannel imagej .tif as input file')
    if 'ome.tiff' != output_fp[-8:]:
        raise RuntimeError('output filepath must have .ome.tiff extension')

    if platform == 'raw':
        fps = sorted(utils.listfiles(input_tif, regex=r'.tif[f]*$'))
        generate_ome_from_tifs(fps, output_fp, platform=platform, bbox=bbox)
    elif platform == 'codex':
        generate_ome_from_codex_imagej_tif(input_tif, output_fp, bbox=bbox)
    elif platform == 'phenocycler':
        generate_ome_from_qptiff(input_tif, output_fp, bbox=bbox)
    logging.info(f'ome.tiff written to {output_fp}')


def segment_ome(input_tif, output_prefix, split_size, nuclei_markers, membrane_markers):
    logging.info(f'starting segmentation for {input_tif}')
    labeled_cells, labeled_nuclei = segment_cells(
        input_tif, split_size=split_size, nuclei_channels=nuclei_markers,
        membrane_channels=membrane_markers)
    logging.info(f'finished segmentation')
    
    logging.info(f'writing {output_prefix}_nuclei_segmentation.tif')
    tifffile.imwrite(f'{output_prefix}_nuclei_segmentation.tif', labeled_nuclei, compression='LZW')
    logging.info(f'writing {output_prefix}_cell_segmentation.tif')
    tifffile.imwrite(f'{output_prefix}_cell_segmentation.tif', labeled_cells, compression='LZW')


def run_generate_spatial_features(labeled_fp, ome_fp, output_prefix='output', thresholds_fp=None):
    if thresholds_fp is not None:
        df = pd.read_csv(thresholds_fp, sep='\t', header=None)
        thresholds = {k:v for k, v in zip(df.iloc[:, 0], df.iloc[:, 1])}
    else:
        thresholds = None
    df, a = get_spatial_features(labeled_fp, ome_fp, output_prefix=output_prefix, thresholds=thresholds)
    logging.info(f'spatial features written to {output_prefix}.h5ad')
    a.write_h5ad(f'{output_prefix}.h5ad')
    logging.info(f'spatial features written to {output_prefix}.txt')
    df.to_csv(f'{output_prefix}.txt', sep='\t', index=False)


# def run_generate_region_features():
#     df = pd.read_csv(args.spatial_features, sep='\t', index_col=0)
#     metadata_cols = list(df.columns[2:])
#     cols = ['x', 'y']
#     cols += metadata_cols
#     df.columns = cols

#     channel_df = pd.read_csv(args.channel_thresholds_grid, sep='\t')
#     channel_to_thresh_grid = {c:t for c, t in zip(channel_df.iloc[:, 0], channel_df.iloc[:, 1])}

#     channel_df = pd.read_csv(args.channel_thresholds_pixel, sep='\t')
#     channel_to_thresh_pixel = {c:t for c, t in zip(channel_df.iloc[:, 0], channel_df.iloc[:, 1])}
 
#     generate_region_metrics(
#         df, args.ome_tiff, args.regions_mask, args.output_dir,
#         y_col='y', x_col='x', cell_metadata_cols=metadata_cols,
#         boundary_dist=args.boundary_dist,
#         parallel_step=args.parallel_step, perp_steps=args.perp_steps,
#         expansion=args.expansion, grouping_dist=args.breakage_dist,
#         area_thresh=args.area_thresh, group_line_thresh=args.breakage_line_thresh,
#         channel_to_thresh_grid=channel_to_thresh_grid,
#         channel_to_thresh_pixel=channel_to_thresh_pixel,
#         min_region_size=args.min_region_size,
#         max_region_size=args.max_region_size, calculate_grid_metrics=not args.skip_grid_metrics
    # )

def main():
    if args.mode == 'make-ome':
        bbox = [int(x) for x in args.bbox.split(',')] if args.bbox is not None and isinstance(
            args.bbox, str) else None
        make_ome(args.input_tif, args.output_filepath, platform=args.platform, bbox=bbox)
    elif args.mode == 'segment-ome':
        nuclei_markers = args.nuclei_channels.split(',')
        membrane_markers = args.membrane_channels.split(',')
        segment_ome(args.input_tif, args.output_prefix, args.split_size,
                 nuclei_markers=nuclei_markers, membrane_markers=membrane_markers)
    elif args.mode == 'generate-spatial-features':
        run_generate_spatial_features(
            args.labeled_image, args.input_tif, output_prefix=args.output_prefix)
    # elif args.mode == 'generate-region-features':
    #     run_generate_region_features()
    elif args.mode == 'show-channels':
        run_show_channels(args.ome_tiff, args.sep)
    else:
        raise RuntimeError(f'{args.mode} is not a valid mode.')
    


if __name__ == '__main__':
    main()
