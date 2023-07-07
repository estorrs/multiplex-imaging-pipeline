import argparse
import logging

import numpy as np
import pandas as pd
import tifffile
import yaml

import multiplex_imaging_pipeline.utils as utils
from multiplex_imaging_pipeline.ome import generate_ome_from_tifs, generate_ome_from_qptiff, generate_ome_from_codex_imagej_tif
from multiplex_imaging_pipeline.spatial_features import get_spatial_features
from multiplex_imaging_pipeline.region_features import get_region_features
from multiplex_imaging_pipeline.segmentation import segment_cells

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('mode', type=str,
    choices=['make-ome', 'segment-ome', 'generate-spatial-features', 'generate-region-features', 'show-channels'],
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

parser.add_argument('--gating-strategy', type=str,
    help='Filepath to .yaml file containing gating strategy to use while annotating cells.')

###############################
## generate-region-features ##
###############################
# parser.add_argument('--input-tif', type=str,
#     help='ome.tiff file to generate spatial features for')

# parser.add_argument('--output-prefix', type=str, default='output',
#     help='Output prefix to use when writing output files. Two files will be written: *_spatial_features.h5ad and *_spatial_features.txt. Default is "output". For example, if --output-prefix is "path/to/out/directory/sample" then the two output files will be named path/to/out/directory/sample_spatial_features.h5ad and path/to/out/directory/sample_spatial_features.txt.')

parser.add_argument('--mask-tif', type=str,
    help='Filepath of region mask tif.')

parser.add_argument('--mask-markers', type=str, default='Pan-Cytokeratin,E-cadherin',
    help='If --mask-tif is not provided, these markers will be used to generate a region mask.')

parser.add_argument('--spatial-features', type=str,
    help='Filepath to .h5ad file output by generate-spatial-features.')


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


def run_generate_spatial_features(labeled_fp, ome_fp, output_prefix='output',
                                  thresholds_fp=None, gating_strategy_fp=None):
    if thresholds_fp is not None:
        df = pd.read_csv(thresholds_fp, sep='\t', header=None)
        thresholds = {k:v for k, v in zip(df.iloc[:, 0], df.iloc[:, 1])}
    else:
        thresholds = None

    if gating_strategy_fp is not None:
        gating_strategy = yaml.safeload(open(gating_strategy_fp))
    else:
        gating_strategy = None
    
    df, a = get_spatial_features(labeled_fp, ome_fp, output_prefix=output_prefix,
                                 thresholds=thresholds, gating_strategy=gating_strategy)
    logging.info(f'spatial features written to {output_prefix}_spatial_features.h5ad')
    a.write_h5ad(f'{output_prefix}_spatial_features.h5ad')
    logging.info(f'spatial features written to {output_prefix}_spatial_features.txt')
    df.to_csv(f'{output_prefix}_spatial_features.txt', sep='\t', index=False)

    logging.info(f'saving annotated cell segmentation image to {output_prefix}_annotated_cell_types.png')
    fp = f'{output_prefix}_annotated_cell_types.png'
    utils.save_annotated_rgba(tifffile.imread(labeled_fp), a, fp, figsize=(25, 25), dpi=300)


def run_generate_region_features():
    combined, labeled_dict = get_region_features(
        args.spatial_features, args.input_tif, mask_fp=args.mask_tif,
        mask_markers=args.mask_markers.split(','))

    logging.info(f'writing region features table to {args.output_prefix}_region_features.txt')
    combined.to_csv(f'{args.output_prefix}_region_features.txt', sep='\t', index=False)

    logging.info(f'writing region masks to {args.output_prefix}_*.tif')
    for name, img in labeled_dict.items():
        tifffile.imwrite(f'{args.output_prefix}_{name}_mask.tif', img, compression='LZW')


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
            args.labeled_image, args.input_tif, output_prefix=args.output_prefix,
            gating_strategy_fp=args.gating_strategy)
    elif args.mode == 'generate-region-features':
        run_generate_region_features()
    elif args.mode == 'show-channels':
        run_show_channels(args.input_tif, args.sep)
    else:
        raise RuntimeError(f'{args.mode} is not a valid mode.')
    


if __name__ == '__main__':
    main()
