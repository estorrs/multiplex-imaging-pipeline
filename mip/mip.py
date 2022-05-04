import argparse
import logging

import numpy as np
import pandas as pd

import mip.utils as utils
from mip.ome import generate_ome_from_tifs
from mip.spatial_features import save_spatial_features

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('mode', type=str,
    choices=['make-ome', 'generate-spatial-features'],
    help='Which task mip is to execute.')

##############
## make-ome ##
##############
parser.add_argument('--tif-directory', type=str,
    help='Used in make-ome mode. Directory of stitched tif files that will be combined into a single ome.tiff file.')

parser.add_argument('--output-filepath', type=str,
    help='Location to write ome.tiff file')

parser.add_argument('--platform', type=str,
    choices=['codex'], default='codex',
    help='Which platform produced the input images.')

###########################
## generate-spatial-features ##
###########################
parser.add_argument('--label-image', type=str,
    help='Filepath of labeled cell segmentation tif.')

parser.add_argument('--ome-tiff', type=str,
    help='Filepath of ome.tiff')

parser.add_argument('--output-prefix', type=str, default='output',
    help='Filepath of ome.tiff')

args = parser.parse_args()


def run_make_ome(fps, output_fp, platform='codex'):
    generate_ome_from_tifs(fps, output_fp, platform=platform)
    logging.info(f'ome.tiff written to {output_fp}')

def run_generate_spatial_features(label_fp, ome_fp, output_prefix='output'):
    save_spatial_features(label_fp, ome_fp, output_prefix=output_prefix)
    logging.info(f'spatial features written to {output_prefix}')


def main():
    if args.mode == 'make-ome':
        fps = sorted(utils.listfiles(args.tif_directory, regex=r'.tif[f]*$'))
        run_make_ome(fps, args.output_filepath, platform=args.platform)
    elif args.mode == 'generate-spatial-features':
        run_generate_spatial_features(
            args.label_image, args.ome_tiff, output_prefix=args.output_prefix)
    else:
        raise RuntimeError(f'{args.mode} is not a valid mode.')
 

if __name__ == '__main__':
    main()
