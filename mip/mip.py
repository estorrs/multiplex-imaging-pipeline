import argparse
import logging

import numpy as np
import pandas as pd

import mip.utils as utils
import mip.ome as ome

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('mode', type=str,
    choices=['make-ome'],
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

args = parser.parse_args()


def run_make_ome(fps, output_fp, platform='codex'):
    ome.generate_ome_from_tifs(fps, output_fp, platform=platform)
    logging.info(f'ome.tiff written to {output_fp}')


def main():
    if args.mode == 'make-ome':
        fps = sorted(utils.listfiles(args.tif_directory, regex=r'.tif[f]*$'))
        run_make_ome(fps, args.output_filepath, platform=args.platform)
    else:
        raise RuntimeError(f'{args.mode} is not a valid mode.')
 

if __name__ == '__main__':
    main()
