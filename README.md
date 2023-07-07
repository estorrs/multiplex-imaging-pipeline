# multiplex-imaging-pipeline

A pipeline for multiplex imaging analysis


## Installation

```bash
pip install multiplex-imaging-pipeline
```

## Usage

### ome.tiff generation

```bash
mip make-ome --input-tif INPUT_TIF --platform PLATFORM --bbox BBOX --output-filepath OUTPUT_FILEPATH
```

###### Arguments
+ --input-tif
  + Input tif or directory of tifs to be converted into HTAN compatible ome.tiff format. Default is "phenocycler". If --platform is "phenocycler", then a file with a .qptiff extension output by the phenocycler platform is expected. If --platform is "codex", then a multichannel imagej .tif output by the first-generation codex machine is expected. If --platform is "raw" then a directory of tifs is expected, where the images will be combined into a multichannel .ome.tiff and the channel names will be named based on the files (for example if one of the files is named protein1.tif, then that channel in the ome.tiff will be named protein1).

+ --platform
  + Can be one of ["phenocycler", "codex", "raw"]. For Akoya images, use "phenocycler" for their new platform (.qptiff file extension), and use "codex" for the first-generation machine (multichannel imagej .tif file). To convert a directory of single channel .tif files into a multichannel ome.tiff, use "raw". See --input-tif for more details.

+ --bbox
  + Crops image to given pixel coordinates. Useful for phenocycler images that image multiple pieces of tissue on the same slide. --bbox is specified as a string and is in the format "TOP,BOTTOM,LEFT,RIGHT".

+ --output-filepath
  + Where to write output .ome.tiff. Default is output.ome.tiff.
 
###### Examples

**Making an ome.tiff from a phenocycler .qptiff file**

```bash
mip make-ome --input-tif /path/to/file.qptiff --platform phenocycler --output-filepath output.ome.tiff
```

**Making an ome.tiff from a phenocycler .qptiff file where the output image is cropped to the given pixel coordinates ("TOP,BOTTOM,LEFT,RIGHT").**

```bash
mip make-ome --input-tif /path/to/file.qptiff --platform phenocycler --output-filepath output.ome.tiff --bbox "0,10000,2000,12000" 
```

**Making an ome.tiff from a codex multichannel imagej .tif file**

```bash
mip make-ome --input-tif /path/to/file.tif --platform codex --output-filepath output.ome.tiff
```

**Making an ome.tiff from a directory of .tif files**

```bash
mip make-ome --input-tif /path/to/dir/of/tifs/ --platform raw --output-filepath output.ome.tiff
```

### segmentation of ome.tiff file

```bash
mip segment-ome --input-tif INPUT_TIF --output-prefix OUTPUT_PREFIX --split-size SPLIT_SIZE --nuclei-channels NUCLEI_CHANNELS --membrane-channels MEMBRANE_CHANNELS
```

###### Arguments

+ --input-tif
  + Path to ome.tiff to be segmented.

+ --output-prefix
  + Output prefix where labeled cell and nuclei segmentation files will be written. --output-prefix can include directory paths. Default is "output". Two files will be written by segment-ome: 1) {OUTPUT_PREFIX}_labeled_cells.tif and 2) {OUTPUT_PREFIX}_labeled_nuclei.tif. In the labeled .tif files, pixels belong to a cell will have that cells integer ID, background pixels have a value of zero.

+ --split-size
  + If image is larger than --split-size, image will be split into tiles for segmentation and then recombined afterwards. Default is 25000. If you run into memory issues, consider decreasing this value. 

+ --nuclei-channels
  + Which channels in --input-tif to use for nuclei segmentation. If multiple channels are provided, then those channels are merged into a single image that is used for segmentation. Channels are provided as a string in the following format: "channel1,channel2,.....,channelz". Default is "DAPI". Note: if channel is not in image, than it is attempted to be converted to a more standard name using the CHANNEL_MAPPING [here](https://github.com/estorrs/multiplex-imaging-pipeline/blob/main/multiplex_imaging_pipeline/utils.py#L19).

+ --membrane-channels
  + Which channels in --input-tif to use for membrane segmentation. If multiple channels are provided, then those channels are merged into a single image that is used for segmentation. Channels are provided as a string in the following format: "channel1,channel2,.....,channelz". Default is "Pan-Cytokeratin,E-cadherin,CD45,CD8,CD3e,Vimentin,SMA,CD31,C20". Note: if channel is not in image, than it is attempted to be converted to a more standard name using the CHANNEL_MAPPING [here](https://github.com/estorrs/multiplex-imaging-pipeline/blob/main/multiplex_imaging_pipeline/utils.py#L19).

###### Examples

**Segmentation of ome.tiff using default parameters**

```bash
mip segment-ome --input-tif /path/to/file.ome.tiff --output-prefix output
```

**Segmentation of ome.tiff using custom nuclei and membrane markers**

```bash
mip segment-ome --input-tif /path/to/file.ome.tiff --output-prefix output --nuclei-markers "DAPI" --membrane-markers "Pan-Cytokeratin,E-cadherin,Vimentin".
```



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

