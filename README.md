# multiplex-imaging-pipeline

A pipeline for multiplex imaging analysis.

Under active development


## Installation

```bash
pip install multiplex-imaging-pipeline
```

## Usage

### ome.tiff generation

ome.tiff files can be created with the `make-ome` command.

```bash
mip make-ome --input-tif INPUT_TIF --platform PLATFORM --bbox BBOX --output-filepath OUTPUT_FILEPATH
```

#### Arguments
+ --input-tif
  + Input tif or directory of tifs to be converted into HTAN compatible ome.tiff format. If --platform is "phenocycler", then a file with a .qptiff extension output by the phenocycler platform is expected. If --platform is "codex", then a multichannel imagej .tif output by the first-generation codex machine is expected. If --platform is "raw" then a directory of tifs is expected, where the images will be combined into a multichannel .ome.tiff and the channel names will be named based on the files (for example if one of the files is named protein1.tif, then that channel in the ome.tiff will be named protein1).

+ --platform
  + Can be one of ["phenocycler", "codex", "raw"]. For Akoya images, use "phenocycler" for their new platform (.qptiff file extension), and use "codex" for the first-generation machine (multichannel imagej .tif file). To convert a directory of single channel .tif files into a multichannel ome.tiff, use "raw". See --input-tif for more details. Default is "phenocycler".

+ --bbox
  + Crops image to given pixel coordinates. Useful for phenocycler images that image multiple pieces of tissue on the same slide. --bbox is specified as a string and is in the format "TOP,BOTTOM,LEFT,RIGHT".

+ --output-filepath
  + Where to write output .ome.tiff. Default is output.ome.tiff.
 
#### Examples

###### Making an ome.tiff from a phenocycler .qptiff file

```bash
mip make-ome --input-tif /path/to/file.qptiff --platform phenocycler --output-filepath output.ome.tiff
```

###### Making an ome.tiff from a phenocycler .qptiff file where the output image is cropped to the given pixel coordinates ("TOP,BOTTOM,LEFT,RIGHT").

```bash
mip make-ome --input-tif /path/to/file.qptiff --platform phenocycler --output-filepath output.ome.tiff --bbox "0,10000,2000,12000" 
```

###### Making an ome.tiff from a codex multichannel imagej .tif file

```bash
mip make-ome --input-tif /path/to/file.tif --platform codex --output-filepath output.ome.tiff
```

###### Making an ome.tiff from a directory of .tif files

```bash
mip make-ome --input-tif /path/to/dir/of/tifs/ --platform raw --output-filepath output.ome.tiff
```

### segmentation of ome.tiff file

```bash
mip segment-ome --input-tif INPUT_TIF --output-prefix OUTPUT_PREFIX --split-size SPLIT_SIZE --nuclei-channels NUCLEI_CHANNELS --membrane-channels MEMBRANE_CHANNELS
```

#### Arguments

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

#### Examples

###### Segmentation of ome.tiff using default parameters

```bash
mip segment-ome --input-tif /path/to/file.ome.tiff --output-prefix output
```

###### Segmentation of ome.tiff using custom nuclei and membrane markers

```bash
mip segment-ome --input-tif /path/to/file.ome.tiff --output-prefix output --nuclei-markers "DAPI" --membrane-markers "Pan-Cytokeratin,E-cadherin,Vimentin".
```

### generation of spatial features

```bash
mip generate-spatial-features --input-tif INPUT_TIF --output-prefix OUTPUT_PREFIX --labeled-image LABELED_IMAGE --gating-strategy GATING_STRATEGY --thresholds THRESHOLDS
```

#### Arguments

+ --input-tif
  + Path to ome.tiff to be quantified. Should match --labeled-image segmentation result.

+ --output-prefix
  + Output prefix where feature tables and cell type image will be written. --output-prefix can include directory paths. Default is "output". Three files will be written: 1) {OUTPUT_PREFIX}_spatial_features.h5ad and 2) {OUTPUT_PREFIX}_spatial_features.txt, and 3) {OUPUT_PREFIX}_annotated_cell_types.png. spatial_features.* files contain feature tables describing cell morphology, marker intensities (raw and normalized), positive pixel fractions (if --thresholds is specified), and cell type (based on --gating-strategy). annotated_cell_types.png is an image of cell boundaries colored by cell type annotation.

+ --labeled-image
  + Labeled .tif file containing cell segmentation information used when generating spatial features. In the labeled .tif files, pixels belong to a cell will have that cells integer ID, background pixels have a value of zero.
 
+ --gating-strategy
  + A .yaml file specifying the gating strategy to use when identifying cell types. By default, the gating strategy [here](https://github.com/estorrs/multiplex-imaging-pipeline/blob/main/multiplex_imaging_pipeline/spatial_features.py#L15) will be used.
 
+ --thresholds
  + If provided, the given manually defined thresholds will be used to calculate "positive pixel fraction" for specified markers, which is the % of positive pixels for that marker in a given cell. Useful when cell typing and wanting to be certain about eliminating batch effects. By default thresholds are automatically calculated. --thresholds is a tab-seperated .txt file where the first column is a channel name, and the second is the value to use as a threshold (i.e. pixels with an intensity above this threshold will be considered positive). No header should be included in the file.
 
#### Examples

###### Spatial feature generation using default parameters

```bash
mip generate-spatial-features --input-tif /path/to/file.ome.tiff --labeled-image /path/to/segmentation/labeled/image.tif --output-prefix output
```

###### Spatial feature generation using custom gating strategy

```bash
mip generate-spatial-features --input-tif /path/to/file.ome.tiff --labeled-image /path/to/segmentation/labeled/image.tif --output-prefix output --gating-strategy /path/to/file.yaml
```

###### Spatial feature generation using thresholds file

```bash
mip generate-spatial-features --input-tif /path/to/file.ome.tiff --labeled-image /path/to/segmentation/labeled/image.tif --output-prefix output --thresholds /path/to/file.txt
```

### generation of region features

```bash
mip generate-region-features --input-tif INPUT_TIF --output-prefix OUTPUT_PREFIX --spatial-features SPATIAL_FEATURES --mask-tif MASK_TIF --mask-markers MASK_MARKERS
```

#### Arguments

+ --input-tif
  + Path to ome.tiff to be segmented.

+ --output-prefix
  + Output prefix where feature tables and cell type image will be written. --output-prefix can include directory paths. Default is "output". Two types of files will be written: 1) `{OUTPUT_PREFIX}_region_features.txt` and 2) `{OUTPUT_PREFIX}_{NAME}_mask.tif`. `region_features.txt` is a tab-seperated file containing feature describing each region in the image. There are three main features: `{NAME}_cell_type_fraction_*` (fraction of each cell type in region type with given NAME), `{NAME}_marker_intensity_*` (raw intensity in region type with given NAME), `{NAME}_marker_intensity_scaled_*` (intensity scaled by std in region type with given NAME), `{NAME}_marker_fraction_*` (positive pixel fraction in region type with given NAME, thresholds are whatever was used in generate-spatial-features). `*.tif` files are labeled images delineating each region type.

+ --spatial-features
  + `.h5ad` object output by generate-spatial-features
 
+ --mask-tif
  + Boolean tif defining regions to use. If none is provided, channels in --mask-markers will be used to auto-generate masks.
 
+ --mask-markers
  + when --mask-tif is not provided, will be used to auto-generate masks. Defaults to 'Pan-Cytokeratin,E-cadherin'.
  
#### Examples

###### Region feature generation using default parameters

```bash
mip generate-region-features --input-tif /path/to/file.ome.tiff --spatial-features /path/to/spatial/features.h5ad --output-prefix output
```

###### Region feature generation using user provided mask

```bash
mip generate-region-features --input-tif /path/to/file.ome.tiff --spatial-features /path/to/spatial/features.h5ad --output-prefix output --mask-tif /path/to/mask.tif
```

###### Region feature generation using mask generated from custom markers

```bash
mip generate-region-features --input-tif /path/to/file.ome.tiff --spatial-features /path/to/spatial/features.h5ad --output-prefix output --mask-markers "CD4,CD45,CD8"
```


## Docker

```bash
docker pull estorrs/multiplex-imaging-pipeline:0.0.1
```

Example of ome.tiff generation with docker.

Note that you'll need to map input and output data folders with the -v flag.

```bash
docker run -v /path/to/input/dir:/inputs -v /path/to/output/dir:/outputs estorrs/multiplex-imaging-pipeline:0.0.1 mip make-ome --input-tif /inputs/file.qptiff --platform phenocycler --output-filepath /outputs/output.ome.tiff
```


