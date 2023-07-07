class: Workflow
cwlVersion: v1.0
id: full_imaging_workflow
inputs:
- id: specimen_id
  type: string
- id: input_tif
  type: File
- default: phenocycler
  id: platform
  type: string?
- id: bbox
  type: string?
- default: DAPI
  id: nuclei_channels
  type: string?
- default: Pan-Cytokeratin,E-cadherin,CD45,CD8,CD3e,CD4,Vimentin,SMA,CD31
  id: membrane_channels
  type: string?
- default: Pan-Cytokeratin,E-cadherin
  id: mask_markers
  type: string?
- default: Multiplex_Imaging
  id: project
  type: string?
- default: htan-imaging.wucon.wustl.edu
  id: host
  type: string?
- default: '4064'
  id: port
  type: string?
- default: HTAN
  id: group
  type: string?
label: full_imaging_workflow
outputs:
- id: ome_tiff
  outputSource: image_preprocessing/ome_tiff
  type: File
- id: labeled_nuclei
  outputSource: image_processing_to_thresholds/labeled_nuclei
  type: File
- id: labeled_cells
  outputSource: image_processing_to_thresholds/labeled_cells
  type: File
- id: spatial_features_txt
  outputSource: image_processing_to_thresholds/spatial_features_txt
  type: File
- id: spatial_features_h5ad
  outputSource: image_processing_to_thresholds/spatial_features_h5ad
  type: File
- id: region_features_txt
  outputSource: image_processing_to_thresholds/region_features_txt
  type: File
- id: labeled_regions
  outputSource: image_processing_to_thresholds/labeled_regions
  type: File
requirements: []
steps:
- id: image_preprocessing
  in:
  - id: specimen_id
    source: specimen_id
  - id: input_tif
    source: input_tif
  - id: platform
    source: platform
  - id: bbox
    source: bbox
  - id: group
    source: group
  - id: project
    source: project
  - id: port
    source: port
  - id: host
    source: host
  label: image_preprocessing
  out:
  - id: ome_tiff
  run: image_preprocessing_workflow.cwl
- id: image_processing_to_thresholds
  in:
  - id: specimen_id
    source: specimen_id
  - id: ome_tiff
    source: image_preprocessing/ome_tiff
  - id: nuclei_channels
    source: nuclei_channels
  - id: membrane_channels
    source: membrane_channels
  - id: mask_markers
    source: mask_markers
  - id: group
    source: group
  - id: project
    source: project
  - id: port
    source: port
  - id: host
    source: host
  label: image_processing_to_thresholds
  out:
  - id: labeled_nuclei
  - id: labeled_cells
  - id: spatial_features_txt
  - id: spatial_features_h5ad
  - id: region_features_txt
  - id: labeled_regions
  run: image_processing_no_thresholds_workflow.cwl
