class: Workflow
cwlVersion: v1.0
id: image_processing_no_thresholds_workflow
inputs:
- id: specimen_id
  type: string
- id: ome_tiff
  type: File
- id: mask_tif
  type: File?
- default: DAPI
  id: nuclei_channels
  type: string?
- default: Pan-Cytokeratin,E-cadherin,CD45,CD8,CD3e,CD4,Vimentin,SMA,CD31
  id: membrane_channels
  type: string?
- default: Pan-Cytokeratin,E-cadherin
  id: mask_markers
  type: string?
label: image_processing_no_thresholds_workflow
outputs:
- id: labeled_nuclei
  outputSource: segment_ome/labeled_nuclei
  type: File
- id: labeled_cells
  outputSource: segment_ome/labeled_cells
  type: File
- id: spatial_features_txt
  outputSource: generate_spatial_features/output_txt
  type: File
- id: spatial_features_h5ad
  outputSource: generate_spatial_features/output_h5ad
  type: File
- id: cell_type_image
  outputSource: generate_spatial_features/cell_type_image
  type: File
- id: region_features_txt
  outputSource: generate_region_features/output_txt
  type: File
- id: labeled_regions
  outputSource: generate_region_features/region_mask
  type: File
requirements: []
steps:
- id: segment_ome
  in:
  - id: input_tif
    source: ome_tiff
  - id: nuclei_channels
    source: nuclei_channels
  - id: membrane_channels
    source: membrane_channels
  label: segment_ome
  out:
  - id: labeled_cells
  - id: labeled_nuclei
  run: segmentation.cwl
- id: generate_spatial_features
  in:
  - id: input_tif
    source: ome_tiff
  - id: labeled_image
    source: segment_ome/labeled_cells
  label: generate_spatial_features
  out:
  - id: output_txt
  - id: output_h5ad
  - id: cell_type_image
  run: spatial_features.cwl
- id: generate_region_features
  in:
  - id: input_tif
    source: ome_tiff
  - id: spatial_features
    source: generate_spatial_features/output_h5ad
  - id: mask_markers
    source: mask_markers
  label: generate_region_features
  out:
  - id: output_txt
  - id: region_mask
  run: region_features.cwl
