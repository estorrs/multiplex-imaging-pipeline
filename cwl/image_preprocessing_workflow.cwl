class: Workflow
cwlVersion: v1.0
id: image_preprocessing_workflow
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
label: image_preprocessing_workflow
outputs:
- id: ome_tiff
  outputSource: make_ome/output_ome_tiff
  type: File
requirements: []
steps:
- id: make_ome
  in:
  - id: platform
    source: platform
  - id: bbox
    source: bbox
  - id: input_tif
    source: input_tif
  label: make_ome
  out:
  - id: output_ome_tiff
  run: ome_generation.cwl
