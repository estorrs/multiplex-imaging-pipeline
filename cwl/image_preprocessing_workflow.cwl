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
- id: upload_ome
  in:
  - id: dataset
    source: specimen_id
  - id: image_name
    source: specimen_id
  - id: group
    source: group
  - id: project
    source: project
  - id: port
    source: port
  - id: host
    source: host
  - id: filepath
    source: make_ome/output_ome_tiff
  label: upload_ome
  out: []
  run: ../submodules/omero-wrapper/cwl/omero_wrapper_upload.cwl
