$namespaces:
  sbg: https://www.sevenbridges.com/
arguments:
- position: 1
  valueFrom: segment-ome
baseCommand:
- python
- /multiplex-imaging-pipeline/mip/mip.py
class: CommandLineTool
cwlVersion: v1.0
id: mesmer_segmentation
inputs:
- id: input_tif
  inputBinding:
    position: '0'
    prefix: --input-tif
  type: File
- default: output
  id: output_prefix
  inputBinding:
    position: '0'
    prefix: --output-prefix
  type: string?
- default: DAPI
  id: nuclei_channels
  inputBinding:
    position: '0'
    prefix: --nuclei-channels
  type: string?
- default: Pan-Cytokeratin,E-cadherin,CD45,CD8,CD3e,Vimentin,SMA,CD31,C20
  id: membrane_channels
  inputBinding:
    position: '0'
    prefix: --membrane-channels
  type: string?
- default: 25000
  id: split_size
  inputBinding:
    position: '0'
    prefix: --split-size
  type: int?
- default: /miniconda/envs/mip/bin:$PATH
  id: environ_PATH
  type: string?
label: mesmer_segmentation
outputs:
- id: labeled_cells
  outputBinding:
    glob: '*_cell_segmentation.tif'
  type: File
- id: labeled_nuclei
  outputBinding:
    glob: '*_nuclei_segmentation.tif'
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/multiplex-imaging-pipeline:0.0.1
- class: ResourceRequirement
  ramMin: 100000
- class: EnvVarRequirement
  envDef:
    PATH: $(inputs.environ_PATH)
