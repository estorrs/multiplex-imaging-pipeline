$namespaces:
  sbg: https://www.sevenbridges.com/
arguments:
- position: 1
  valueFrom: generate-region-features
baseCommand:
- mip
class: CommandLineTool
cwlVersion: v1.0
id: region_features
inputs:
- id: input_tif
  inputBinding:
    position: '0'
    prefix: --input-tif
  type: File
- id: spatial_features
  inputBinding:
    position: '0'
    prefix: --spatial-features
  type: File
- id: mask_tif
  inputBinding:
    position: '0'
    prefix: --mask-tif
  type: File?
- default: Pan-Cytokeratin,E-cadherin
  id: mask_markers
  inputBinding:
    position: '0'
    prefix: --mask-markers
  type: string?
- default: output
  id: output_prefix
  inputBinding:
    position: '0'
    prefix: --output-prefix
  type: string?
- default: /miniconda/envs/mip/bin:$PATH
  id: environ_PATH
  type: string?
label: region_features
outputs:
- id: output_txt
  outputBinding:
    glob: '*.txt'
  type: File
- id: region_mask
  outputBinding:
    glob: '*region_mask.tif'
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/multiplex-imaging-pipeline:0.0.1
- class: ResourceRequirement
  ramMin: 100000
- class: EnvVarRequirement
  envDef:
    PATH: $(inputs.environ_PATH)
