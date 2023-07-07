$namespaces:
  sbg: https://www.sevenbridges.com/
arguments:
- position: 1
  valueFrom: generate-spatial-features
baseCommand:
- mip
class: CommandLineTool
cwlVersion: v1.0
id: spatial_features
inputs:
- id: input_tif
  inputBinding:
    position: '0'
    prefix: --input-tif
  type: File
- id: labeled_image
  inputBinding:
    position: '0'
    prefix: --labeled-image
  type: File
- id: thresholds
  inputBinding:
    position: '0'
    prefix: --thresholds
  type: File?
- default: output
  id: output_prefix
  inputBinding:
    position: '0'
    prefix: --output-prefix
  type: string?
- default: /miniconda/envs/mip/bin:$PATH
  id: environ_PATH
  type: string?
label: spatial_features
outputs:
- id: output_txt
  outputBinding:
    glob: '*.txt'
  type: File
- id: output_h5ad
  outputBinding:
    glob: '*.h5ad'
  type: File
- id: cell_type_image
  outputBinding:
    glob: '*.png'
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/multiplex-imaging-pipeline:0.0.1
- class: ResourceRequirement
  ramMin: 60000
- class: EnvVarRequirement
  envDef:
    PATH: $(inputs.environ_PATH)
