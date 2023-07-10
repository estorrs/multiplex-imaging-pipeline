$namespaces:
  sbg: https://www.sevenbridges.com/
arguments:
- position: 1
  valueFrom: make-ome
baseCommand:
- mip
class: CommandLineTool
cwlVersion: v1.0
id: ome_generation
inputs:
- id: input_tif
  inputBinding:
    position: '0'
    prefix: --input-tif
  type: File
- default: output.ome.tiff
  id: output_filepath
  inputBinding:
    position: '0'
    prefix: --output-filepath
  type: string?
- default: phenocycler
  id: platform
  inputBinding:
    position: '0'
    prefix: --platform
  type: string?
- id: bbox
  inputBinding:
    position: '0'
    prefix: --bbox
  type: string?
- default: /miniconda/envs/mip/bin:$PATH
  id: environ_PATH
  type: string?
label: ome_generation
outputs:
- id: output_ome_tiff
  outputBinding:
    glob: '*.ome.tiff'
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/multiplex-imaging-pipeline:0.0.1
- class: ResourceRequirement
  ramMin: 200000
- class: EnvVarRequirement
  envDef:
    PATH: $(inputs.environ_PATH)
