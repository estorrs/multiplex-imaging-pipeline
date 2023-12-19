FROM ubuntu:18.04

RUN apt-get update && apt-get install -y vim wget curl git build-essential

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

RUN conda create -n mip -c conda-forge -c defaults python==3.9 -y
ENV PATH="/miniconda/envs/mip/bin:$PATH"

# for now installing requirements with pip up here so debug goes faster, will remove eventually
RUN pip install --no-cache-dir deepcell
RUN pip install numpy pandas einops matplotlib scanpy seaborn tifffile ome-types scikit-image scikit-learn imagecodecs>=2022.7.27 torch torchvision "pydantic<=1.10.12" timm pytorch-lightning

COPY . /multiplex-imaging-pipeline

RUN pip install git+https://github.com/estorrs/wombat

RUN pip install /multiplex-imaging-pipeline[segmentation]

CMD /bin/bash
