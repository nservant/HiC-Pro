FROM ubuntu:latest
LABEL authors="Nicolas Servant" \
      description="Docker image containing all requirements for the HiC-Pro pipeline"

## Install gcc for pip iced install
RUN apt-get update && apt-get install -y gcc g++ && apt-get clean -y

COPY conda.yml /
RUN conda env create -f /conda.yml && conda clean -a
ENV PATH /opt/conda/envs/HiC-Pro_v3.0.0/bin:$PATH
