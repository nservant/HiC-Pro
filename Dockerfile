FROM ubuntu:16.04
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

LABEL authors="Nicolas Servant" \
      description="Docker image containing all requirements for the HiC-Pro pipeline"

## Install system tools
RUN apt-get update \
  && apt-get install -y build-essential \
  wget \
  unzip \
  bzip2 \
  gcc \
  g++ && apt-get clean


## Install miniconda.
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh -O ~/anaconda.sh
RUN bash ~/anaconda.sh -b -p /usr/local/anaconda
RUN rm ~/anaconda.sh
ENV PATH /usr/local/anaconda/bin:$PATH


## Install all dependencies using conda
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/HiC-Pro_v3.0.0/bin:$PATH

## Install HiCPro
RUN cd /tmp && \
    echo "devel_py3.zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && \
    unzip hicpro_latest.zip && \
    cd HiC-Pro-devel_py3  && \ 
    make configure prefix=/ && \
    make install && \
    cd .. && \
    rm -fr HiC-Pro*

RUN /HiC-Pro_3.0.0/bin/HiC-Pro -h