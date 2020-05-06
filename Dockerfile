FROM ubuntu:latest
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

LABEL authors="Nicolas Servant" \
      description="Docker image containing all requirements for the HiC-Pro pipeline"


# Install miniconda. Copied from https://hub.docker.com/r/continuumio/miniconda/~/dockerfile/ (because we need debian:stretch to install gcc 6 which is needed for HiCPro)
RUN apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.3.27-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

## Install gcc
RUN apt-get update && apt-get install -y gcc && apt-get clean

ENV PATH /opt/conda/bin:$PATH

## Install all dependencies using conda
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/HiC-Pro_v3.0.0/bin:$PATH

# Install HiCPro
RUN cd /tmp && \
    echo "devel_py3.zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && \
    unzip hicpro_latest.zip && \
    cd HiC-Pro-devel_py3  && \ 
    make configure prefix=/ && \
    make install && \
    cd .. && \
    rm -fr HiC-Pro* && \
    ln -s /HiC-Pro_3.0.0 /HiC-Pro