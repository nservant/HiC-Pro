BootStrap: docker
From: ubuntu:latest

%labels
    AUTHOR Nicolas Servant

%files
    environment.yml /opt/environment.yml

%environment
    export LANG=en_US.utf-8
    export INSTALLED_VERSION=$(find /usr/local/bin -name HiC-Pro | xargs dirname)
    export PATH=$PATH:${INSTALLED_VERSION}:${INSTALLED_VERSION}/utils/
    export PATH=/usr/local/conda/bin:$PATH
    source /opt/etc/bashrc

%post
    apt-get update && \
    apt-get install -y \
               wget \
               debootstrap \
               gzip \
               bzip2 \
               curl \
               unzip \
               g++ \

    # use bash as /bin/sh instead of dash (otherwise source /opt/etc/bashrc will not work)
    echo "dash dash/sh boolean false" | debconf-set-selections
    DEBIAN_FRONTEND=noninteractive dpkg-reconfigure dash

    # install miniconda
    if [ ! -d /usr/local/miniconda ]; then
       wget https://repo.continuum.io/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh -O miniconda.sh && \
       sh miniconda.sh -b -p /usr/local/conda && \
       rm miniconda.sh && \
       ln -s /usr/local/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
       find /usr/local/conda/ -follow -type f -name '*.a' -delete && \
       find /usr/local/conda/ -follow -type f -name '*.js.map' -delete && \
       /usr/local/conda/bin/conda clean -afy
    fi

    # set miniconda path
    export PATH=/usr/local/conda/bin:$PATH
    conda update conda
    conda clean --packages -y

    #bashrc
    echo "#! /bin/bash\n\n" > ~/.bashrc
    conda init bash
    
    # external tools
    echo "Installing conda environment ... "
    conda env create -f /opt/environment.yml -n hicpro
    echo "conda activate hicpro" >> ~/.bashrc

    # Install HiC-pro
    export PATH=/usr/local/conda/envs/hicpro/bin:$PATH
    echo "Installing latest HiC-Pro release ..."
    #VERSION=$(curl -s https://github.com/nservant/HiC-Pro/releases/latest | egrep -o 'v3.[0-9]*.[0-9]*')
    VERSION="3.1.0"
    echo "v"$VERSION".zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && unzip hicpro_latest.zip
    
    cd $(echo HiC-Pro-$VERSION)
    make configure
    make install
 
    # Let us save some space
    conda clean --packages -y
    conda clean --all -y
    rm -rf /usr/local/miniconda/pkgs
    apt-get clean

    mkdir /opt/etc/
    cp ~/.bashrc /opt/etc/bashrc

%test
    INSTALLED_HICPRO_VERSION=$(find /usr/local/bin -name HiC-Pro | xargs dirname)
    ${INSTALLED_HICPRO_VERSION}/HiC-Pro -h

