BootStrap: docker
From: ubuntu:latest

%labels
    AUTHOR Nicolas Servant

%pre
    apt-get install -y debootstrap

%files
    environment.yml /opt/environment.yml

%post
    apt-get update
    apt-get install -y wget
    apt-get install -y gzip
    apt-get install -y bzip2
    apt-get install -y curl
    apt-get install -y unzip

    ## g++
    apt-get install -y build-essential
    
    # install anaconda
    if [ ! -d /usr/local/anaconda ]; then
       wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh \
       	    -O ~/anaconda.sh && \
	    bash ~/anaconda.sh -b -p /usr/local/anaconda && \
	    rm ~/anaconda.sh
    fi

    # set anaconda path
    export PATH=$PATH:/usr/local/anaconda/bin
    conda update conda

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    
    # Let us save some space
    conda clean --packages -y

    # external tools
    echo "Installing conda environment ... "
    conda env create -f /opt/environment.yml
    echo "source activate HiC-Pro_v3.1.0" > ~/.bashrc

    # Install HiC-pro
    echo "Installing latest HiC-Pro release ..."
    VERSION=$(curl -s https://github.com/nservant/HiC-Pro/releases/latest | egrep -o 'v3.[0-9]*.[0-9]*')
    echo "v"$VERSION".zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && unzip hicpro_latest.zip
    #VERSION="devel_py3"
    #echo $VERSION".zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && unzip hicpro_latest.zip
    
    cd $(echo HiC-Pro-$VERSION)
    make configure
    make install
 
    # Let us save some space
    conda clean --packages -y
    conda clean --all -y
    rm -rf /usr/local/anaconda/pkgs

%test
    INSTALLED_HICPRO_VERSION=$(find /usr/local/bin -name HiC-Pro | xargs dirname)
    $INSTALLED_HICPRO_VERSION/HiC-Pro -h

%environment
    export PATH=$PATH:/usr/local/anaconda/bin
    INSTALLED_VERSION=$(find /usr/local/bin -name HiC-Pro | xargs dirname)
    export PATH=$PATH:$INSTALLED_VERSION 
    export LANG=C
