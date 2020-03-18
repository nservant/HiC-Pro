BootStrap: docker
From: ubuntu:latest

%labels
    AUTHOR Nicolas Servant

%pre
    apt-get install -y debootstrap

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
    echo "Installing dependancies ... "
    conda install -y bowtie2
    conda install -y samtools

    # Python (>2.7) with *pysam (>=0.8.3)*, *bx(>=0.5.0)*, *numpy(>=1.8.2)*, and *scipy(>=0.15.1)* libraries
    conda install -y -c conda-forge python=2.7.15
    conda install -y -c anaconda scipy=1.2.1 
    conda install -y -c anaconda numpy=1.16.3
    conda install -y -c bcbio bx-python=0.8.2
    conda install -y -c bioconda pysam=0.15.2 

    # Install R
    conda update readline
    #conda install -c conda-forge readline=6.2
    conda install -c r r-base 
    conda install -c r r-ggplot2=2.2.1
    conda install -c r r-rcolorbrewer
    conda install -c r r-gridbase	

    # Install MultiQC
    conda install -c bioconda multiqc=1.7 
   
    # Install HiC-pro
    echo "Installing latest HiC-Pro release ..."
    #VERSION=$(curl -s https://github.com/nservant/HiC-Pro/releases/latest | egrep -o '2.[0-9]*.[0-9]*')
    #echo "v"$VERSION".zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && unzip hicpro_latest.zip
    VERSION="devel"
    echo $VERSION".zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && unzip hicpro_latest.zip
    
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
