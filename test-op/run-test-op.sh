#!/bin/bash

export PATH=$(pwd)/../bin/:$PATH
##export PATH=/bioinfo/local/curie/ngs-data-analysis/centos/HiC-Pro_2.10.1/bin/:$PATH

function die
{
    echo -e "Unit testing failed !"
    exit 1
}

## HiC-Pro should be install and in the PATH
which HiC-Pro > /dev/null;
if [ $? != "0" ]; then
    echo -e "Can not proceed without HiC-Pro, please install and/or set the PATH variable"
    exit 1;
fi
version=$(HiC-Pro -v | awk '{print $NF}')


RES_PREFIX=HiC-Pro_testop_${version}

RUN_STANDALONE_MODE=1
RUN_CLUSTER_MODE=0
RUN_AS_MODE=0
RUN_CAP_MODE=0

if [[ $RUN_STANDALONE_MODE == 1 ]]; then
    CONFIG=config_test_latest.txt

    ## Get data
    if [[ $(ls bowtie2_indexes_hg19 | wc -l) == 0 ]]; then
	echo -e "bowtie2_indexes_hg19 folder looks empty. Please provide the hg19 bowtie2 indexes or download them using "
	echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip"
	exit 1
    fi

    if [[ ! -e test_data ]]; then
	wget https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata.tar.gz && tar -zxvf HiCPro_testdata.tar.gz
	/bin/rm -f HiCPro_testdata.tar.gz
    fi

    ## standalone complete pipeline
    cmd="time HiC-Pro -i test_data/ -o ${RES_PREFIX}_all -c ${CONFIG}"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi

    ## standalone stepwise pipeline
    cmd="time HiC-Pro -i test_data -o ${RES_PREFIX} -c ${CONFIG} -s mapping -s quality_checks"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi
    
    cmd="time HiC-Pro -i ${RES_PREFIX}/bowtie_results/bwt2 -o ${RES_PREFIX} -c ${CONFIG} -s proc_hic -s quality_checks"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi
    
    cmd="time HiC-Pro -i ${RES_PREFIX}/hic_results/data/ -o ${RES_PREFIX} -c ${CONFIG} -s merge_persample -s build_contact_maps -s ice_norm"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi
    
    cmd="time HiC-Pro -i ${RES_PREFIX}/hic_results/data/ -o ${RES_PREFIX} -c ${CONFIG} -s build_contact_maps"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi
    
    cmd="time HiC-Pro -i  ${RES_PREFIX}/hic_results/matrix/ -o ${RES_PREFIX} -c ${CONFIG} -s ice_norm"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi
fi

## Cluster complete pipeline
if [[ $RUN_CLUSTER_MODE == 1 ]]; then
    CONFIG=config_test_latest.txt
    time HiC-Pro -i test_data/ -o ${RES_PREFIX}_cluster -c ${CONFIG} -p
    file2sub=$(find ${RES_PREFIX}_cluster -name "*step1*")
    qsub ${file2sub}
fi

## as test
if [[ $RUN_AS_MODE == 1 ]]; then
    CONFIG=config_test_as.txt
    
    if [[ ! -e test_data_as || ! -e bowtie2_indexes_as_test ]]; then
	wget https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata_as.tar.gz && tar -zxvf HiCPro_testdata_as.tar.gz
	/bin/rm -f HiCPro_testdata_as.tar.gz
    fi
    export CONFIG=config_test_as.txt
    cmd="time HiC-Pro -i test_data_as -o ${RES_PREFIX}_as -c ${CONFIG}"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi
fi

if [[ $RUN_CAP_MODE == 1 ]]; then
    CONFIG=config_test_cap.txt

    ## Get data
    if [[ $(ls bowtie2_indexes_mm9 | wc -l) == 0 ]]; then
        echo -e "bowtie2_indexes_mm9 folder looks empty. Please provide the mm9 bowtie2 indexes or download them using "
        echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm9.zip"
        exit 1
    fi

    if [[ ! -e test_data_cap ]]; then
        wget https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata_cap.tar.gz && tar -zxvf HiCPro_testdata_cap.tar.gz
        /bin/rm -f HiCPro_testdata_cap.tar.gz
    fi
    cmd="time HiC-Pro -i test_data_cap -o ${RES_PREFIX}_cap -c ${CONFIG}"
    echo $cmd
    eval $cmd
    if [[ $? != 0 ]]; then die; fi
fi



#/bin/rm -rf HiC_Pro_test-op*
