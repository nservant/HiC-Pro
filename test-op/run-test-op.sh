#!/bin/bash

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

## Get data
if [[ $(ls bowtie2_indexes | wc -l) == 0 ]]; then
    echo -e "bowtie2_indexes folder looks empty. Please provide the hg19 bowtie2 indexes or download them using "
    echo -e "wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip"
    exit 1
fi

#wget https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata.tar.gz && tar -zxvf HiCPro_testdata.tar.gz
#/bin/rm -f HiCPro_testdata.tar.gz

RES_PREFIX=HiC_Pro_test-op
CONFIG=config_test_latest.txt
RUN_CLUSTER_MODE=0

## standalone complete pipeline

cmd="time HiC-Pro -i test_data/ -o ${RES_PREFIX} -c ${CONFIG}"
echo $cmd
eval $cmd
if [[ $? != 0 ]]; then die; fi

## standalone stepwise pipeline

cmd="time HiC-Pro -i test_data -o ${RES_PREFIX}_stepwise1 -c ${CONFIG} -s mapping -s quality_checks"
echo $cmd
eval $cmd
if [[ $? != 0 ]]; then die; fi

cmd="time HiC-Pro -i ${RES_PREFIX}_stepwise1/bowtie_results/bwt2 -o ${RES_PREFIX}_stepwise2 -c ${CONFIG} -s proc_hic -s quality_checks"
echo $cmd
eval $cmd
if [[ $? != 0 ]]; then die; fi

cmd="time HiC-Pro -i ${RES_PREFIX}_stepwise2/hic_results/data/ -o ${RES_PREFIX}_stepwise3 -c ${CONFIG} -s merge_persample -s build_contact_maps -s ice_norm"
echo $cmd
eval $cmd
if [[ $? != 0 ]]; then die; fi

cmd="time HiC-Pro -i ${RES_PREFIX}_stepwise3/hic_results/data/ -o ${RES_PREFIX}_stepwise4 -c ${CONFIG} -s build_contact_maps"
echo $cmd
eval $cmd
if [[ $? != 0 ]]; then die; fi

cmd="time HiC-Pro -i  ${RES_PREFIX}_stepwise4/hic_results/matrix/ -o ${RES_PREFIX}_stepwise5 -c ${CONFIG} -s ice_norm"
echo $cmd
eval $cmd
if [[ $? != 0 ]]; then die; fi

## Cluster complete pipeline
if [[ $RUN_CLUSTER_MODE == 1 ]]; then
    time HiC-Pro -i test_data/ -o ${RES_PREFIX}_cluster -c ${CONFIG} -p
    file2sub=$(find ${RES_PREFIX}_cluster -name "*step1*")
    qsub ${file2sub}
fi

## as test
#export CONFIG=config_test_as.txt

#time HiC-Pro -i test_data_as -o ${RES_PREFIX}_5 -c config_test_as.txt

## dnase test
#export CONFIG=config_test_dnase.txt

#time HiC-Pro -i test_data/ -o ${RES_PREFIX}_6 -c config_test_dnase.txt


#/bin/rm -rf HiC_Pro_test-op*
