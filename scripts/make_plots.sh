#!/bin/bash
## HiC-Pro
## Copyleft 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Launcher for all plotting function in R
##

dir=$(dirname $0)


################### Initialize ###################
#set -- $(getopt c:i:g:b:s:h "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

################### Read the config file ###################

#read_config $ncrna_conf
CONF=$conf_file . $dir/hic.inc.sh

################### Define Variables ###################

DATA_DIR=${MAPC_OUTPUT}/data/

################### Combine Bowtie mapping ###################

for RES_FILE_NAME in ${DATA_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    PIC_DIR=${MAPC_OUTPUT}/pic/${RES_FILE_NAME}
    DATA_DIR=${MAPC_OUTPUT}/data/${RES_FILE_NAME}

    ## Check if output directory exists
    if [ ! -d ${PIC_DIR} ]; then mkdir -p ${PIC_DIR}; fi

    ## make plots
    echo "Plot mapping results ..."
    cmd="${R_PATH}/R --no-save CMD BATCH \"--args picDir='${PIC_DIR}' bwtDir='${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}' sampleName='${RES_FILE_NAME}' r1tag='${PAIR1_EXT}' r2tag='${PAIR2_EXT}'\" ${SCRIPTS}/plot_mapping_portion.R ${LOGS_DIR}/plot_mapping_portion.Rout"
    echo $cmd
    eval $cmd
    
    echo "Plot pairing results ..."
    cmd="${R_PATH}/R --no-save CMD BATCH \"--args picDir='${PIC_DIR}' bwtDir='${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}' sampleName='${RES_FILE_NAME}' rmMulti='${RM_MULTI}' rmSingle='${RM_SINGLETON}'\" ${SCRIPTS}/plot_pairing_portion.R ${LOGS_DIR}/plot_pairing_portion.Rout"
    echo $cmd
    eval $cmd
    
    echo "Plot Hi-C processing ..."
    cmd="${R_PATH}/R --no-save CMD BATCH \"--args picDir='${PIC_DIR}' hicDir='${DATA_DIR}' sampleName='${RES_FILE_NAME}'\" ${SCRIPTS}/plot_hic_fragment.R ${LOGS_DIR}/plot_hic_fragment.Rout"
    echo $cmd
    eval $cmd
done
