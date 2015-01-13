#!/bin/bash
## Eric Viara
## plotMappingPortion.sh
## Launcher for plotMappingPortiont.R script

dir=$(dirname $0)

##. $dir/hic.inc.sh

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
    ${R_PATH}/R --no-save CMD BATCH "--args picDir='${PIC_DIR}' bwtDir='${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}' sampleName='${RES_FILE_NAME}' r1tag='${PAIR1_EXT}' r2tag='${PAIR2_EXT}'" ${SCRIPTS}/plotMappingPortion.R ${LOGS_DIR}/plotMappingPortion.Rout
    
    echo "Plot pairing results ..."
    ${R_PATH}/R --no-save CMD BATCH "--args picDir='${PIC_DIR}' bwtDir='${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}' sampleName='${RES_FILE_NAME}' rmMulti='${RM_MULTI}' rmSingle='${RM_SINGLETON}'" ${SCRIPTS}/plotPairingPortion.R ${LOGS_DIR}/plotPairingPortion.Rout
    
    echo "Plot Hi-C processing ..."
    ${R_PATH}/R --no-save CMD BATCH "--args picDir='${PIC_DIR}' hicDir='${DATA_DIR}' sampleName='${RES_FILE_NAME}'" ${SCRIPTS}/plotHiCFragment.R ${LOGS_DIR}/plotHiCFragment.Rout  
    wait
done
