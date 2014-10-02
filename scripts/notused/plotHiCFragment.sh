#!/bin/bash
## Eric Viara
## plotHiCFragment.sh
## Launcher for plotHiCFragment.R script

dir=$(dirname $0)

. $dir/hic.inc.sh

################### Initialize ###################
#set -- $(getopt c:i:g:b:s:h "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

################### Read the config file ###################

read_config $ncrna_conf

DATA_DIR=${MAPC_OUTPUT}/data/
PIC_DIR=${MAPC_OUTPUT}/pic

for RES_FILE_NAME in ${DATA_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    echo ${R_PATH}/R --no-save CMD BATCH "--args picDir='${PIC_DIR}' dataDir='${DATA_DIR}/${RES_FILE_NAME}' bwtDir='${BOWTIE2_OUTPUT_DIR}' sampleName='${RES_FILE_NAME}'" ${SCRIPTS}/plotHiCFragment.R ${LOGS_DIR}/plotHiCFragment.Rout
    ${R_PATH}/R --no-save CMD BATCH "--args picDir='${PIC_DIR}' dataDir='${DATA_DIR}/${RES_FILE_NAME}' bwtDir='${BOWTIE2_OUTPUT_DIR}' sampleName='${RES_FILE_NAME}'" ${SCRIPTS}/plotHiCFragment.R ${LOGS_DIR}/plotHiCFragment.Rout
done

