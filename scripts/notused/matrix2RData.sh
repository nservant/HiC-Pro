#!/bin/bash
## Nicolas Servant & Eric Viara
## Launcher for matrix2RData.R script

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

RDATA_DIR=${MAPC_OUTPUT}/rdata
DATA_DIR=${MAPC_OUTPUT}/data/

for RES_FILE_NAME in ${DATA_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    RES_FILE_NAME_OBJ=$(echo $RES_FILE_NAME | sed -e 's/-/_/g')
    MATRIX_DIR=${MAPC_OUTPUT}/matrix/${RES_FILE_NAME}/raw
    for bsize in ${BIN_SIZE}
    do
	mkdir -p ${MATRIX_DIR}/${BSIZE}
	echo "Launching" ${R_PATH}/R --no-save CMD BATCH "--args matDir='${MATRIX_DIR}/${bsize}' obj='${RES_FILE_NAME_OBJ}_${bsize}' cpu='${N_CPU}' rdataDir='${RDATA_DIR}' org='${ORGANISM}'" ${SCRIPTS}/matrix2RData.R ${LOGS_DIR}/matrix2RData.Rout
	${R_PATH}/R --no-save CMD BATCH "--args matDir='${MATRIX_DIR}/${bsize}' obj='${RES_FILE_NAME_OBJ}_${bsize}' cpu='${N_CPU}' rdataDir='${RDATA_DIR}/${RES_FILE_NAME}' org='${ORGANISM}'" ${SCRIPTS}/matrix2RData.R ${LOGS_DIR}/matrix2RData.Rout
    done
done

