#!/bin/bash
## Nicolas Servant
## Eric Viara
## assignRead2bins.sh
## Launcher for assignRead2bins.pl script

dir=$(dirname $0)

. $dir/hic.inc.sh

################### Initialize ###################

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

################### Define Variables ###################

DATA_DIR=${MAPC_OUTPUT}/data/

################### Combine Bowtie mapping ###################

for RES_FILE_NAME in ${DATA_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    if [ -d ${DATA_DIR}/${RES_FILE_NAME} ]; then
	echo RES: $RES_FILE_NAME $bsize
	MATRIX_DIR=${MAPC_OUTPUT}/matrix/${RES_FILE_NAME}/raw
	for bsize in ${BIN_SIZE}
	do
	    mkdir -p ${MATRIX_DIR}/${bsize}
	    cat ${DATA_DIR}/${RES_FILE_NAME}/*.interaction | ${SCRIPTS}/build_matrix --binsize ${bsize} --chrsizes $genome --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize} &
	done
    fi
    wait
done
