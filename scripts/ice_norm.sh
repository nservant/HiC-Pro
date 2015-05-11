#!/bin/bash
## HiC-Pro
## Copyleft 2015 Institut Curie                               
## Author(s): Nicolas Servant, Nelle Varoquaux
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Launcher for ICE normalization scripts
##

dir=$(dirname $0)

#. $dir/hic.inc.sh

################### Initialize ###################

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
    ## Logs
    LDIR=${LOGS_DIR}/${RES_FILE_NAME}
    mkdir -p ${LDIR}

    if [ -d ${DATA_DIR}/${RES_FILE_NAME} ]; then
	NORM_DIR=${MAPC_OUTPUT}/matrix/${RES_FILE_NAME}/iced
	for bsize in ${BIN_SIZE}
	do
	    mkdir -p ${NORM_DIR}/${bsize}
	    INPUT_MATRIX=${MAPC_OUTPUT}/matrix/${RES_FILE_NAME}/raw/${bsize}/${RES_FILE_NAME}_${bsize}.matrix
	    ln -f -s ../../raw/${bsize}/${RES_FILE_NAME}_${bsize}_abs.bed ${NORM_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_abs.bed
	    ln -f -s ../../raw/${bsize}/${RES_FILE_NAME}_${bsize}_ord.bed  ${NORM_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_ord.bed
	    ${SCRIPTS}/ice --results_filename ${NORM_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_iced.matrix --filtering_perc ${SPARSE_FILTERING} --max_iter ${MAX_ITER} --eps ${EPS} --verbose 1 ${INPUT_MATRIX} > ${LDIR}/ice.log

	done
    fi
    wait
done
