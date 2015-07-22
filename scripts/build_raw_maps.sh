#!/bin/bash
## HiC-Pro
## Copyleft 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## assignRead2bins.sh
## Launcher for assignRead2bins.pl script
##

dir=$(dirname $0)

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

GENOME_SIZE_FILE=`abspath $GENOME_SIZE`
if [[ ! -e $GENOME_SIZE_FILE ]]; then
    GENOME_SIZE_FILE=$ANNOT_DIR/$GENOME_SIZE
    if [[ ! -e $GENOME_SIZE_FILE ]]; then
	echo "$GENOME_SIZE not found. Exit"
	exit -1
    fi
fi

################### Create Matrix files ###################

for RES_FILE_NAME in ${DATA_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    ## Logs
    LDIR=${LOGS_DIR}/${RES_FILE_NAME}
    mkdir -p ${LDIR}

    if [ -d ${DATA_DIR}/${RES_FILE_NAME} ]; then
	MATRIX_DIR=${MAPC_OUTPUT}/matrix/${RES_FILE_NAME}/raw
	for bsize in ${BIN_SIZE}
	do
	    mkdir -p ${MATRIX_DIR}/${bsize}

	    ## Build haplotype contact maps if specified
	    if [[ ! -z ${ALLELE_SPECIFIC_SNP} ]]; then
		#awk -F"\t" '$9~/1-0|0-1|1-1/{print}' ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs | ${SCRIPTS}/build_matrix --matrix-format ${MATRIX_FORMAT} --binsize ${bsize} --chrsizes $ANNOT_DIR/$GENOME_SIZE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_G1 2> ${LDIR}/build_raw_maps_G1.log
		#awk -F"\t" '$9~/2-0|0-2|2-2/{print}' ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs | ${SCRIPTS}/build_matrix --matrix-format ${MATRIX_FORMAT} --binsize ${bsize} --chrsizes $ANNOT_DIR/$GENOME_SIZE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_G2 2> ${LDIR}/build_raw_maps_G2.log
	    	cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs_G1 | ${SCRIPTS}/build_matrix --matrix-format ${MATRIX_FORMAT} --binsize ${bsize} --chrsizes $ANNOT_DIR/$GENOME_SIZE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_G1 2> ${LDIR}/build_raw_maps_G1.log
		cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs_G2 | ${SCRIPTS}/build_matrix --matrix-format ${MATRIX_FORMAT} --binsize ${bsize} --chrsizes $ANNOT_DIR/$GENOME_SIZE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_G2 2> ${LDIR}/build_raw_maps_G2.log

	    else
		## Build normal diploid maps
		cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs | ${SCRIPTS}/build_matrix --matrix-format ${MATRIX_FORMAT} --binsize ${bsize} --chrsizes $ANNOT_DIR/$GENOME_SIZE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize} 2> ${LDIR}/build_raw_maps.log
	    fi
	done
    fi
    wait
done
