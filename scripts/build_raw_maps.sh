#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## build_raw_maps.sh
## Launcher for build_matrix C++ code
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

################### Define Input Directory ###################                                                                                                                                                 
input_data_type=$(get_data_type)

if [[ $input_data_type == "allvalid" ]]
then
    DATA_DIR=${RAW_DIR}
else
    DATA_DIR=${MAPC_OUTPUT}/data/
fi

nbf=$(find -L ${DATA_DIR} -mindepth 1 | wc -l)
if [[ $nbf == 0 ]]; then die "Error : empty ${DATA_DIR} folder."; fi

################### Define Variables ###################

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
	    if [[ ${bsize} == -1 ]]; then
		GENOME_FRAGMENT_FILE=`abspath $GENOME_FRAGMENT`
		if [[ $GENOME_FRAGMENT == "" || ! -f $GENOME_FRAGMENT_FILE ]]; then
		    GENOME_FRAGMENT_FILE=$ANNOT_DIR/$GENOME_FRAGMENT
		    if [[ ! -f $GENOME_FRAGMENT_FILE ]]; then
			echo "GENOME_FRAGMENT not found. Cannot generate fragment level Hi-C map"
			exit 1
		    fi
		fi
		bsize="rfbin"
	    fi
	    
	    mkdir -p ${MATRIX_DIR}/${bsize}

	    ## Build haplotype contact maps if specified
	    if [[ ! -z ${ALLELE_SPECIFIC_SNP} ]]; then
	    	cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs_G1 | ${SCRIPTS}/build_matrix --matrix-format ${MATRIX_FORMAT} --binsize ${bsize} --chrsizes $GENOME_SIZE_FILE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_G1 2> ${LDIR}/build_raw_maps_G1.log
		cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs_G2 | ${SCRIPTS}/build_matrix --matrix-format ${MATRIX_FORMAT} --binsize ${bsize} --chrsizes $GENOME_SIZE_FILE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}_G2 2> ${LDIR}/build_raw_maps_G2.log

	    else
		## Build normal diploid maps
		## RS resolution
		if [[ ${bsize} == "rfbin" ]]; then
		    cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs | ${SCRIPTS}/build_matrix --binfile ${GENOME_FRAGMENT_FILE} --chrsizes $GENOME_SIZE_FILE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize} --matrix-format ${MATRIX_FORMAT} 2> ${LDIR}/build_raw_maps.log
		## Bin resolution
		else
		    #echo " cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs | ${SCRIPTS}/build_matrix --binsize ${bsize} --chrsizes $GENOME_SIZE_FILE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize} --matrix-format ${MATRIX_FORMAT} 2> ${LDIR}/build_raw_maps.log"
		    cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs | ${SCRIPTS}/build_matrix --binsize ${bsize} --chrsizes $GENOME_SIZE_FILE --ifile /dev/stdin --oprefix ${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize} --matrix-format ${MATRIX_FORMAT} 2> ${LDIR}/build_raw_maps.log
		fi
	    fi
	done
    fi
    wait
done
