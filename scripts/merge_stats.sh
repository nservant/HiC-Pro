#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## Merge multiple stat files
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

##read_config $ncrna_conf
CONF=$conf_file . $dir/hic.inc.sh

DATA_DIR=${MAPC_OUTPUT}/data/

################### Combine Bowtie mapping ###################
for RES_FILE_NAME in ${RAW_DIR}/*/
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    ldir=${LOGS_DIR}/${RES_FILE_NAME}
    echo "Logs: ${ldir}/merge_stats.log"

    odir=${MAPC_OUTPUT}/stats/${RES_FILE_NAME}
    mkdir -p ${odir}
    
    ##mapstat
    if [[ -d ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME} ]]; then
	nb_map_r1=$(find -L ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.bam" -and -name "*${PAIR1_EXT}*" | wc -l)
	nb_map_r2=$(find -L ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.bam" -and -name "*${PAIR2_EXT}*" | wc -l)

	if [[ $nb_map_r1 -gt 0 && $nb_map_r2 -gt 0 ]]; then
	    echo "Merge mapstat files ..." >> ${ldir}/merge_stats.log
	    cmd="${PYTHON_PATH}/python ${SCRIPTS}/merge_statfiles.py -d ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}/ -p "*${PAIR1_EXT}*.mapstat" -v > ${odir}/${RES_FILE_NAME}${PAIR1_EXT}.mmapstat"
	    exec_cmd $cmd > ${ldir}/merge_stats.log 2>&1
	    cmd="${PYTHON_PATH}/python ${SCRIPTS}/merge_statfiles.py -d ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}/ -p "*${PAIR2_EXT}*.mapstat" -v > ${odir}/${RES_FILE_NAME}${PAIR2_EXT}.mmapstat"
	     exec_cmd $cmd >> ${ldir}/merge_stats.log 2>&1
	fi  
    fi

    ##pairstat
    if [[ -d ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME} ]]; then
	nb_pairs=$(find -L ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.pairstat" | wc -l)
	if [[ $nb_pairs -gt 0 ]]; then
	    echo "Merge pairstat files ..." >> ${ldir}/merge_stats.log
	    cmd="${PYTHON_PATH}/python ${SCRIPTS}/merge_statfiles.py -d ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}/ -p \"*.pairstat\" -v> ${odir}/${RES_FILE_NAME}.mpairstat"
	    exec_cmd $cmd >> ${ldir}/merge_stats.log 2>&1
	fi
    fi

   ##RSstat 
    if [[ -d ${DATA_DIR}/${RES_FILE_NAME} ]]; then
	nb_rs=$(find -L ${DATA_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.RSstat" | wc -l)
	if [[ $nb_rs -gt 0 ]]; then
	    echo "Merge RSstat files ..." >> ${ldir}/merge_stats.log
            cmd="${PYTHON_PATH}/python ${SCRIPTS}/merge_statfiles.py -d ${DATA_DIR}/${RES_FILE_NAME}/ -p \"*.RSstat\" -v> ${odir}/${RES_FILE_NAME}.mRSstat"
	    exec_cmd $cmd >> ${ldir}/merge_stats.log 2>&1
	fi
    fi
done


