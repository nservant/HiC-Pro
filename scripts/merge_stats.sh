#!/bin/bash
## HiC-Pro
## Copyleft 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

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
for RES_FILE_NAME in ${BOWTIE2_FINAL_OUTPUT_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    LDIR=${LOGS_DIR}/${RES_FILE_NAME}

    ##mapstat
    nb_map=$(find -L ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.mapstat" | wc -l)
    if [[ $nb_map -gt 0 ]]; then
	echo "Merge mapstat files ..." >> ${LDIR}/merge_stat.log
	${PYTHON_PATH}/python ${SCRIPTS}/merge_statfiles.py -d ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}/ -p "*.mapstat" > ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.mapstat 2>> ${LDIR}/merge_stat.log
    fi

    ##pairstat
    nb_pairs=$(find -L ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.pairstat" | wc -l)
    if [[ $nb_pairs -gt 0 ]]; then
	echo "Merge pairstat files ..." >> ${LDIR}/merge_stat.log
	${PYTHON_PATH}/python ${SCRIPTS}/merge_statfiles.py -d ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}/ -p"*.pairstat" > ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.pairstat 2>> ${LDIR}/merge_stat.log
    fi

   ##RSstat                                                                                                                                                            
    nb_rs=$(find -L ${DATA_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.RSstat" | wc -l)
    if [[ $nb_rs -gt 0 ]]; then
	echo "Merge RSstat files ..." >> ${LDIR}/merge_stat.log
        ${PYTHON_PATH}/python ${SCRIPTS}/merge_statfiles.py -d ${DATA_DIR}/${RES_FILE_NAME}/ -p"*.RSstat" > ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.RSstat 2>> ${LDIR}/merge_stat.log
    fi


