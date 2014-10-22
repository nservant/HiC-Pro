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
CONF=$conf_file . $dir/hic.inc.sh
#read_config $ncrna_conf

DATA_DIR=${MAPC_OUTPUT}/data/
DOC_DIR=${MAPC_OUTPUT}/doc
PIC_DIR=${MAPC_OUTPUT}/pic

## picdir
## docdir
## bwtdir
## samplename

sample_dir=$(get_sample_dir $r)
$(get_stat_file)

#for RES_FILE_NAME in ${DATA_DIR}/*
#do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    echo Launching ${R_PATH}/R --no-save CMD BATCH "--args picDir='${PIC_DIR}' docDir='${DOC_DIR}' bwtDir='${BOWTIE2_OUTPUT_DIR}' sampleName='${RES_FILE_NAME}'" ${SCRIPTS}/plotMappingPortion.R ${LOGS_DIR}/plotMappingPortion.Rout

    ${R_PATH}/R --no-save CMD BATCH "--args picDir='${PIC_DIR}' docDir='${DOC_DIR}' bwtDir='${BOWTIE2_OUTPUT_DIR}' sampleName='${RES_FILE_NAME}'" ${SCRIPTS}/plotMappingPortion.R ${LOGS_DIR}/plotMappingPortion.Rout
#done

