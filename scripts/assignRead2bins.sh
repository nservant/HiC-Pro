#!/bin/bash
## Nicolas Servant
## Launcher for assignRead2bins.pl scripts

################### Initialize ###################
set -- $(getopt c:i:g:b:s:h "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf=$2; shift;;
	(-i) map=$2; shift;;
	(-g) genome=$2; shift;;
	(-b) bsize=$2; shift;;
##	(-s) step=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

################### Read the config file ###################

while read curline_read; do
    curline=`echo ${curline_read} | sed -e 's/ = /=/'`

    if [[ $curline != \#* && ! -z $curline ]]; then
	var=`echo $curline | awk -F= '{print $1}'`
	val=`echo $curline | awk -F= '{print $2}'`
	export ${var}="${val}"
    fi
done < $conf

################### Define Variables ###################

RES_FILE_NAME=`basename ${RAW_DIR}`
DATA_DIR=${MAPC_OUTPUT}/data/${RES_FILE_NAME}
MATRIX_DIR=${MAPC_OUTPUT}/matrix/${RES_FILE_NAME}/raw

mkdir -p ${MATRIX_DIR}/${bsize}

################### Combine Bowtie mapping ###################

for i in $GENOME_CHRMS
do
    for j in $GENOME_CHRMS
    do
	if [[ "$i" == "$j" && (${MAPS} == "cis" || ${MAPS} == "all") ]]
	then
	    perl ${SCRIPTS}/assignRead2bins.pl --mapTxt=${DATA_DIR}/${RES_FILE_NAME}.${ORGANISM}.interaction --genomeDesc=$genome --binSize=${bsize} --step=${BIN_STEP} ${BIN_OPTS} --chrchr=chr${i}_${ORGANISM}.chr${j}_${ORGANISM} --output=${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}
	fi
	if [[ "$i" != "$j" && (${MAPS} == "trans" || ${MAPS} == "all") ]]
	then
	    perl ${SCRIPTS}/assignRead2bins.pl --mapTxt=${DATA_DIR}/${RES_FILE_NAME}.${ORGANISM}.interaction --genomeDesc=$genome --binSize=${bsize} --step=${BIN_STEP} ${BIN_OPTS} --chrchr=chr${i}_${ORGANISM}.chr${j}_${ORGANISM} --output=${MATRIX_DIR}/${bsize}/${RES_FILE_NAME}_${bsize}
	fi
    done
done
