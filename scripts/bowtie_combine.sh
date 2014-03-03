#!/bin/bash
## Nicolas Servant
##


NORMAL="\\033[0;39m" 
RED="\\033[1;31m"
BLUE="\\033[0;34m"


################### Initialize ###################
set -- $(getopt c:l:g:o:h "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-l) BOWTIE2_LOCAL_OUTPUT_DIR=$2; shift;;
	(-g) BOWTIE2_GLOBAL_OUTPUT_DIR=$2; shift;;
	(-o) BOWTIE2_FINAL_OUTPUT_DIR=$2; shift;;
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
done < $ncrna_conf

################### Define Variables ###################

RES_FILE_NAME=`basename ${RAW_DIR}`

################### Combine Bowtie mapping ###################

## Merge Local and Global alignment
for r in ${BOWTIE2_GLOBAL_OUTPUT_DIR}/*.aln
do
    prefix=`basename ${r} | sed -e 's/.bwt2glob.aln//'`
    echo ${prefix} >> ${LOGFILE}

    perl ${SCRIPTS}/mergeBwt2GlobLoc.pl -a ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.aln -b ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.aln -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}_bowtie_final.aln;
done

if [ -e ${BOWTIE2_OUTPUT_DIR}/${RES_FILE_NAME}_R1.final.aln ];then
    /bin/rm -f ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}_R1.final.aln
fi
if [ -e ${BOWTIE2_OUTPUT_DIR}/${RES_FILE_NAME}_R2.final.aln ];then
    /bin/rm -f ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}_R2.final.aln
fi

## Merge aln files from multiple samples/input files
for r in ${BOWTIE2_FINAL_OUTPUT_DIR}/*R1*bowtie_final.aln
do
    cat $r >> ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}_R1.final.aln
done

for r in ${BOWTIE2_FINAL_OUTPUT_DIR}/*R2*bowtie_final.aln
do
    cat $r >> ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}_R2.final.aln
done

${SCRIPTS}/mappingstat.sh -i ${BOWTIE2_FINAL_OUTPUT_DIR} -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${RES_FILE_NAME}
