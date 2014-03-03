#!/bin/bash
## Nicolas Servant
##
NORMAL="\\033[0;39m" 
RED="\\033[1;31m"
BLUE="\\033[0;34m"


################### Initialize ###################
mode='global'
samplelist=""
set -- $(getopt c:p:l "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-l) mode='local'; shift;;
	(-p) samplelist=$2; shift;;
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
BOWTIE2_IDX=${BOWTIE2_IDX_PATH}/${ORGANISM}
BOWTIE2_LOCAL_OUTPUT_DIR=${BOWTIE2_OUTPUT_DIR}/bwt2_local/${RES_FILE_NAME}
BOWTIE2_GLOBAL_OUTPUT_DIR=${BOWTIE2_OUTPUT_DIR}/bwt2_global/${RES_FILE_NAME}
UNMAP_READ_DIR=${BOWTIE2_OUTPUT_DIR}/unmap/${RES_FILE_NAME}

################### Run Bowtie on all fastq files on batch mode ###################
echo "SL=$samplelist"
if [[ ${samplelist} == "" ]]; then
    echo "BATCH MODE"
## Global mapping
    if [[ ${mode} == 'global' ]]; then
	for r in ${RAW_DIR}/*.fastq
	do
	    echo ${r} >> ${LOGFILE}
	    prefix=`basename ${r} | sed -e 's/.fastq//'`
	    
	## Align reads
 ##--un ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}_${ORGANISM}_bunmap.fastq
	    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${r} -S ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}_${ORGANISM}.bwt2glob.sam 2>>${LOGS_DIR}/bowtie_${prefix}_global_${ORGANISM}.log"
	    echo $cmd
	    eval $cmd
	    
	## Generate BAM files
	    cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}_${ORGANISM}.bwt2glob.sam > ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}_${ORGANISM}.bwt2glob.bam"
	    echo $cmd
	    eval $cmd
	done
    else
## Local mapping
	for r in ${BOWTIE2_GLOBAL_OUTPUT_DIR}/*.fastq
	do
	    echo ${r} >> ${LOGFILE}
	    prefix=`basename ${r} | sed -e 's/.fastq//'`
	    
	## Align reads
	    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${r} -S ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}_bwt2loc.sam 2>>${LOGS_DIR}/bowtie_${prefix}_local.log"
	    echo $cmd
	    eval $cmd
	    
	## Generate BAM files
	    cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}_bwt2loc.sam > ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}_bwt2loc.bam"
	    echo $cmd
	    eval $cmd
	done
    fi

################### Run Bowtie on all fastq files on PBS mode ###################

else
    echo "PBS MODE"
    cd $PBS_O_WORKDIR
## Global mapping
    if [[ ${mode} == 'global' ]]; then
	
    ## Make cvs file
	r=$(grep "^${PBS_ARRAYID};" ${samplelist} | cut -d';' -f2)
	echo ${r} >> ${LOGFILE}
	prefix=`basename ${r} | sed -e 's/.fastq//'`
	
    ## Align reads
	cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${r} -S ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}_${ORGANISM}.bwt2glob.sam 2>>${LOGS_DIR}/bowtie_${prefix}_global_${ORGANISM}.log"
	echo $cmd
	eval $cmd
	
    ## Generate BAM files
	cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}_${ORGANISM}.bwt2glob.sam > ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}_${ORGANISM}.bwt2glob.bam"
	echo $cmd
	eval $cmd
	
    else
## Local mapping
	
	echo ${r} >> ${LOGFILE}
	prefix=`basename ${r} | sed -e 's/.fastq//'`
	
	## Align reads
	cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${r} -S ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}_bwt2loc.sam 2>>${LOGS_DIR}/bowtie_${prefix}_local.log"
	echo $cmd
	eval $cmd
	
	## Generate BAM files
	cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}_bwt2loc.sam > ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}_bwt2loc.bam"
	echo $cmd
	eval $cmd
    fi
fi
