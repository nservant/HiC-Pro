#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## Init
dir=$(dirname $0)
MODE='global'

## Get args
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(-l) MODE='local'; shift;;
	(-p) FASTQLIST=$2; shift;;
	(-u) UNMAP=1; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

## Read configuration files
CONF=$conf_file . $dir/hic.inc.sh

## Bowtie2 wrapper
## Global Alignment
global_align()
{
    local sample_dir="$1"
    local file="$2"
    local prefix="$3"
    local unmap="$4"
    local filtunmap="$5"
    local cmd
    echo ${file} >> ${LOGFILE}

    echo "Running bowtie on full length reads ..."

    ## Check mapping options
    if [[ -z ${BOWTIE2_GLOBAL_OPTIONS} ]]; then
	echo "Mapping step1 options not defined. Exit"
	exit -1
    fi

    ## Output
    mkdir -p ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}
    
    ## Logs
    LDIR=${LOGS_DIR}/${sample_dir}
    mkdir -p ${LDIR}

    ## Unmapped reads
    if [[ $unmap == 1 ]]; then
	BOWTIE2_GLOBAL_OPTIONS=${BOWTIE2_GLOBAL_OPTIONS}" --un ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
    fi

    ## Check for gz files
    if [[ ! -e $file &&  -e $file.gz ]]; then
	    file=$file.gz
    fi

    ## Run bowtie
    if [[ $N_CPU -lt 2 ]]; then
	echo -e "Warning : HiC-Pro need at least 2 CPUs to run the mapping !!"
	bwt_cpu=1
    else
	bwt_cpu=$(( $N_CPU / 2 ))
    fi
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --rg-id BMG --rg SM:${prefix} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${file} 2> ${LDIR}/bowtie_${prefix}_global_${REFERENCE_GENOME}.log "
    if [[ $filtunmap == 1 ]]; then
	cmd=$cmd"| ${SAMTOOLS_PATH}/samtools view -F 4 -bS - > ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.bam"
    else
	cmd=$cmd"> ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.bam"
    fi

    exec_cmd $cmd
}

## Trimmed reads Alignment
local_align()
{
    local sample_dir="$1"
    local file="$2"
    local prefix="$3"
    local unmap="$4"
    local cmd

    echo "Running bowtie on trimmed reads ..."

    ## Check mapping options
    if [[ -z ${BOWTIE2_LOCAL_OPTIONS} ]]; then
	echo "Mapping step2 options not defined. Exit"
	exit -1
    fi

    mkdir -p ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}
    
    ## Logs
    LDIR=${LOGS_DIR}/${sample_dir}
    mkdir -p ${LDIR}

    ## Starts trimming reads from the ligation site
    tfile=`basename $file | sed -e s/.fastq$/_trimmed.fastq/`
    if [[ ${RM_LOCAL_NO_CUTSITE} == 1 ]]; then
	${SCRIPTS}/cutsite_trimming --fastq $file --cutsite ${LIGATION_SITE} --out ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/$tfile --rmuntrim > ${LDIR}/readsTrimming.log 2>&1
    else
	${SCRIPTS}/cutsite_trimming --fastq $file --cutsite ${LIGATION_SITE} --out ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/$tfile  > ${LDIR}/readsTrimming.log 2>&1
    fi

    ## Unmapped reads
    if [[ $unmap == 1 ]]; then
	BOWTIE2_LOCAL_OPTIONS=${BOWTIE2_LOCAL_OPTIONS}" --un ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
    fi

    ## Run bowtie
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg SM:${prefix} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${tfile} 2>${LDIR}/bowtie_${prefix}_local.log | ${SAMTOOLS_PATH}/samtools view -bS - > ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.bam"
    exec_cmd "$cmd"
}


## If the step 2 is not run - do not filter the BAM
if [[ -z $LIGATION_SITE ]]; then
    FILT_UNMAP=0
else
    FILT_UNMAP=1
fi

if [[ ${MODE} == 'global' ]]; then
    for r in $(get_fastq_for_bowtie_global)
    do
	wasrun=1
       	R1=$r
	R2=$(echo $r | get_R2)
	
	echo "-----------"
	echo $R1
	echo $R2
	echo "-------------"

	#if [[ ! -e $R1 || ! -e $R2 ]]; then
	#    echo "error - input files not found." >&2
	#    exit 1
	#fi
	
	sample_dir=$(get_sample_dir $r)
	prefix1=$(basename ${R1} | sed -e 's/.fastq\(.gz\)*$//')
	prefix2=$(basename ${R2} | sed -e 's/.fastq\(.gz\)*$//')
	
	global_align "$sample_dir" "$R1" "$prefix1" "$UNMAP" "$FILT_UNMAP"&
	pid1=$!
	global_align "$sample_dir" "$R2" "$prefix2" "$UNMAP" "$FILT_UNMAP"&
	pid2=$!

	wait $pid1 $pid2 || die "Error in Bowtie alignment"
    done
    if [[ $wasrun != 1 ]]; then
	echo "Nothing to align ! Please check input files and R1/R2 extension." >&2
	exit 1
    fi
elif [[ ${MODE} == 'local' ]]; then
    if [[ ! -z $LIGATION_SITE ]]; then
	for r in $(get_fastq_for_bowtie_local)
	do
	    wasrun=1
	    R1=$r
	    R2=$(echo $r | get_R2)

	    echo "-----------"
	    echo $R1
	    echo $R2
	    echo "-------------"

	    sample_dir=$(get_sample_dir $r)
	    prefix1=$(basename ${R1} | sed -e 's/.fastq\(.gz\)*$//')
	    prefix2=$(basename ${R2} | sed -e 's/.fastq\(.gz\)*$//')
	    
	    local_align "$sample_dir" "$R1" "$prefix1" "$UNMAP"&
	    pid1=$!
	    local_align "$sample_dir" "$R2" "$prefix2" "$UNMAP"&
	    pid2=$!
	    
	    wait $pid1 $pid2 || die "Error in Bowtie alignment"
	done
	if [[ $wasrun != 1 ]]; then
	    echo "Nothing to align ! Please check input files and R1/R2 extension." >&2
	    exit 1
	fi
    fi
else
    die "Error: Unknown mapping mode !"
fi
