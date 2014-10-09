#!/bin/bash
## Nicolas Servant 
## Eric Viara updated 2014-04-28
##

## Init
dir=$(dirname $0)
. $dir/hic.inc.sh
MODE='global'

## Get args
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) CONF=$2; shift;;
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
read_config $CONF

## Bowtie2 wrapper
## Global Alignment
global_align()
{
    local sample_dir="$1"
    local file="$2"
    local prefix="$3"
    local unmap="$4"
    local cmd
    echo ${file} >> ${LOGFILE}
    mkdir -p ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}

    ## Unmapped reads
    if [[ $unmap == 1 ]]; then
	BOWTIE2_GLOBAL_OPTIONS=${BOWTIE2_GLOBAL_OPTIONS}" --un ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.unmap.fastq"
    fi

    ## Run bowtie
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --rg-id ${prefix} --rg BM:G --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${file} -S ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.sam 2>>${LOGS_DIR}/bowtie_${prefix}_global_${ORGANISM}.log"
    exec_cmd $cmd

    # Generate BAM files with map reads only
    cmd="${SAMTOOLS_PATH}/samtools view -F 4 -bS ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.sam > ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.bam"

    exec_cmd $cmd
}

## Local Alignment
local_align()
{
    local sample_dir="$1"
    local file="$2"
    local prefix="$3"
    local unmap="$4"
    local cmd

    echo ${file} >> ${LOGFILE}
    mkdir -p ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}

    ## Unmapped reads
    if [[ $unmap == 1 ]]; then
	BOWTIE2_LOCAL_OPTIONS=${BOWTIE2_LOCAL_OPTIONS}" --un ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.unmap.fastq"
    fi

    ## Run bowtie
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --rg-id ${prefix} --rg BM:L  --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${file} -S ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.sam 2>>${LOGS_DIR}/bowtie_${prefix}_local.log"

    exec_cmd "$cmd"

    ## Generate BAM files with all reads so that the sum of global + local reads = total reads
    cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.sam > ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.bam"

    exec_cmd "$cmd"
}




echo "BOWTIE_FASTQ_WRAP mode $MODE"

if [[ ${MODE} == 'global' ]]; then
    for r in $(get_fastq_for_bowtie_global)
    do
	R1=$r
	R2=$(echo $r | get_R2)
	sample_dir=$(get_sample_dir $r)
	prefix1=$(basename ${R1} | sed -e 's/.fastq//')
	prefix2=$(basename ${R2} | sed -e 's/.fastq//')
	
	global_align "$sample_dir" "$R1" "$prefix1" "$UNMAP"&
	global_align "$sample_dir" "$R2" "$prefix2" "$UNMAP"&

	wait
    done
else
    for r in $(get_fastq_for_bowtie_local)
    do
	R1=$r
	R2=$(echo $r | get_R2)
	sample_dir=$(get_sample_dir $r)
	prefix1=$(basename ${R1} | sed -e 's/.fastq//')
	prefix2=$(basename ${R2} | sed -e 's/.fastq//')
	
	local_align "$sample_dir" "$R1" "$prefix1" "$UNMAP"&
	local_align "$sample_dir" "$R2" "$prefix2" "$UNMAP"&

	wait
    done
fi
