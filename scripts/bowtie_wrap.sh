#!/bin/bash
## Nicolas Servant 
## Eric Viara updated 2014-04-28
##

dir=$(dirname $0)

. $dir/hic.inc.sh

mode='global'

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-l) mode='local'; shift;;
	(-p) FASTQLIST=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

read_config $ncrna_conf

global_align()
{
    local sample_dir="$1"
    local file="$2"
    local prefix="$3"
    local cmd

    echo ${file} >> ${LOGFILE}

    mkdir -p ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}

    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${file} -S ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.sam 2>>${LOGS_DIR}/bowtie_${prefix}_global_${ORGANISM}.log"
	
    exec_cmd $cmd

    # Generate BAM files
    cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.sam > ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${ORGANISM}.bwt2glob.bam"

    exec_cmd $cmd
}

local_align()
{
    local sample_dir="$1"
    local file="$2"
    local prefix="$3"
    local cmd

    echo ${file} >> ${LOGFILE}

    mkdir -p ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}

    ## Align reads
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${file} -S ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.sam 2>>${LOGS_DIR}/bowtie_${prefix}_local.log"

    exec_cmd "$cmd"

    ## Generate BAM files
    cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.sam > ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.bam"

    exec_cmd "$cmd"
}

echo "BOWTIE_FASTQ_WRAP mode $mode"

if [[ ${mode} == 'global' ]]; then
    for r in $(get_fastq_for_bowtie_global)
    do
	R1=$r
	R2=$(echo $r | get_R2)
	sample_dir=$(get_sample_dir $r)
	prefix1=$(basename ${R1} | sed -e 's/.fastq//')
	prefix2=$(basename ${R2} | sed -e 's/.fastq//')
	
	global_align "$sample_dir" "$R1" "$prefix1" &
	global_align "$sample_dir" "$R2" "$prefix2" &

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
	
	local_align "$sample_dir" "$R1" "$prefix1" &
	local_align "$sample_dir" "$R2" "$prefix2" &

	wait
    done
fi
