#!/bin/bash
## Nicolas Servant 
## Eric Viara updated 2014-04-28
##

## Init
dir=$(dirname $0)
##. $dir/hic.inc.sh
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
##read_config $CONF
CONF=$conf_file . $dir/hic.inc.sh


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
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --rg-id BMG --rg SM:${prefix} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${file} 2> ${LDIR}/bowtie_${prefix}_global_${REFERENCE_GENOME}.log | ${SAMTOOLS_PATH}/samtools view -F 4 -bS - > ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.bam"
    exec_cmd $cmd

    # Generate BAM files with map reads only
    #cmd="${SAMTOOLS_PATH}/samtools view -F 4 -bS ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.sam > ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.bam 2>>${LDIR}/bowtie_${prefix}_global_${REFERENCE_GENOME}.log"
    #exec_cmd $cmd
}

## Local Alignment
local_align()
{
    local sample_dir="$1"
    local file="$2"
    local prefix="$3"
    local unmap="$4"
    local cmd

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
    tfile=`echo $file | sed -e s/.fastq/_trimmed.fastq/`
    if [[ ${RM_LOCAL_NO_CUTSITE} == 1 ]]; then
	${SCRIPTS}/cutsite_trimming --fastq $file --cutsite ${LIGATION_SITE} --out $tfile --rmuntrim > ${LDIR}/readsTrimming.log 2>&1
    else
	${SCRIPTS}/cutsite_trimming --fastq $file --cutsite ${LIGATION_SITE} --out $tfile  > ${LDIR}/readsTrimming.log 2>&1
    fi

    ## Unmapped reads
    if [[ $unmap == 1 ]]; then
	BOWTIE2_LOCAL_OPTIONS=${BOWTIE2_LOCAL_OPTIONS}" --un ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
    fi

    ## Run bowtie
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg SM:${prefix} --${FORMAT}-quals -p ${N_CPU} -x ${BOWTIE2_IDX} -U ${tfile} 2>${LDIR}/bowtie_${prefix}_local.log | ${SAMTOOLS_PATH}/samtools view -bS - > ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.bam"
    exec_cmd "$cmd"

    ## Generate BAM files with all reads so that the sum of global + local reads = total reads
    #cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.sam > ${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}/${prefix}_bwt2loc.bam 2>>${LDIR}/bowtie_${prefix}_local.log"
    #exec_cmd "$cmd"
}


echo -e "\nBOWTIE_FASTQ_WRAP mode $MODE\n"

if [[ ${MODE} == 'global' ]]; then
    for r in $(get_fastq_for_bowtie_global)
    do
       	R1=$r
	R2=$(echo $r | get_R2)
	
	sample_dir=$(get_sample_dir $r)
	prefix1=$(basename ${R1} | sed -e 's/.fastq\(.gz\)*//')
	prefix2=$(basename ${R2} | sed -e 's/.fastq\(.gz\)*//')
	
	global_align "$sample_dir" "$R1" "$prefix1" "$UNMAP"&
	pid1=$!
	global_align "$sample_dir" "$R2" "$prefix2" "$UNMAP"&
	pid2=$!

	wait $pid1 $pid2 || die "Error in Bowtie alignment"
    done
else
    for r in $(get_fastq_for_bowtie_local)
    do
	R1=$r
	R2=$(echo $r | get_R2)
	sample_dir=$(get_sample_dir $r)
	prefix1=$(basename ${R1} | sed -e 's/.fastq\(.gz\)*//')
	prefix2=$(basename ${R2} | sed -e 's/.fastq\(.gz\)*//')
	
	local_align "$sample_dir" "$R1" "$prefix1" "$UNMAP"&
	pid1=$!
	local_align "$sample_dir" "$R2" "$prefix2" "$UNMAP"&
	pid2=$!

	wait $pid1 $pid2 || die "Error in Bowtie alignment"
    done
fi
