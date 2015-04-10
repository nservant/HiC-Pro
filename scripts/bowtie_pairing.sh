#!/bin/bash
## HiC-Pro
## Copyleft 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Launcher of mergeSAM script
## Merge R1 and R2 SAM files into one paired-end SAM file
##

dir=$(dirname $0)

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

CONF=$conf_file . $dir/hic.inc.sh

##
## Merge SE file into one PE file and filter the reads
##
merge_pairs()
{
    local sample_dir="$1"
    local file_r1="$2"
    local file_r2="$3"

    local prefix_r1=$(echo ${sample_dir}/$(basename $file_r1) | sed -e 's/.bwt2merged.bam//')
    local prefix_r2=$(echo ${sample_dir}/$(basename $file_r2) | sed -e 's/.bwt2merged.bam//')
    local prefix_out=$(echo $prefix_r1 | get_pairs)

    ## Merge two SAM files into 1 paired-end SAM file / removed unmapped and multihits reads
    OPTS="-q ${MIN_MAPQ} -t -v"
    if [[ ${RM_SINGLETON} == 0 ]]; then
	OPTS=$OPTS" -s"
    fi
    if [[ ${RM_MULTI} == 0 ]]; then
	OPTS=$OPTS" -m"
    fi

    ## Logs
    LDIR=${LOGS_DIR}/${sample_dir}
    mkdir -p ${LDIR}

    ## Index BAM
    #cmd="${SAMTOOLS_PATH}/samtools index ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r1}.bwt2merged.bam"
    #exec_cmd $cmd  
    #cmd="${SAMTOOLS_PATH}/samtools index ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r2}.bwt2merged.bam"
    #exec_cmd $cmd  

    cmd="${PYTHON_PATH}/python ${SCRIPTS}/mergeSAM.py ${OPTS} -f ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r1}.bwt2merged.bam -r ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r2}.bwt2merged.bam -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.bam > ${LDIR}/mergeSAM.log"
    exec_cmd $cmd

    ## Generate BAM file
    ##cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.sam > ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.bam"
    ##exec_cmd $cmd

    ## Generate index file
    ## cmd="${SAMTOOLS_PATH}/samtools index ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.bam"
    ## exec_cmd $cmd
}

## Combine R1/R2 tags in a single BAM file
for r in $(get_sam_for_merge)
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir $r)

    merge_pairs $sample_dir $R1 $R2 
done
