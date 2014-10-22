#!/bin/bash
## Nicolas Servant updated 2014-08-13
## Eric Viara updated 2014-04-28
##

dir=$(dirname $0)
##. $dir/hic.inc.sh

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
##read_config $CONF

##
## Combine Global and Local Bowtie2 mapping
##
mapping_combine()
{
    local sample_dir="$1"
    local file="$2"
    local prefix=$(echo ${sample_dir}/$(basename $file) | sed -e 's/.bwt2glob.sam//')

    echo ${prefix} >> ${LOGFILE}

    mkdir -p ${BOWTIE2_FINAL_OUTPUT_DIR}/${sample_dir}    

    ## Merge local and global alignment
    cmd="${SAMTOOLS_PATH}/samtools merge -n -f ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.bam ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.bam "
    exec_cmd $cmd

    ## Sort merge file. In theory, should be perform by "merge -n", but do not work in some cases ... depending on read name ?
    cmd="${SAMTOOLS_PATH}/samtools sort -n ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.sorted"
    exec_cmd $cmd
    cmd="mv ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.sorted.bam ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam"
    exec_cmd $cmd

    ## Generate SAM files with header
    cmd="${SAMTOOLS_PATH}/samtools view -h ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam > ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.sam"
    exec_cmd $cmd 
}


## Combine local and global alignments
for r in $(get_sam_for_combine)
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir $r)

    mapping_combine $sample_dir $R1 &
    mapping_combine $sample_dir $R2 &
    
    wait
done
